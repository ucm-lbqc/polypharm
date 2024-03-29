import enum
import glob
import os
import shutil
from pathlib import Path
from typing import Dict, List

import pandas as pd

from .concurrency import Command, async_run, concurrent_subprocess
from .helpers import (
    PathLike,
    get_schrodinger_path,
    get_script_path,
    render_template_to_file,
)


class RankingCriterion(enum.Enum):
    """Criterion for docking pose ranking."""

    BINDING_ENERGY = "DGBIND"
    """Binding free energy. Usually computed via MM/GBSA for a
    protein-ligand complex. Lower is better."""
    CONTACTS = "INT"
    """Number of residue contacts. Encodes possible interactions between
    the protein's residue and ligand atoms. Larger is better."""
    NORMALIZED_BINDING_ENERGY = "DGBIND_NORM"
    """Normalized binding free energy (values within 0 and 1) computed
    from energies of multiple docking poses. Larger is better."""
    NORMALIZED_CONTACTS = "INT_NORM"
    """Normalized residue contacts (values within 0 and 1) computed from
    number of contacts of multiple docking poses. Larger is better."""
    TOTAL_SCORE = "NORMT"
    """Total score is defined as (NORMALIZED_BINDING_ENERGY +
    NORMALIZED_CONTACTS) / 2 (values within 0 and 1). Larger is
    better."""

    def ascending(self) -> bool:
        """Returns whether the criterion is in ascending order."""
        return self == self.BINDING_ENERGY


def rank_poses(
    df: pd.DataFrame,
    # FIXME: Use total score as default
    criteria: List[RankingCriterion] = [
        RankingCriterion.NORMALIZED_CONTACTS,
        RankingCriterion.TOTAL_SCORE,
    ],
) -> pd.DataFrame:
    """Rank docking poses by the given criteria.

    Columns for the normalized energies and number of contacts
    (``DGBIND_NORM`` and ``INT_NORM``), the total score (``NORMT``), and
    docking ranking (``RANK``) per protein are computed.

    :param pd.DataFrame df: Table containing the protein and ligand
        names, energies, and residue contacts (generated by
        :py:func:`report`).
    :param List[RankingCriterion] criteria: List of ranking criteria
        (See :py:class:`RankingCriterion`). Defaults to
        ``NORMALIZED_CONTACTS`` and ``TOTAL_SCORE``.
    :return pd.DataFrame: Table with the entries sorted by the given
        criteria.
    """

    def normalize(series: pd.Series) -> pd.Series:
        min_value = series.min()
        return (series - min_value) / (series.max() - min_value)

    df = df.copy(deep=True)

    if "PROTEIN" in df.columns:
        for _, sub_df in df.groupby("PROTEIN"):
            df.loc[sub_df.index, "INT_NORM"] = normalize(sub_df["INT"])
            df.loc[sub_df.index, "DGBIND_NORM"] = 1 - normalize(sub_df["DGBIND"])
    else:
        df["INT_NORM"] = normalize(df["INT"])
        df["DGBIND_NORM"] = 1 - normalize(df["DGBIND"])
    df["NORMT"] = (df["INT_NORM"] + df["DGBIND_NORM"]) / 2

    sort_fields = [criterion.value for criterion in criteria]
    ascending = [criterion.ascending() for criterion in criteria]
    if "PROTEIN" in df.columns:
        sort_fields = ["PROTEIN"] + sort_fields
        ascending = [True] + ascending
    df.sort_values(
        by=sort_fields,
        ascending=ascending,
        inplace=True,
    )
    df.reset_index(drop=True, inplace=True)

    if "PROTEIN" in df.columns:
        for _, sub_df in df.groupby("PROTEIN"):
            df.loc[sub_df.index, "RANK"] = list(range(len(sub_df)))
    else:
        df["RANK"] = list(range(len(df)))

    return df


def rank_molecules(
    results: pd.DataFrame,
    # FIXME: Use total score as default
    criteria: List[RankingCriterion] = [
        RankingCriterion.NORMALIZED_CONTACTS,
        RankingCriterion.TOTAL_SCORE,
    ],
) -> pd.DataFrame:
    """Rank molecules by the given criteria across multiple proteins.

    The molecules are ranked according to the docking ranking (see
    :py:func:`rank_poses`) obtained in multiple proteins, such that the
    top-ranked molecules are those that can potentially bind to multiple
    receptors determined by either the binding energies or residues
    contacts, or both. Therefore, molecule without poses for one or more
    proteins are excluded from the analysis.

    Global scores and ranking are derived from the docking scores (see
    :py:class:`RankingCriterion.TOTAL_SCORE`) and ranking per protein.

    The following columns are dinamically generated for each protein:
    ``'<protein>_INT'``, ``'<protein>_DGBIND'``, ``'<protein>_NORMT'``,
    and ``'<protein>_RANK'``.

    :param pd.DataFrame df: Table containing the protein and ligand
        names, energies, residue contacts, and docking ranking
        (generated by :py:func:`rank_poses`).
    :param List[RankingCriterion] criteria: List of ranking criteria
        (See :py:class:`RankingCriterion`). Defaults to
        ``NORMALIZED_CONTACTS`` and ``TOTAL_SCORE``.
    :return pd.DataFrame: Table containing the molecule name, energies,
        and total contacts per protein, and global scores and ranking
        for each molecule.
    """
    results = rank_poses(results, criteria)

    common_names: List[str] = set.intersection(
        *[set(names) for names in results.groupby("PROTEIN")["NAME"].unique()]
    )
    results = results[results["NAME"].isin(common_names)]
    n_proteins = results["PROTEIN"].nunique()

    rows = []
    for name, df in results.groupby("NAME"):
        if len(df) < n_proteins:  # missing poses for some proteins
            continue
        row_data = dict(NAME=name)
        rank_sum = 0
        total_score_sum = 0
        # TODO: test other poses besides top-ranked per protein
        for row in [sub_df.iloc[0] for _, sub_df in df.groupby("PROTEIN")]:
            prot = row["PROTEIN"]
            row_data[f"{prot}_INDEX"] = row["INDEX"]
            row_data[f"{prot}_INT"] = row["INT"]
            row_data[f"{prot}_DGBIND"] = row["DGBIND"]
            row_data[f"{prot}_NORMT"] = row["NORMT"]
            row_data[f"{prot}_RANK"] = row["RANK"]
            rank_sum += row["RANK"]
            total_score_sum += row["NORMT"]
        # FIXME: use averages, not sums
        row_data["GLOBAL_RANK"] = rank_sum
        row_data["GLOBAL_NORMT"] = total_score_sum
        rows.append(row_data)
    results = (
        pd.DataFrame(rows)
        .sort_values("GLOBAL_RANK", ascending=True)
        .reset_index(drop=True)
    )
    results.index +=1
    results.index.name = "POSITION"
    sorted_columns = sorted(results.columns[1:-2], key=lambda x: (x.split("_")[-1], x))
    sorted_columns = ["NAME"] + sorted_columns + ["GLOBAL_RANK", "GLOBAL_NORMT"]
    return results[sorted_columns]


def report(
    maefiles: List[PathLike],
    bs_residues: Dict[str, List[str]],
    contact_cutoff: float = 5,
    use_existing: bool = True,
    tasks: int = 1,
    quiet: bool = False,
) -> pd.DataFrame:
    """Generate a report for the protein-ligand complexes.

    Each protein-ligand complex is inspected for noncovalent contacts
    (controlled by *contact_cutoff*) between the protein and ligand
    structures. The contacts are checked for the given binding site
    residues (*bs_residues*) only. Computed energies by IFD and MM/GBSA,
    and found contacts are encoded as a :py:class:`pd.DataFrame`
    instance, where each row corresponds to a protein/ligand pair, and
    columns correspond to the protein and ligand name, energies, and
    inspected residues. Contacts are reported as either 1 or 0 whether
    they were found or not, respectively. Contact columns are
    dinamically generated for each given residue.

    Inspection of protein-ligand complexes is run as multiple
    subprocesses by invoking the ``$SCHRODINGER_PATH/run`` program to
    execute the ``report.py`` script, which uses the `Schrödinger's
    Python API <https://www.schrodinger.com/pythonapi>`_. The script
    generates a comma-separated value (CSV) file per Maestro file, which
    is read by this function. The CSV files can be re-used for future
    calls (e.g., previous molecules are not re-analyzed after adding new
    ones) unless *use_existing* is False.

    The subprocesses are run as parallel tasks (controlled by the
    *tasks* argument) via the
    :py:func:`concurrency.concurrent_subprocess` function.

    :param List[str | Path] maefiles: List of Maestro files
        (\\*.mae[gz]) containing the protein-ligand complexes (e.g.,
        obtained from MM/GBSA runs).
    :param Dict[str, List[str]] bs_residues: Dictionary of ``'protein':
        list-of-residues``, where the list of the residues (e.g.,
        'A:81') comprises the binding site. Note that the keys
        (proteins) must be equal to the basename (filename without
        extension) of the protein files. For instance, ``'protA'`` for a
        file named ``/path/to/protA.mae``.
    :param float contact_cutoff: Distance cutoff in Å to detect residue
        contacts. Defaults to 5 Å.
    :param bool use_existing: If True, read any existing individual
        output CSV file from previous run, else perform the report again
        overwriting existing files. Defaults to True.
    :param int tasks: Number of parallel tasks to run. Defaults to 1.
    :param bool quiet: Do not print progress to standard output.
        Defaults to False.
    :return pd.DataFrame: Table containing the protein and ligand names,
        energies, and residue contacts.
    """
    commands: List[Command] = []
    for maefile in maefiles:
        path = Path(maefile)
        prot_name = Path(maefile).parent.name
        lig_name = path.stem.replace("-out", "")
        csvfile = path.parent / f"{lig_name}-report.csv"
        args = [
            os.path.join(get_schrodinger_path(), "run"),
            get_script_path("report.py"),
            str(maefile),
            "--cutoff",
            str(contact_cutoff),
            "--output",
            str(csvfile),
            "--residues",
            ",".join(bs_residues[prot_name]),
        ]
        data = dict(csvfile=csvfile, maefile=maefile, prot_name=prot_name)
        jobid = f"report/{prot_name}/{lig_name}"
        commands.append(Command(jobid, args, data=data))

    commands_to_run = [
        cmd
        for cmd in commands
        if not use_existing or not os.path.exists(cmd.data["csvfile"])
    ]
    if commands_to_run:
        async_run(concurrent_subprocess(commands_to_run, tasks, quiet))

    results: List[pd.DataFrame] = []
    for cmd in commands:
        if not os.path.exists(cmd.data["csvfile"]):  # may have failed
            continue
        df = pd.read_csv(cmd.data["csvfile"])
        df.insert(0, "PROTEIN", [cmd.data["prot_name"] for _ in range(len(df))])
        results.append(df)
    return pd.concat(results).reset_index(drop=True)


def cross_dock(
    prot_files: List[PathLike],
    lig_files: List[PathLike],
    bs_residues: Dict[str, List[str]],
    workdir: PathLike = "ifd",
    glide_cpus: int = 1,
    prime_cpus: int = 1,
    tasks: int = 1,
    quiet: bool = False,
) -> None:
    """Run cross-docking using the induced-fit protocol.

    The induced-fit docking (IFD) protocol developed by Schrödinger_,
    based on Glide_ and the Refinement module in Prime_, accurately
    predicts ligand binding modes and concomitant structural changes in
    the receptor.

    IFD is run for each molecule in every protein, effectively
    performing a cross-docking. The calculations are run within the
    given working directory, where each IFD is exectured at the
    ``<protein>/<molecule>`` subdirectory (created on the fly).

    IFD calculations are executed as multiple subprocesses by invoking
    the ``$SCHRODINGER_PATH/ifd`` program. The subprocesses are run as
    parallel tasks (controlled by the *tasks* argument) via the
    :py:func:`concurrency.concurrent_subprocess` function.

    :param List[str | Path] prot_files: List of Maestro files
        (\\*.mae[gz]) containing the prepared protein structures.
    :param List[str | Path] lig_files: List of Maestro files
        (\\*.mae[gz]) containing the prepared structures of the
        molecules to be docked.
    :param Dict[str, List[str]] bs_residues: Dictionary of ``'protein':
        list-of-residues``, where the list of the residues (e.g.,
        'A:81') comprises the binding site. Note that the keys
        (proteins) must be equal to the basename (filename without
        extension) of the protein files. For instance, ``'protA'`` for a
        file named ``/path/to/protA.mae``.
    :param str | Path workdir: Working directory to run the IFD
        calculations. It will be created if missing. Defaults to "ifd".
    :param int glide_cpus: Number of processors to be used by Glide.
        Defaults to 1.
    :param int prime_cpus: Number of processors to be used by Prime.
        Defaults to 1.
    :param int tasks: Number of parallel tasks to run. Defaults to 1.
    :param bool quiet: Do not print progress to standard output.
        Defaults to False.

    .. note::
        Requires a Schrodinger Suite working installation (2018-4 or
        greater) including the Glide and Prime modules. The
        ``SCHRODINGER_PATH`` environment variable must be set.

    .. _Glide: https://schrodinger.com/products/glide
    .. _Prime: https://schrodinger.com/products/prime
    .. _Schrödinger: https://schrodinger.com
    """
    commands: List[Command] = []
    for prot_file in map(Path, prot_files):
        prot_name = prot_file.stem

        prot_workdir = Path(workdir, prot_name)
        prot_workdir.mkdir(parents=True, exist_ok=True)
        shutil.copy(prot_file, prot_workdir)

        for lig_file in map(Path, lig_files):
            lig_name = lig_file.stem
            jobid = f"ifd/{prot_name}/{lig_name}"

            lig_workdir = prot_workdir / lig_name
            if lig_workdir.exists():
                print(f"Skipping {jobid}... already exists")
                continue
            else:
                lig_workdir.mkdir()
            shutil.copy(lig_file, lig_workdir)

            inp_file = f"{lig_workdir}/{lig_file.stem}.inp"
            render_template_to_file(
                "ifd.inp",
                inp_file,
                protfile=os.path.join("..", prot_file.name),
                ligfile=lig_file.name,
                resids=",".join(bs_residues[prot_name]),
            )

            args = [
                os.path.join(get_schrodinger_path(), "ifd"),
                inp_file,
                "-NGLIDECPU",
                str(glide_cpus),
                "-NPRIMECPU",
                str(prime_cpus),
                "-REMOVE_WORKDIR",
                "-HOST",
                "localhost",
                "-SUBHOST",
                "localhost",
                "-TMPLAUNCHDIR",
                "-WAIT",
            ]
            commands.append(Command(jobid, args, workdir=lig_workdir))
    async_run(concurrent_subprocess(commands, tasks, quiet))


def rescore_docking(
    ifd_files: List[PathLike],
    workdir: PathLike = "mmgbsa",
    cpus: int = 1,
    tasks: int = 1,
    quiet: bool = False,
) -> None:
    """Run MM/GBSA to rescore IFD docking poses.

    The MM/GBSA method computes the binding free energy of
    protein-ligand complexes. The MM/GBSA implementation distributed by
    Schrödinger_ within Prime_ is used here, which employs the OPLS_2005
    force field [banks2005]_ and VSGB 2.0 solvation model [li2011]_.
    Note that the complex structures are energy minimized before
    computing the binding energy.

    The calculations are run within the given working directory, where
    each MM/GBSA is exectured at the ``<protein>`` subdirectory (created
    on the fly).

    MM/GBSA calculations are executed as multiple subprocesses by
    invoking the ``$SCHRODINGER_PATH/prime_mmgbsa`` program. The
    subprocesses are run as parallel tasks (controlled by the *tasks*
    argument) via the :py:func:`concurrency.concurrent_subprocess`
    function.

    :param List[str | Path] ifd_files: List of Maestro files
        (\\*.mae[gz]) containing the protein-ligand complexes (e.g.,
        obtained from IFD runs).
    :param str | Path workdir: Working directory to run the MM/GBSA
        calculations. It will be created if missing. Defaults to
        "mmgbsa".
    :param int cpus: Number of processors to be used by Prime. Defaults
        to 1.
    :param int tasks: Number of parallel tasks to run. Defaults to 1.
    :param bool quiet: Do not print progress to standard output.
        Defaults to False.

    .. note::
        Requires a Schrodinger Suite working installation (2018-4 or
        greater) including the Prime module. The ``SCHRODINGER_PATH``
        environment variable must be set.

    .. [banks2005] Banks, J. L., Beard, H. S., Cao, Y., Cho, A. E.,
        Damm, W., Farid, R., ... & Levy, R. M. (2005). Integrated
        modeling program, applied chemical theory (IMPACT). Journal of
        computational chemistry, 26(16), 1752-1780.
        https://doi.org/10.1002/jcc.20292
    .. [li2011] Li, J., Abel, R., Zhu, K., Cao, Y., Zhao, S., &
        Friesner, R. A. (2011). The VSGB 2.0 model: a next generation
        energy model for high resolution protein structure modeling.
        Proteins: Structure, Function, and Bioinformatics, 79(10),
        2794-2812. https://dx.doi.org/10.1002/prot.23106

    .. _Prime: https://schrodinger.com/products/prime
    .. _Schrödinger: https://schrodinger.com
    """
    commands: List[Command] = []
    for ifd_file in map(Path, ifd_files):
        ifd_file = ifd_file.absolute()
        prot_name = ifd_file.parent.parts[-2]
        prot_workdir = Path(workdir, prot_name)

        lig_name = ifd_file.stem.replace("-out", "")
        jobid = f"{workdir.name}/{prot_name}/{lig_name}"

        if os.path.exists(prot_workdir / f"{lig_name}-out.maegz"):
            print(f"Skipping {jobid}... already exists")
            continue

        args = [
            os.path.join(get_schrodinger_path(), "prime_mmgbsa"),
            str(ifd_file),
            "-ligand",
            "(res.pt UNK)",
            "-job_type",
            "REAL_MIN",
            "-out_type",
            "COMPLEX",
            "-csv_output",
            "yes",
            "-JOBNAME",
            Path(ifd_file).stem.replace("-out", ""),
            "-HOST",
            f"localhost:{cpus}",
            "-WAIT",
        ]
        commands.append(Command(jobid, args, workdir=prot_workdir))
    async_run(concurrent_subprocess(commands, tasks, quiet))
