import contextlib
import enum
import glob
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Generator, List, Optional, Tuple, Union

import jinja2
import pandas as pd

PathLike = Union[str, Path]


SCHRODINGER_PATH = os.getenv("SCHRODINGER_PATH")
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_ENV = jinja2.Environment(
    loader=jinja2.FileSystemLoader(os.path.join(SCRIPT_DIR, "templates"))
)

if not SCHRODINGER_PATH:
    print("error: Environment variable SCHRODINGER_PATH is not set", file=sys.stderr)
    sys.exit(1)


class RankingCriterion(enum.Enum):
    CONTACTS = 1
    NORMALIZED_CONTACTS = 2
    BINDING_ENERGY = 3
    NORMALIZED_BINDING_ENERGY = 4
    TOTAL_SCORE = 5

    def ascending(self) -> bool:
        return self == self.BINDING_ENERGY


RANKING_COLUMN_MAP = {
    RankingCriterion.CONTACTS: "INT",
    RankingCriterion.NORMALIZED_CONTACTS: "INT_NORM",
    RankingCriterion.BINDING_ENERGY: "DGBIND",
    RankingCriterion.NORMALIZED_BINDING_ENERGY: "DGBIND_NORM",
    RankingCriterion.TOTAL_SCORE: "NORMT",
}


def normalize(series: pd.Series) -> pd.Series:
    min_value = series.min()
    return (series - min_value) / (series.max() - min_value)


def rank_poses(
    df: pd.DataFrame,
    criteria: List[RankingCriterion] = [
        RankingCriterion.NORMALIZED_CONTACTS,
        RankingCriterion.TOTAL_SCORE,
    ],
) -> pd.DataFrame:
    df = df.copy(deep=True)
    df["INT_NORM"] = normalize(df["INT"])
    df["DGBIND_NORM"] = 1 - normalize(df["DGBIND"])
    df["NORMT"] = (df["INT_NORM"] + df["DGBIND_NORM"]) / 2

    sort_fields = [RANKING_COLUMN_MAP[criterion] for criterion in criteria]
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


def rank_poses_cross(
    results: pd.DataFrame,
    criteria: List[RankingCriterion] = [
        RankingCriterion.NORMALIZED_CONTACTS,
        RankingCriterion.TOTAL_SCORE,
    ],
) -> pd.DataFrame:
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
    return pd.DataFrame(rows).sort_values("GLOBAL_NORMT", ascending=True)


def report(
    output_dir: PathLike, resids: str, contact_cutoff: float, use_existing: bool = True
) -> pd.DataFrame:
    csvfile = os.path.join(output_dir, f"{os.path.basename(output_dir)}.csv")
    if not use_existing or not os.path.exists(csvfile):
        run_silent(
            os.path.join(SCHRODINGER_PATH, "run.exe"),
            os.path.join(SCRIPT_DIR, "scripts", "report.py"),
            str(output_dir),
            cutoff=contact_cutoff,
            output=csvfile,
            residues=resids,
        )
    return pd.read_csv(csvfile)


def report_cross(
    output_dir: PathLike,
    bs_residues: Dict[str, str],
    contact_cutoff: float,
    use_existing: bool = True,
) -> pd.DataFrame:
    results: List[pd.DataFrame] = []
    for prot_dir in glob.glob(os.path.join(output_dir, "*")):
        if not os.path.isdir(prot_dir):
            continue
        prot_name = os.path.basename(prot_dir)
        df = report(
            output_dir=prot_dir,
            resids=bs_residues[prot_name],
            contact_cutoff=contact_cutoff,
            use_existing=use_existing,
        )
        df.insert(0, "PROTEIN", [prot_name for _ in range(len(df))])
        results.append(df)
    return pd.concat(results).reset_index(drop=True)


def run_silent(
    program: str,
    *args: str,
    use_single_dash: bool = False,
    **kwargs: Union[bool, int, float, str],
) -> None:
    cmd = [program] + list(args)
    for option, value in kwargs.items():
        prefix = "-" if use_single_dash else ("--" if len(option) > 1 else "-")
        option = prefix + option.replace("_", "-")

        cmd.append(option)
        if not isinstance(value, bool):
            cmd.append(str(value))
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)


def run_ifd(
    prot_file: PathLike,
    lig_file: PathLike,
    resids: str,
    glide_cpus: int = 2,
    prime_cpus: int = 2,
) -> None:
    inp_file = f"{Path(lig_file).stem}.inp"
    with open(inp_file, "w") as io:
        template = TEMPLATE_ENV.get_template("ifd.inp")
        vars = dict(
            protfile=prot_file,
            ligfile=Path(lig_file).name,
            resids=resids,
        )
        io.write(template.render(vars))

    run_silent(
        os.path.join(SCHRODINGER_PATH, "ifd"),
        inp_file,
        NGLIDECPU=glide_cpus,
        NPRIMECPU=prime_cpus,
        HOST="localhost",
        SUBHOST="localhost",
        TMPLAUNCHDIR=True,
        WAIT=True,
        use_single_dash=True,
    )


def run_ifd_cross(
    prot_files: List[PathLike],
    lig_files: List[PathLike],
    bs_residues: Dict[str, str],
    workdir: PathLike = "ifd",
    glide_cpus: int = 2,
    prime_cpus: int = 2,
) -> None:
    for prot_file in map(Path, prot_files):
        prot_name = prot_file.stem

        prot_workdir = Path(workdir, prot_name)
        prot_workdir.mkdir(parents=True, exist_ok=True)
        shutil.copy(prot_file, prot_workdir)

        for i, lig_file in enumerate(map(Path, lig_files)):
            lig_name = lig_file.stem
            jobid = f"{prot_name}/{lig_name}"

            lig_workdir = prot_workdir / lig_name
            if lig_workdir.exists():
                print(f"Skipping IFD for {jobid}... already exists")
                continue
            else:
                lig_workdir.mkdir()
            shutil.copy(lig_file, lig_workdir)

            print(f"Running IFD for {jobid} [{i}/{len(lig_files)}]...")
            with transient_dir(lig_workdir):
                run_ifd(
                    os.path.join("..", prot_file.name),
                    lig_file.name,
                    bs_residues[prot_file.name],
                    glide_cpus,
                    prime_cpus,
                )


def run_mmgbsa(ifd_file: PathLike, cpus: int = 2) -> None:
    run_silent(
        os.path.join(SCHRODINGER_PATH, "prime_mmgbsa"),
        str(ifd_file),
        ligand="(res.pt UNK)",
        job_type="REAL_MIN",
        out_type="COMPLEX",
        csv_output="yes",
        JOBNAME=Path(ifd_file).stem.replace("-out", ""),
        HOST=f"localhost:{cpus}",
        WAIT=True,
        use_single_dash=True,
    )


def run_mmgbsa_cross(
    ifd_files: List[PathLike],
    workdir: PathLike,
    cpus: int = 2,
) -> None:
    for i, ifd_file in enumerate(map(Path, ifd_files)):
        ifd_file = ifd_file.absolute()
        prot_name = ifd_file.parent.parts[-1]

        prot_workdir = Path(workdir, prot_name)
        prot_workdir.mkdir(parents=True, exist_ok=True)

        lig_name = ifd_file.stem.replace("-out", "")
        jobid = f"{prot_name}/{lig_name}"

        if os.path.exists(prot_workdir / f"{lig_name}-out.maegz"):
            print(f"Skipping MM/GBSA for {jobid}... already exists")
            continue

        print(f"Running MM/GBSA for {jobid} [{i}/{len(ifd_files)}]...")
        with transient_dir(prot_workdir):
            run_mmgbsa(ifd_file, cpus)


@contextlib.contextmanager
def transient_dir(path: PathLike) -> Generator[None, None, None]:
    cwd = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(cwd)


def analysis(
    workdir: str,
    bs_residues: Dict[str, str],
    radius: float,
    rank_criteria: List[RankingCriterion] = [
        RankingCriterion.NORMALIZED_CONTACTS,
        RankingCriterion.TOTAL_SCORE,
    ],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    results = report_cross(workdir, bs_residues, radius)
    cross_results = rank_poses_cross(results, rank_criteria)
    return results, cross_results
