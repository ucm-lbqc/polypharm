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


def rank_poses(
    df: pd.DataFrame,
    # FIXME: Use total score as default
    criteria: List[RankingCriterion] = [
        RankingCriterion.NORMALIZED_CONTACTS,
        RankingCriterion.TOTAL_SCORE,
    ],
) -> pd.DataFrame:
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
    # FIXME: Use total score as default
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
    maefiles: List[PathLike],
    bs_residues: Dict[str, List[str]],
    contact_cutoff: float,
    use_existing: bool = True,
    tasks: int = 1,
    quiet: bool = False,
) -> pd.DataFrame:
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


def run_ifd_cross(
    prot_files: List[PathLike],
    lig_files: List[PathLike],
    bs_residues: Dict[str, List[str]],
    workdir: PathLike = "ifd",
    glide_cpus: int = 2,
    prime_cpus: int = 2,
    tasks: int = 1,
    quiet: bool = False,
) -> None:
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

            inp_file = f"{lig_file.stem}.inp"
            render_template_to_file(
                "ifd.inp",
                inp_file,
                protfile=os.path.join("..", prot_file.name),
                ligfile=lig_file.name,
                resids=",".join(bs_residues[prot_file.name]),
            )

            args = [
                os.path.join(get_schrodinger_path(), "ifd"),
                inp_file,
                "-NGLIDECPU",
                str(glide_cpus),
                "-NPRIMECPU",
                str(prime_cpus),
                "-HOST",
                "localhost",
                "-SUBHOST",
                "localhost",
                "-TMPLAUNCHDIR",
                "-WAIT",
            ]
            commands.append(Command(jobid, args, workdir=lig_workdir))
    async_run(concurrent_subprocess(commands, tasks, quiet))


def run_mmgbsa_cross(
    ifd_files: List[PathLike],
    workdir: PathLike = "mmgbsa",
    cpus: int = 2,
    tasks: int = 1,
    quiet: bool = False,
) -> None:
    commands: List[Command] = []
    for ifd_file in map(Path, ifd_files):
        ifd_file = ifd_file.absolute()
        prot_name = ifd_file.parent.parts[-1]
        prot_workdir = Path(workdir, prot_name)

        lig_name = ifd_file.stem.replace("-out", "")
        jobid = f"mmgbsa/{prot_name}/{lig_name}"

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
