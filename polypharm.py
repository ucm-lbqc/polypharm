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

pd.options.display.max_rows = 100

PathLike = Union[str, Path]


SCHRODINGER_PATH = os.getenv("SCHRODINGER_PATH")
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_ENV = jinja2.Environment(
    loader=jinja2.FileSystemLoader(os.path.join(SCRIPT_DIR, "templates"))
)

if not SCHRODINGER_PATH:
    print("error: Environment variable SCHRODINGER_PATH is not set", file=sys.stderr)
    exit(1)


class RankingCriterion(enum.Enum):
    CONTACTS = 1
    NORMALIZED_CONTACTS = 2
    BINDING_ENERGY = 3
    NORMALIZED_BINDING_ENERGY = 4
    TOTAL_SCORE = 5

    def ascending(self) -> bool:
        return self in [
            self.NORMALIZED_CONTACTS,
            self.NORMALIZED_BINDING_ENERGY,
            self.BINDING_ENERGY,
            self.TOTAL_SCORE,
        ]


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
    limit: Optional[int] = None,
) -> pd.DataFrame:
    df = df.copy(deep=True)
    df["INT_NORM"] = normalize(df["INT"])
    df["DGBIND_NORM"] = 1 - normalize(df["DGBIND"])
    df["NORMT"] = df["INT_NORM"] + df["DGBIND_NORM"]

    df.sort_values(
        by=[RANKING_COLUMN_MAP[criterion] for criterion in criteria],
        ascending=[criterion.ascending() for criterion in criteria],
        inplace=True,
    )
    df.reset_index(drop=True, inplace=True)
    df["RANK"] = list(range(len(df)))

    if limit:
        top_molecules = df["NAME"].unique()[:limit]
        df = df[df["NAME"].isin(top_molecules)]

    return df


def rank_poses_cross(
    results: pd.DataFrame,
    criteria: List[RankingCriterion] = [
        RankingCriterion.NORMALIZED_CONTACTS,
        RankingCriterion.TOTAL_SCORE,
    ],
    limit: Optional[int] = None,
) -> pd.DataFrame:
    results = pd.concat(
        rank_poses(df, criteria, limit) for _, df in results.groupby("PROTEIN")
    )

    common_names: List[str] = set.intersection(
        *[set(names) for names in results.groupby("PROTEIN")["NAME"].unique()]
    )
    results = results[results["NAME"].isin(common_names)]

    rows = []
    for [name, index], df in results.groupby(["NAME", "INDEX"]):
        row_data = dict(NAME=name, INDEX=index)
        for row in df.iterrows():
            prot = row["PROTEIN"]
            row_data[f"{prot}_INT"] = row["INT"]
            row_data[f"{prot}_DGBIND"] = row["DGBIND"]
            row_data[f"{prot}_NORMT"] = row["NORMT"]
            row_data[f"{prot}_RANK"] = row["RANK"]
        row_data["SUM_RANK"] = df["RANK"].sum()
        row_data["SUM_NORMT"] = df["NORMT"].sum()
        rows.append(row_data)
    return pd.DataFrame(rows).sort_values("SUM_NORMT", ascending=True)


def report(output_dir: PathLike, resids: str, contact_cutoff: float) -> pd.DataFrame:
    csvfile = os.path.join(output_dir, f"{os.path.basename(output_dir)}.csv")
    run_silent(
        os.path.join(SCHRODINGER_PATH, "run"),
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
) -> pd.DataFrame:
    results: List[pd.DataFrame] = []
    for prot_dir in glob.glob(os.path.join(output_dir, "*")):
        prot_name = os.path.basename(prot_dir)
        df = report(
            output_dir=prot_dir,
            resids=bs_residues[f"{prot_name}.mae"],
            contact_cutoff=contact_cutoff,
        )
        df["PROTEIN"] = prot_name
        results.append(df)
    return pd.concat(results)


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
            cmd.append(repr(value))

    with open(os.devnull, "w") as io:
        try:
            subprocess.run(cmd, check=True, stdout=io, stderr=io)
        except subprocess.CalledProcessError:
            print("Error occured executing {}".format(" ".join(cmd)), file=sys.stderr)


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
