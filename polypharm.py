import asyncio
import contextlib
import dataclasses
import enum
import glob
import os
import shutil
import subprocess
import sys
import threading
from pathlib import Path
from typing import Any, Coroutine, Dict, Generator, List, Optional, Tuple, Union

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
    output_dir: PathLike,
    bs_residues: Dict[str, str],
    contact_cutoff: float,
    use_existing: bool = True,
    tasks: int = 1,
) -> pd.DataFrame:
    commands: List[_Command] = []
    for maefile in glob.glob(os.path.join(output_dir, "**", "*-out.maegz")):
        csvfile = Path(maefile)
        csvfile = csvfile.parent / (csvfile.stem.replace("-out", "") + "-report.csv")
        prot_name = Path(maefile).parent.name
        args = [
            os.path.join(SCHRODINGER_PATH, "run"),
            os.path.join(SCRIPT_DIR, "scripts", "report.py"),
            str(maefile),
            "--cutoff",
            str(contact_cutoff),
            "--output",
            str(csvfile),
            "--residues",
            bs_residues[prot_name],
        ]
        data = dict(csvfile=csvfile, maefile=maefile, prot_name=prot_name)
        commands.append(_Command(args, data=data))

    commands_to_run = [
        cmd
        for cmd in commands
        if not use_existing or not os.path.exists(cmd.data["csvfile"])
    ]
    _async_run(_concurrent_subprocess(commands_to_run, tasks))

    results: List[pd.DataFrame] = []
    for cmd in commands:
        df = pd.read_csv(cmd.data["csvfile"])
        df.insert(0, "PROTEIN", [cmd.data["prot_name"] for _ in range(len(df))])
        results.append(df)
    return pd.concat(results).reset_index(drop=True)


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

    _run_silent(
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
    tasks: int = 1,
) -> None:
    commands: List[_Command] = []
    for prot_file in map(Path, prot_files):
        prot_name = prot_file.stem

        prot_workdir = Path(workdir, prot_name)
        prot_workdir.mkdir(parents=True, exist_ok=True)
        shutil.copy(prot_file, prot_workdir)

        for lig_file in map(Path, lig_files):
            lig_name = lig_file.stem
            jobid = f"{prot_name}/{lig_name}"

            lig_workdir = prot_workdir / lig_name
            if lig_workdir.exists():
                print(f"Skipping IFD for {jobid}... already exists")
                continue
            else:
                lig_workdir.mkdir()
            shutil.copy(lig_file, lig_workdir)

            inp_file = f"{lig_file.stem}.inp"
            with open(inp_file, "w") as io:
                template = TEMPLATE_ENV.get_template("ifd.inp")
                vars = dict(
                    protfile=os.path.join("..", prot_file.name),
                    ligfile=lig_file.name,
                    resids=bs_residues[prot_file.name],
                )
                io.write(template.render(vars))

            args = [
                os.path.join(SCHRODINGER_PATH, "ifd"),
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
            commands.append(_Command(args, workdir=lig_workdir))
    _async_run(_concurrent_subprocess(commands, tasks))


def run_mmgbsa(ifd_file: PathLike, cpus: int = 2) -> None:
    _run_silent(
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


# TODO: use _concurrent_subprocess
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
        with _transient_dir(prot_workdir):
            run_mmgbsa(ifd_file, cpus)


@dataclasses.dataclass
class _Command:
    args: List[str]
    workdir: Optional[PathLike] = None
    data: Dict[str, Any] = dataclasses.field(default_factory=dict)


# hack to run on a Jupyter Notebook (use existing run loop)
def _async_run(coro: Coroutine[Any, Any, None]) -> None:
    class RunThread(threading.Thread):
        def run(self):
            asyncio.run(coro)

    try:
        loop = asyncio.get_running_loop()
    except RuntimeError:
        loop = None

    if loop and loop.is_running():
        thread = RunThread()
        thread.start()
        thread.join()
    else:
        return asyncio.run(coro)


async def _concurrent_subprocess(commands: List[_Command], tasks: int = 1) -> None:
    async def worker():
        while True:
            cmd = await queue.get()
            with _transient_dir(cmd.workdir or os.getcwd()):
                proc = await asyncio.create_subprocess_exec(
                    cmd.args[0],
                    *cmd.args[1:],
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.STDOUT,
                )
                stdout, stderr = await proc.communicate()
                if proc.returncode != 0:
                    raise subprocess.CalledProcessError(
                        proc.returncode or 0, cmd.args, stdout, stderr
                    )
                queue.task_done()

    queue: asyncio.Queue[_Command] = asyncio.Queue()
    for cmd in commands:
        queue.put_nowait(cmd)

    workers: List[asyncio.Task[Any]] = [
        asyncio.create_task(worker()) for _ in range(tasks)
    ]
    await queue.join()

    for task in workers:
        task.cancel()
    await asyncio.gather(*workers, return_exceptions=True)


def _run_silent(
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


@contextlib.contextmanager
def _transient_dir(path: PathLike) -> Generator[None, None, None]:
    Path(path).mkdir(parents=True, exist_ok=True)
    cwd = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(cwd)
