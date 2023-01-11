import contextlib
import enum
import glob
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

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


def report(
    csvfile: PathLike, output_dir: PathLike, resids: str, contact_cutoff: float
) -> pd.DataFrame:
    run_silent(
        os.path.join(SCHRODINGER_PATH, "run"),
        os.path.join(SCRIPT_DIR, "scripts", "report.py"),
        str(output_dir),
        cutoff=contact_cutoff,
        output=str(csvfile),
        residues=resids,
    )
    return pd.read_csv(csvfile)


def ranking_compounds_old(
    common_compounds, df_task1_rank, df_nav_rank, df_kv_cc_rank, df_kv_sp_rank, order
):
    new_df = pd.DataFrame()
    for c in common_compounds:
        task1_row = df_task1_rank.loc[df_task1_rank["NAME"].str.contains(f"{c}-out")]
        nav_row = df_nav_rank.loc[df_nav_rank["NAME"].str.contains(f"{c}-out")]
        kv_cc_row = df_kv_cc_rank.loc[df_kv_cc_rank["NAME"].str.contains(f"{c}-out")]
        kv_sp_row = df_kv_sp_rank.loc[df_kv_sp_rank["NAME"].str.contains(f"{c}-out")]
        new_rank = (
            task1_row.RANK.iloc[0]
            + nav_row.RANK.iloc[0]
            + kv_cc_row.RANK.iloc[0]
            + kv_sp_row.RANK.iloc[0]
        )
        new_norm_total = (
            task1_row.NORMT.iloc[0]
            + nav_row.NORMT.iloc[0]
            + kv_cc_row.NORMT.iloc[0]
            + kv_sp_row.NORMT.iloc[0]
        )
        task1_normt = task1_row.NORMT.iloc[0]
        nav_normt = nav_row.NORMT.iloc[0]
        kv_cc_normt = kv_cc_row.NORMT.iloc[0]
        kv_sp_normt = kv_sp_row.NORMT.iloc[0]
        task1_dgbind = task1_row.DGBIND.iloc[0]
        nav_dgbind = nav_row.DGBIND.iloc[0]
        kv_cc_dgbind = kv_cc_row.DGBIND.iloc[0]
        kv_sp_dgbind = kv_sp_row.DGBIND.iloc[0]
        task1_int = task1_row.INT.iloc[0]
        nav_int = nav_row.INT.iloc[0]
        kv_cc_int = kv_cc_row.INT.iloc[0]
        kv_sp_int = kv_sp_row.INT.iloc[0]
        new_row = {
            "NAME": c,
            "TASK1_INT": task1_int,
            "NAV_INT": nav_int,
            "KV_CC_INT": kv_cc_int,
            "KV_SP_INT": kv_sp_int,
            "TASK1_DGBIND": task1_dgbind,
            "NAV_DGBIND": nav_dgbind,
            "KV_CC_DGBIND": kv_cc_dgbind,
            "KV_SP_DGBIND": kv_sp_dgbind,
            "TASK1_NORMT": task1_normt,
            "NAV_NORMT": nav_normt,
            "KV_CC_NORMT": kv_cc_normt,
            "KV_SP_NORMT": kv_sp_normt,
            "SUM_NORMT": new_norm_total,
            "SUM_RANK": new_rank,
        }
        new_row = pd.DataFrame(new_row, index=[0])
        new_df = pd.concat([new_df, new_row], ignore_index=True)
        if order == "SUM_NORM_TOTAL":
            new_df = new_df.sort_values(by=["SUM_NORM_TOTAL"], ascending=False)
        if order == "SUM_RANK":
            new_df = new_df.sort_values(by=["SUM_RANK"], ascending=True)
    new_df.reset_index(drop=True, inplace=True)
    return new_df


def ranking_compounds(
    common_compounds: list,
    dict_rank_set: dict,
    order: str = "SUM_RANK",
) -> pd.DataFrame:
    new_row = {}
    new_df = pd.DataFrame()
    for c in common_compounds:
        rows_list = list(
            v.loc[v["NAME"].str.contains(f"{c}-out")]
            for k, v in dict_rank_set.items()
            if k.endswith("rank")
        )
        # print(rows_list)
        prot_list = list(
            k.replace("_rank", "")
            for k, v in dict_rank_set.items()
            if k.endswith("rank")
        )
        new_rank = sum(list(v.RANK.iloc[0] for v in rows_list))
        new_norm_total = sum(list(v.NORMT.iloc[0] for v in rows_list))
        normt_list = list(v.NORMT.iloc[0] for v in rows_list)
        dgbind_list = list(v.DGBIND.iloc[0] for v in rows_list)
        int_list = list(v.INT.iloc[0] for v in rows_list)
        new_row[f"NAME"] = c
        for i in range(len(prot_list)):
            new_row[f"{prot_list[i]}_INT"] = int_list[i]
            new_row[f"{prot_list[i]}_DGBIND"] = dgbind_list[i]
            new_row[f"{prot_list[i]}_NORMT"] = normt_list[i]
        new_row[f"SUM_NORMT"] = new_norm_total
        new_row[f"SUM_RANK"] = new_rank

        new_row = pd.DataFrame(new_row, index=[0])
        new_df = pd.concat([new_df, new_row], ignore_index=True)
        if order == "SUM_NORMT":
            new_df = new_df.sort_values(by=["SUM_NORMT"], ascending=False)
        if order == "SUM_RANK":
            new_df = new_df.sort_values(by=["SUM_RANK"], ascending=True)

    new_df.reset_index(drop=True, inplace=True)
    new_df.index = new_df.index + 1
    new_df.index.name = "POSITION"
    return new_df


def common_compounds_old(set1: set, set2: set, set3: set, set4: set) -> list:
    common = list(set.intersection(set1, set2, set3, set4))
    return common


def common_compounds(dict_rank_set: dict) -> list:
    common = list(
        set.intersection(
            *(set(v) for k, v in dict_rank_set.items() if k.endswith("set"))
        )
    )
    return common


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
):
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
):
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
def transient_dir(path: PathLike):
    cwd = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(cwd)


def analysis(
    working_folder: str,
    mmgbsa_output_path: str,
    bs_residues: Dict[str, str],
    radius: float,
    output_folder_name: str,
    order: List[RankingCriterion] = [
        RankingCriterion.NORMALIZED_CONTACTS,
        RankingCriterion.TOTAL_SCORE,
    ],
):
    os.chdir(working_folder)
    Path(f"{output_folder_name}").mkdir(parents=True, exist_ok=True)
    mmgbsa_output_path = os.path.abspath(mmgbsa_output_path)
    proteins = glob.glob(f"{mmgbsa_output_path}/**")

    # Create dict of every system
    csv_data = {}
    for protein in proteins:
        os.chdir(f"{working_folder}/{output_folder_name}")
        protein_name = os.path.basename(f"{protein}")

        Path(f"{protein_name}").mkdir(parents=True, exist_ok=True)
        os.chdir(f"{working_folder}/{output_folder_name}/{protein_name}")

        print(f"Analizing {protein_name}")
        df = report(
            csvfile=os.path.join(
                working_folder, output_folder_name, f"{protein_name}.csv"
            ),
            output_dir=os.path.join(mmgbsa_output_path, protein_name),
            resids=bs_residues[f"{protein_name}.mae"],
            contact_cutoff=radius,
        )

        data_rank = rank_poses(df, order)
        data_set: Set[str] = set(data_rank["NAME"].unique())
        csv_data[f"{protein_name}_rank"] = data_rank
        csv_data[f"{protein_name}_set"] = data_set
        print(" ")

    # Ranking compounds
    common = common_compounds(csv_data)
    print("Common compounds for the systems:", common)
    print("There are", len(common), "common compounds")
    ranking = ranking_compounds(common_compounds=common, dict_rank_set=csv_data)
    return ranking
