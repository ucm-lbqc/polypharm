import contextlib
import glob
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Union

import pandas as pd
import jinja2

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


def ranking_poses(df: pd.DataFrame, max_rank: int) -> Tuple[pd.DataFrame, set]:
    df_rank = pd.DataFrame()
    compound_set = set()
    n = 0
    for i in range(len(df)):
        compound = df.loc[i, "NAME"].split("-out")[0]
        if len(compound_set) == max_rank:
            break
        if compound not in compound_set:
            n += 1
            compound_set.add(compound)
            name = df.loc[i, "NAME"]
            int = df.loc[i, "INT"]
            int_norm = df.loc[i, "INT_NORM"]
            norm_dgbind = df.loc[i, "DGBIND_NORM"]
            dgbind = df.loc[i, "DGBIND"]
            norm_total = df.loc[i, "NORMT"]
            new_row = {
                "NAME": name,
                "DGBIND": dgbind,
                "INT": int,
                "DGBIND_NORM": norm_dgbind,
                "INT_NORM": int_norm,
                "NORMT": norm_total,
                "RANK": n,
            }
            new_row = pd.DataFrame(new_row, index=[0])
            df_rank = pd.concat([df_rank, new_row], ignore_index=True)
            df_rank.index = df_rank.index + 1
            df_rank.index.name = "POSITION"
    return df_rank, compound_set


def reading_raw_data(
    csvfile: str,
    order: str = "INT_NORM-NORMT",
) -> pd.DataFrame:
    print("Reading data...")
    df: pd.DataFrame = pd.read_csv(csvfile)
    df["INT"] = df.iloc[:, 3:].sum(axis=1)
    df = df.sort_values(by=["INT"], ascending=False)

    # Interactions and energy normalization
    best_dgbind = df.DGBIND.min()
    worst_dgbind = df.DGBIND.max()
    worst_int = df.INT.min()
    best_int = df.INT.max()
    print(f"Better dGbind: {df.DGBIND.min()}")
    print(f"Worst dGbind: {df.DGBIND.max()}")
    print(f"Better interactions number: {best_int}")
    print(f"Worst interactions number: {worst_int}")

    df["INT_NORM"] = df.apply(
        lambda row: (row.INT - worst_int) / (best_int - worst_int), axis=1
    )
    if worst_dgbind < 0 and best_dgbind < 0:
        worst_dgbind = abs(worst_dgbind)
        best_dgbind = abs(best_dgbind)
        df["DGBIND_NORM"] = df.apply(
            lambda row: (abs(row.DGBIND) - worst_dgbind) / (best_dgbind - worst_dgbind),
            axis=1,
        )

    elif best_dgbind < 0 and worst_dgbind > 0:
        worst_dgbind = worst_dgbind
        best_dgbind = abs(best_dgbind)
        df["DGBIND_NORM"] = df.apply(
            lambda row: (abs(row.DGBIND) - worst_dgbind) / (best_dgbind - worst_dgbind),
            axis=1,
        )

    elif best_dgbind >= 0 and worst_dgbind > 0:
        worst_dgbind = worst_dgbind
        best_dgbind = best_dgbind
        df["DGBIND_NORM"] = df.apply(
            lambda row: (
                ((row.DGBIND) - best_dgbind) / (best_dgbind - worst_dgbind) * -1 + 1
            ),
            axis=1,
        )

    df["NORMT"] = df.apply(lambda row: row.INT_NORM + row.DGBIND_NORM, axis=1)

    if order == "INT-NORMT":
        df = df.sort_values(by=["INT", "NORMT"], ascending=False)
    if order == "INT_NORM-NORMT":
        df = df.sort_values(by=["INT_NORM", "NORMT"], ascending=False)
    if order == "DGBIND_NORM":
        df = df.sort_values(by=["DGBIND_NORM"], ascending=True)
    if order == "DGBIND":
        df = df.sort_values(by=["DGBIND"], ascending=True)
    if order == "INT-DGBIND_NORM":
        df = df.sort_values(by=["INT", "DGBIND_NORM"], ascending=True)

    # Replace -out-out by _out_out to ignore the "-out-out" words produced by the ifd and MMGBSA output
    # df["NAME"] = df["NAME"].str.replace("-out", "*out")
    df.reset_index(drop=True, inplace=True)
    return df


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
            run_silent(
                os.path.join(SCHRODINGER_PATH, "prime_mmgbsa"),
                str(ifd_file),
                ligand="(res.pt UNK)",
                job_type="REAL_MIN",
                out_type="COMPLEX",
                csv_output="yes",
                j=lig_name,
                HOST=f"localhost:{cpus}",
                WAIT=True,
                use_single_dash=True,
            )


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
    order: str = "INT_NORM-NORMT",
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
        resids = bs_residues[f"{protein_name}.mae"]

        Path(f"{protein_name}").mkdir(parents=True, exist_ok=True)
        os.chdir(f"{working_folder}/{output_folder_name}/{protein_name}")

        print(f"Analizing {protein_name}")
        csvfile = f"{working_folder}/{output_folder_name}/{protein_name}.csv"
        run_silent(
            os.path.join(SCHRODINGER_PATH, "run"),
            os.path.join(SCRIPT_DIR, "scripts", "report.py"),
            os.path.join(mmgbsa_output_path, protein_name),
            cutoff=radius,
            output=csvfile,
            residues=resids,
        )

        # ranking poses of every system.
        df = reading_raw_data(csvfile, order)
        data_rank, data_set = ranking_poses(df=df, max_rank=100)
        csv_data[f"{protein_name}_rank"] = data_rank
        csv_data[f"{protein_name}_set"] = data_set
        print(" ")

    # Ranking compounds
    common = common_compounds(csv_data)
    print("Common compounds for the systems:", common)
    print("There are", len(common), "common compounds")
    ranking = ranking_compounds(common_compounds=common, dict_rank_set=csv_data)
    return ranking
