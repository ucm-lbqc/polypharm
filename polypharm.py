import glob
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import jinja2

pd.options.display.max_rows = 100

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


def run_silent(command: str, basename: str):
    with open(os.devnull, "w") as FNULL:
        try:
            subprocess.run(command, shell=True, stdout=FNULL)
        except subprocess.CalledProcessError:
            print(f"Error occured with system {basename}")


def run_ifd(
    prot_files: List[str],
    lig_files: List[str],
    bs_residues: Dict[str, str],
    workdir: str = "ifd",
    glide_cpu: int = 2,
    prime_cpu: int = 2,
):
    template = TEMPLATE_ENV.get_template("ifd.inp")

    for prot_file in map(Path, prot_files):
        prot_name = prot_file.stem

        prot_workdir = Path(workdir, prot_name)
        prot_workdir.mkdir(parents=True, exist_ok=True)
        shutil.copy(prot_file, prot_workdir)

        for i, lig_file in enumerate(map(Path, lig_files)):
            lig_name = lig_file.stem
            lig_workdir = prot_workdir / lig_name

            if lig_workdir.exists():
                print(f"Skipping IFD for {prot_name}/{lig_name}... already exists")
                continue
            else:
                lig_workdir.mkdir()

            shutil.copy(lig_file, lig_workdir)
            with open(lig_workdir / f"{lig_name}.inp", "w") as io:
                vars = dict(
                    protfile=os.path.join("..", prot_file.name),
                    ligfile=lig_file.name,
                    resids=bs_residues[prot_file.name],
                )
                io.write(template.render(vars))

            os.chdir(lig_workdir)
            cmd = f"{SCHRODINGER_PATH}/ifd -NGLIDECPU {glide_cpu} -NPRIMECPU {prime_cpu} {lig_name}.inp -HOST localhost -SUBHOST localhost -TMPLAUNCHDIR -WAIT > /dev/null 2>&1"
            print(f"Running IFD for{prot_name}/{lig_name} ({i}/{len(lig_files)})...")
            run_silent(cmd, lig_name)


def run_mmgbsa(
    working_folder: str,
    ifd_output_path: str,
    output_folder_name: str,
    mmgbsa_cpu=2,
):
    os.chdir(working_folder)
    mmgbsa_exec = f"{SCHRODINGER_PATH}/prime_mmgbsa"
    Path(f"{output_folder_name}").mkdir(parents=True, exist_ok=True)
    ifd_output_path = os.path.abspath(ifd_output_path)
    proteins = glob.glob(f"{ifd_output_path}/**")
    UNK = "'UNK '"

    for protein in proteins:
        os.chdir(f"{working_folder}/{output_folder_name}")
        protein_name = os.path.basename(f"{protein}")
        print(f"Protein: {protein_name}")
        Path(f"{protein_name}").mkdir(parents=True, exist_ok=True)
        os.chdir(f"{working_folder}/{output_folder_name}/{protein_name}")
        ligands = glob.glob(f"{protein}/**/*-out.maegz")
        n_ligands = len(ligands)
        lig_id = 0
        for ligand in ligands:
            lig_id += 1
            ligand_name = os.path.basename(f"{ligand}")
            ligand_basename = os.path.splitext(ligand_name)[0]
            if os.path.exists(f"{ligand_basename}-out.maegz"):
                print(f"{ligand_basename} already exists. Skipping ligand...")
                continue

            mmgbsa_calculation = f'{mmgbsa_exec} -HOST "localhost:{mmgbsa_cpu}" "{ligand}" -j {ligand_basename} -ligand "(res.pt  {UNK})" -csv_output yes -out_type COMPLEX -job_type REAL_MIN -WAIT > /dev/null 2>&1'
            print(
                f"Running MMGBSA for ligand {ligand_basename} in {protein_name} ({lig_id}/{n_ligands})..."
            )
            run_silent(mmgbsa_calculation, ligand_basename)


def analysis(
    working_folder: str,
    mmgbsa_output_path: str,
    bs_residues: Dict[str, str],
    radius: float,
    output_folder_name: str,
    order: str = "INT_NORM-NORMT",
):
    os.chdir(working_folder)
    run_exec = f"{SCHRODINGER_PATH}/run"
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
        cmd = f"{run_exec} {SCRIPT_DIR}/scripts/report.py -c {radius} -o {csvfile} -r {resids} {mmgbsa_output_path}/{protein_name}"
        run_silent(cmd, protein_name)

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
