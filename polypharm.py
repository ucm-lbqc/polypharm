import pandas as pd
import subprocess
import os
import glob
from pathlib import Path
import shutil
from typing import Tuple

pd.options.display.max_rows = 100
def run_silent(command: str, basename: str):
    with open(os.devnull, "w") as FNULL:
        try:
            subprocess.run(command, shell=True, stdout=FNULL)
        except subprocess.CalledProcessError:
            print(f"Error occured with system {basename}")


def run_ifd(
    working_folder: str,
    schrodinger_path: str,
    ligands_path: str,
    proteins_path: str,
    output_folder_name: str,
    bs_residues: list,
    nglide_cpu=2,
    nprime_cpu=2,
):
    os.chdir(working_folder)
    schrod_path = os.path.abspath(schrodinger_path)
    ifd_exec = f"{schrod_path}/ifd"
    Path(f"{output_folder_name}").mkdir(parents=True, exist_ok=True)
    ligands_path = os.path.abspath(ligands_path)
    proteins_path = os.path.abspath(proteins_path)

    ligands = glob.glob("{}/*.{}".format(ligands_path, "maegz"))
    proteins = glob.glob("{}/*.{}".format(proteins_path, "mae"))
    jump = "\n"
    tab = "  "
    n_ligands = len(ligands)

    for protein in proteins:
        os.chdir(f"{working_folder}/{output_folder_name}")
        protein_name = os.path.basename(f"{protein}")
        protein_basename = os.path.splitext(protein_name)[0]
        Path(f"{protein_basename}").mkdir(parents=True, exist_ok=True)
        shutil.copy(f"{protein}", f"{protein_basename}/")
        lig_id = 0
        for ligand in ligands:
            lig_id += 1
            os.chdir(f"{working_folder}/{output_folder_name}/{protein_basename}")
            ligand_name = os.path.basename(f"{ligand}")
            ligand_basename = os.path.splitext(ligand_name)[0]
            if os.path.exists(f"{ligand_basename}"):
                print(
                    f"Folder {ligand_basename} already exists, calculation skipped ..."
                )
                continue
            else:
                Path(f"{ligand_basename}").mkdir(parents=True, exist_ok=True)

            shutil.copy(f"{ligand}", f"{ligand_basename}/")
            os.chdir(
                f"{working_folder}/{output_folder_name}/{protein_basename}/{ligand_basename}"
            )
            # Write input for ifd by every ligand
            with open(f"{ligand_basename}.inp", "w") as f:
                f.write(f"INPUT_FILE ../{protein_name} {jump}")
                f.write(f"STAGE GLIDE_DOCKING2 {jump}")
                f.write(
                    f'{tab}BINDING_SITE residues {bs_residues[f"{protein_name}"]}  {jump}'
                )
                f.write(f"{tab}INNERBOX 10.0 {jump}")
                f.write(f"{tab}OUTERBOX auto {jump}")
                f.write(f"{tab}LIGAND_FILE {ligand_name} {jump}")
                f.write(f"{tab}LIGANDS_TO_DOCK all {jump}")
                f.write(f"{tab}GRIDGEN_RECEP_CCUT  0.25 {jump}")
                f.write(f"{tab}GRIDGEN_RECEP_VSCALE 0.50 {jump}")
                f.write(f"{tab}DOCKING_PRECISION SP {jump}")
                f.write(f"{tab}DOCKING_LIG_CCUT  0.15 {jump}")
                f.write(f"{tab}DOCKING_CV_CUTOFF  100.0 {jump}")
                f.write(f"{tab}DOCKING_LIG_VSCALE 0.50 {jump}")
                f.write(f"{tab}DOCKING_POSES_PER_LIG 10 {jump}")
                f.write(f"{tab}DOCKING_RINGCONFCUT 2.5 {jump}")
                f.write(f"{tab}DOCKING_AMIDE_MODE penal {jump}")
                f.write(f"{jump}")
                f.write(f"STAGE COMPILE_RESIDUE_LIST {jump}")
                f.write(f"{tab}DISTANCE_CUTOFF 5.0 {jump}")
                f.write(f"{jump}")
                f.write(f"STAGE PRIME_REFINEMENT {jump}")
                f.write(f"{tab}NUMBER_OF_PASSES 1 {jump}")
                f.write(f"{tab}USE_MEMBRANE no {jump}")
                f.write(f"{tab}OPLS_VERSION OPLS_2005 {jump}")
                f.write(f"{jump}")
                f.write(f"STAGE SORT_AND_FILTER {jump}")
                f.write(f"{tab}POSE_FILTER r_psp_Prime_Energy {jump}")
                f.write(f"{tab}POSE_KEEP 30.0 {jump}")
                f.write(f"{jump}")
                f.write(f"STAGE SORT_AND_FILTER {jump}")
                f.write(f"{tab}POSE_FILTER r_psp_Prime_Energy {jump}")
                f.write(f"{tab}POSE_KEEP 10# {jump}")
                f.write(f"{jump}")
                f.write(f"STAGE GLIDE_DOCKING2 {jump}")
                f.write(f"{tab}BINDING_SITE ligand Z:999 {jump}")
                f.write(f"{tab}INNERBOX 10.0 {jump}")
                f.write(f"{tab}OUTERBOX auto {jump}")
                f.write(f"{tab}LIGAND_FILE {ligand_name} {jump}")
                f.write(f"{tab}LIGANDS_TO_DOCK self {jump}")
                f.write(f"{tab}GRIDGEN_RECEP_CCUT  0.25 {jump}")
                f.write(f"{tab}GRIDGEN_RECEP_VSCALE 1.00 {jump}")
                f.write(f"{tab}DOCKING_PRECISION SP {jump}")
                f.write(f"{tab}DOCKING_LIG_CCUT  0.15 {jump}")
                f.write(f"{tab}DOCKING_CV_CUTOFF  0.0 {jump}")
                f.write(f"{tab}DOCKING_LIG_VSCALE 0.80 {jump}")
                f.write(f"{tab}DOCKING_POSES_PER_LIG 1 {jump}")
                f.write(f"{tab}DOCKING_RINGCONFCUT 2.5 {jump}")
                f.write(f"{tab}DOCKING_AMIDE_MODE penal {jump}")
                f.write(f"{jump}")
                f.write(f"STAGE SCORING {jump}")
                f.write(f"{tab}SCORE_NAME r_psp_IFDScore {jump}")
                f.write(f"{tab}TERM 1.0,r_i_glide_gscore,0 {jump}")
                f.write(f"{tab}TERM 0.05,r_psp_Prime_Energy,1 {jump}")
                f.write(f"{tab}REPORT_FILE report.csv {jump}")
            # Running IFD
            ifd_calculation = f"{ifd_exec} -NGLIDECPU {nglide_cpu} -NPRIMECPU {nprime_cpu} {ligand_basename}.inp -HOST localhost -SUBHOST localhost -TMPLAUNCHDIR -WAIT > /dev/null 2>&1"
            print(
                f"Running IFD for ligand {ligand_basename} in {protein_name} ({lig_id}/{n_ligands})..."
            )
            run_silent(ifd_calculation, ligand_basename)
