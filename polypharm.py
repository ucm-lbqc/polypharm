import pandas as pd
import subprocess
import os
import glob
from pathlib import Path
import shutil
from typing import Tuple

pd.options.display.max_rows = 100

def reading_raw_data(
    mmgbsa_file: str,
    interactions_file: str,
    residues: list,
    to_filter: str,
    order: str = "INT_NORM-NORMT",
) -> pd.DataFrame:
    print("Reading data...")
    protein_name = os.path.basename(f"{mmgbsa_file}")
    protein_name = protein_name.replace('_mmgbsa.csv', '')
    residues = residues[f'{protein_name}.mae'].split(",")
    residues.append("NAME")
    df_mmgbsa = pd.read_csv(f"{mmgbsa_file}")
    print("MMGBSA data rows:", df_mmgbsa.shape[0])
    df_interactions = pd.read_csv(f"{interactions_file}")
    print("Interactions data rows:", df_interactions.shape[0])
    df_interactions.rename(columns=str.upper, inplace=True)

    # Erase columns that are not part of the binding-site residues
    for column in df_interactions:
        if column not in residues:
            df_interactions.drop([f"{column}"], axis=1, inplace=True)

    df = pd.merge(df_mmgbsa, df_interactions, left_on="NAME", right_on="NAME")
    print("Total rows (merged):", df.shape[0])
    check_rows(df_interactions.shape[0], df_mmgbsa.shape[0], df.shape[0])

    # Filter the poses with the option=to_filter.
    if to_filter != "":
        print("to_Filtering...")
        df = df[df["NAME"].str.contains(to_filter) == False]

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


def run_mmgbsa(
    working_folder: str,
    schrodinger_path: str,
    ifd_output_path: str,
    output_folder_name: str,
    mmgbsa_cpu=2,
):
    os.chdir(working_folder)
    schrod_path = os.path.abspath(schrodinger_path)
    mmgbsa_exec = f"{schrod_path}/prime_mmgbsa"
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


def write_extractor():
    with open("extractor.py", "w") as f:
        f.write(
            """from __future__ import print_function
import argparse
import os
import sys
import glob
from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.utils.fileutils import get_jobname
from schrodinger.structutils.analyze import evaluate_asl

def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", help="Path to the (*.mae[gz]) files")
    parser.add_argument("-o", "--output", help="Output name for structures.")
    parser.add_argument(
        "-asl",
        dest="asl_expr",
        required=False,
        help="Only atoms selected by the given ASL expression will " "be extracted.",
    )
    opts = parser.parse_args(argv)
    return vars(opts), opts.path
def main(argv):
    opts, path = parse_args(argv)
    maefiles = sorted(glob.glob(os.path.join(path, "*.maegz")))
    endname = "%s" % opts["output"]
    for maefile in maefiles:
        output_name = get_jobname(maefile)
        reader = StructureReader(maefile)
        count = 0
        for st in reader:
            count += 1
            st_number = str(count)
            output_name2 = "{}_{}".format(output_name, st_number)
            data1 = "{:<20}".format(output_name2)
            # Adding property
            name_property = st.property.get("s_m_name")
            if name_property is not None:
                pass
            else:
                st.property["s_m_name"] = output_name2
            print(data1)
            # Writing ligands
            if opts["asl_expr"] is not None:
                indices1 = evaluate_asl(st, "(%s)" % opts["asl_expr"])
                st1 = st.extract(indices1, copy_props=True)
                st1.title = str(output_name2)
                writer = StructureWriter(output_name2 + ".mae")
                writer.append(st1)
                writer.close()
            else:
                pass
        reader.close()
        writer.close()
if __name__ == "__main__":
    main(sys.argv[1:])
    """
        )


def write_props_reader():
    with open("props_reader.py", "w") as f:
        f.write(
            """from __future__ import print_function

import argparse
import os
import sys
import glob
from schrodinger.structure import StructureReader

def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", help="Path to the (*.mae[gz]) files")
    parser.add_argument("-o", "--output", help="Output name.")
    opts = parser.parse_args(argv)
    return vars(opts), opts.path

def main(argv):
    opts, path = parse_args(argv)
    maefiles = sorted(glob.glob(os.path.join(path, "*.mae")))
    count = 0
    filename = "%s.csv" % opts["output"]
    data = open(filename, "w+")
    data.write("{},{},{}".format("NAME", "IFDSCORE", "DGBIND"))
    data.write("\\n")
    print("reading properties")
    for maefile in maefiles:
        reader = StructureReader(maefile)
        for st in reader:
            name_property = st.property.get("s_m_name")
            print(name_property)
            ifdscore = st.property["r_psp_IFDScore"]
            dg_bind = st.property["r_psp_MMGBSA_dG_Bind"]
            count += 1
            data1 = f"{name_property},{ifdscore},{dg_bind}"
            data.write(data1)
            data.write("\\n")
    data.close()
    reader.close()

if __name__ == "__main__":
    main(sys.argv[1:])
    """
        )


def write_vmd_script(bs_residues: dict, protein_name: str, radius: float):
    jump = "\n"
    resids_string = bs_residues[protein_name]
    numbers = set()
    chains = set()
    resids_string = resids_string.split(",")
    for pair in resids_string:
        chains.add(pair.split(":")[0])
        numbers.add(pair.split(":")[1])

    chains = sorted(chains)
    numbers = sorted(numbers)
    resids = ""
    for x in numbers:
        resids += " " + str(x)
    resids_number = resids
    resids_number = f"{ { resids_number } }"
    resids_number = resids_number.replace("'", "")

    chs = ""
    for x in chains:
        chs += " " + '"' + x + '"'
    chains = chs
    chains = f"{ { chains } }"
    chains = chains.replace("'", "")

    with open("interactions.tcl", "w+") as f:
        f.write(
            """set folder [lindex $argv 0]
set output [lindex $argv 1]
set pdbfiles [glob $folder/*.mae]
set ligand "UNK"
"""
        )
        f.write(f"set resids {resids_number} {jump}")
        f.write(f"set radius {radius} {jump}")
        f.write(f"set chains {chains} {jump}")
        f.write(
            """
set outfile [open ${output}_interactions.csv w]

puts -nonewline $outfile "Name"
foreach resid $resids {
    foreach chain $chains {
        puts -nonewline $outfile ",$chain:$resid"
    }
}
puts $outfile ""


foreach file $pdbfiles {
    set basename [file rootname [file tail $file]]
    puts -nonewline $outfile $basename
    mol new $file
    foreach resid $resids {
        foreach chain $chains {
            set sel [atomselect top "same residue as (resid $resids and chain $chain and within $radius of resname $ligand)"]
            set current_resids [lsort -unique [$sel get resid]]
            if {[lsearch $current_resids $resid] > -1} {
                # lset mat $i [expr $j*[llength $chains]+$offset] 1
                puts -nonewline $outfile ",1"
            } else {
                puts -nonewline $outfile ",0"
            }
        }
    }
    puts $outfile ""
    mol delete top
}

close $outfile

exit
    """
        )


def analysis(
    working_folder: str,
    schrodinger_path: str,
    mmgbsa_output_path: str,
    bs_residues: str,
    radius: float,
    output_folder_name: str,
    extract_poses: bool= True,
    extract_mmgbsa_info: bool= True,
    calculate_interactions: bool=True,
):
    write_extractor()
    write_props_reader()

    os.chdir(working_folder)
    schrod_path = os.path.abspath(schrodinger_path)
    run_exec = f"{schrod_path}/run"
    Path(f"{output_folder_name}").mkdir(parents=True, exist_ok=True)
    mmgbsa_output_path = os.path.abspath(mmgbsa_output_path)
    proteins = glob.glob(f"{mmgbsa_output_path}/**")

    # Create dict of every system
    csv_data = {}
    for protein in proteins:
        os.chdir(f"{working_folder}/{output_folder_name}")
        protein_name = os.path.basename(f"{protein}")
        write_vmd_script(
            bs_residues=bs_residues, protein_name=f"{protein_name}.mae", radius=radius
        )
        Path(f"{protein_name}").mkdir(parents=True, exist_ok=True)
        os.chdir(f"{working_folder}/{output_folder_name}/{protein_name}")

        print(f"Analizing {protein_name}")
        if extract_poses:
            print(f"Extracting poses...")
            extract_poses_cmd = f'{run_exec} {working_folder}/extractor.py {protein} -asl "all"'
            run_silent(extract_poses_cmd, protein_name)
        if extract_mmgbsa_info:
            print(f"Extracting MMGBSA info...")
            extract_mmgbsa_info_cmd = f"{run_exec} {working_folder}/props_reader.py {working_folder}/{output_folder_name}/{protein_name} -o {protein_name}_mmgbsa"
            run_silent(extract_mmgbsa_info_cmd, protein_name)
            
        if calculate_interactions:
            print(f"Calculating interactions...")
            calculate_interactions_cmd = f"vmd -dispdev none -e {working_folder}/{output_folder_name}/interactions.tcl -args {working_folder}/{output_folder_name}/{protein_name}  {protein_name} > /dev/null 2>&1"
            run_silent(calculate_interactions_cmd, protein_name)

        # ranking poses of every system.
        df = reading_raw_data(mmgbsa_file=f'{protein_name}_mmgbsa.csv', interactions_file=f'{protein_name}_interactions.csv', residues=bs_residues, to_filter="")
        data_rank, data_set = ranking_poses(df=df, max_rank=100)
        csv_data[f'{protein_name}_rank'] = data_rank
        csv_data[f'{protein_name}_set'] = data_set
        print(" ")

    # Ranking compounds
    common = common_compounds(csv_data)
    print("Common compounds for the systems:", common)
    print("There are", len(common), "common compounds")
    ranking = ranking_compounds(common_compounds=common, dict_rank_set=csv_data)
    return ranking

