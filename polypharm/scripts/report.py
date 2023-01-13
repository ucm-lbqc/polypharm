import argparse
import csv
import glob
import os
import sys
from typing import List, Tuple

from schrodinger.structure import Structure, StructureReader, _StructureAtom
from schrodinger.structutils.analyze import evaluate_asl


class Args:
    def __init__(
        self, maefile: str, resids: str, output: str, cutoff: float = 5.0
    ) -> None:
        self.maefile = maefile
        self.cutoff = cutoff
        self.output = output
        tokens = [resid.split(":") for resid in resids.split(",")]
        self.resids = sorted((ch.strip(), int(num)) for ch, num in tokens)


def get_contact_residues(structure: Structure, cutoff: float) -> List[Tuple[str, int]]:
    idxs: List[int] = evaluate_asl(structure, f"within {cutoff} (ligand)")
    atoms: List[_StructureAtom] = [structure.atom[i] for i in idxs]
    return sorted(set((atom.chain, atom.resnum) for atom in atoms))


def parse_args(argv: List[str]) -> Args:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "maefile",
        help="Maestro (*-out.maegz) file generated by Prime MM/GBSA",
    )
    parser.add_argument(
        "-c",
        "--cutoff",
        type=float,
        default=5,
        help="Contact cutoff. Defaults to %(default)s A.",
    )
    parser.add_argument(
        "-r",
        "--residues",
        dest="resids",
        required=True,
        help="""
        Comma-separated list of residue ids (CHAIN:NUMBER) to report, e.g.,
        'A:81,B:205'.
        """,
    )
    parser.add_argument("-o", "--output", required=True, help="Output CSV filename.")

    opts = parser.parse_args(argv)
    return Args(**vars(opts))


def main(argv: List[str]) -> None:
    opts = parse_args(argv)

    with open(opts.output, "w") as io:
        csvfile = csv.writer(io, lineterminator="\n")

        headers = ["NAME", "INDEX", "IFDSCORE", "DGBIND", "INT"]
        headers.extend([f"{ch.upper()}:{resnum}" for ch, resnum in opts.resids])
        csvfile.writerow(headers)

        basename = os.path.basename(opts.maefile)
        basename = basename.replace(".mae", "").replace("-out", "")
        with StructureReader(opts.maefile) as reader:
            for i, structure in enumerate(reader):
                resids = get_contact_residues(structure, opts.cutoff)
                contact_mask = [int(resid in resids) for resid in opts.resids]
                total_contacts = sum(contact_mask)

                name: str = structure.property.get("s_m_title", basename)
                ifdscore: float = structure.property["r_psp_IFDScore"]
                dg_bind: float = structure.property["r_psp_MMGBSA_dG_Bind"]

                row = [name, i, ifdscore, dg_bind, total_contacts]
                row.extend(contact_mask)
                csvfile.writerow(row)


if __name__ == "__main__":
    main(sys.argv[1:])