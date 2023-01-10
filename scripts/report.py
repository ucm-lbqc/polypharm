import argparse
import csv
import glob
import os
import sys
from typing import Any, Dict, List, Tuple

from schrodinger.structure import StructureReader


def parse_args(argv: List[str]) -> Tuple[Dict[str, Any], str]:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dir", help="Folder containing Maestro (*.mae[gz]) files")
    parser.add_argument("-o", "--output", required=True, help="Output name.")
    opts = parser.parse_args(argv)
    return vars(opts), opts.dir


def main(argv: List[str]) -> None:
    opts, dir = parse_args(argv)
    basename = os.path.basename(dir)

    with open(opts["output"], "w") as io:
        csvfile = csv.writer(io)
        csvfile.writerow(["NAME", "INDEX", "IFDSCORE", "DGBIND"])

        maefiles = sorted(glob.glob(os.path.join(dir, "*.mae*")))
        for maefile in maefiles:
            with StructureReader(maefile) as reader:
                for i, st in enumerate(reader):
                    name = st.property.get("s_m_name", basename)
                    ifdscore = st.property["r_psp_IFDScore"]
                    dg_bind = st.property["r_psp_MMGBSA_dG_Bind"]
                    csvfile.writerow([name, i, ifdscore, dg_bind])


if __name__ == "__main__":
    main(sys.argv[1:])
