# Polypharm
 A Python-based library to perform IFD and MMGBSA calculations on different targets using a polypharmacological approach.

## Features

It allows to perform Induced Fit Docking (IFD) and binding free energy calculations using the MMGBSA method available in the Schrödinger suite by means of a wrapper library. The protocol performs the calculations on different targets/proteins and considering the experimentally relevant residues creates a scoring function to evaluate the best candidates with activity on the different targets, i.e. it proposes the molecules with the best polypharmacological activity. 

**NOTE: The initial structures must be previously prepared, with added hydrogens, 
correct bond orders and charges, and assignment of protonation states to the desired 
pH. Verify that everything is correct in the initial .mae and .maegz files.**

## Requirements

* Python 3.0 or higher versions.
* Schrödinger software. 

## Usage

a  [jupyter notebook](examples/) was created with the complete protocol, which performs the calculations and analysis to rank the compounds.
The library can be imported into Python using:
```python
import polypharm as ppm
```

## Options

Environment variables containing file paths, folder output names, and folder names must first 
be defined: 

* ``WORKING_PATH``: < Working directory path >   
Contains the full path to the main directory.

* ``SCHRODINGER_PATH``: < Specify the path to the Desmond installation >   
Acceptable values: any path  
Default values: $SCHRODINGER 

* ``LIGANDS_PATH``: < Specify the path to the ligand files >   
Acceptable values: any path  
**NOTE: The ligands must be in .maegz format. We suggest using ligand names as simple as possible, for example: 1a, 1b, 5f, etc. Do not use the word "-out" in ligand names.**

* ``PROTEINS_PATH``: < Specify the path to the target files >   
Acceptable values: any path  
**NOTE: The targets must be in .mae format**

* ``IFD_OUTPUT_NAME``: < Specify the output folder name for IFD calculations >   
Acceptable values: any string  
**We suggest using the name: "03_ifd"**

* ``MMGBSA_OUTPUT_NAME``: < Specify the output folder name for MMGBSA calculations >   
Acceptable values: any string  
**We suggest using the name: "04_mmgbsa"**

* ``BS_RESIDS``: < Specify the binding site residues for the different targets to use >   
Acceptable values: A python dictionary using as key the name of the target including the .mae extension, followed by the residues specifying the chain and residue number with the format ``chain:resid_number,chain:resid_number...``.  
Example:  
```python
BS_RESIDS = {'protein_name_1.mae': "A:1,A:9,A:11,B:7,B:15,B25,C:1",
             'protein_name_2.mae': "B:122,B:126,B:120"}
```

* ``TO_FILTER``: < Specify the name of a ligand to filter from the analysis of the ranking >   
Acceptable values: any string  

* ``MAX_ELEMENTS_PER_CHANNEL``: < Specify the maximum number of compounds to select for the initial ranking of compounds. >   
Acceptable values: any integer.
Default values: 100 

* ``SCHEME_ORDER``: < Specify the name of the scheme to perform the ranking of the compounds.. >   
Acceptable values: 
``INT_NORM-NORMT``:  Perform the ranking considering the normalization of interactions as the first criterion and the total normalization as the second criterion. This scheme gives more weight to the relevant interactions for each target and the binding free energy is considered as the second weight as it is implicit in the total normalization "NORMT". "NORMT" = normalized interactions + normalized binding free energy.
``DGBIND_NORM``: Perform the ranking considering only the normalization of the binding free energy. This option could be useful when choosing targets for which the relevant residues of the binding site are unknown.
``INT-DGBIND_NORM``: Perform the ranking considering the number of interactios directly as a first criterion and the binding free energy normalization as the second criterion.
Default values: ``INT_NORM-NORMT`` 

**NOTE: The criterion used in the original work was: ``INT_NORM-NORMT``. The other options are experimental and have not been fully tested.** 

## License

Licensed under the MIT license, see the separate LICENSE file.

## Citing
* 



