# polypharm
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


## Options

Environment variables containing file paths, folder output names, and folder names must first 
be defined: 
* [WORKING_PATH]  
Contains global settings.