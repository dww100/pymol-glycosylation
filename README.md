pymol-glycosylation
===========================================================

This software was developed in the Structural Immunology Group at UCL and was created by David W. Wright (dave.william.wright@gmail.com).

The aim of the software is to provide a simple way to add template glycans to PDB structures.
The link residues (THR or SER for O-linked glycans, ASN for N-linked) and glycans must be provided in template PDB files for all glycans to be added.
PDBs containing templates used in the Structural Immunology Group are provided alonside the code.

The target residues are specified in a YAML file. An example of which is also provided in this repository.

Requirements:
-------------

Pymol - with the python modules in the PYTHONPATH

The following packages must be in the PYTHONPATH:
yaml

