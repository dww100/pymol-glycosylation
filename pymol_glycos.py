#!/bin/env python
#
#
 
import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
 
import sys, time, os
import argparse
import pymol

def parse_arguments():
	"""Parse command line arguments and ensure correct combinations present"""
	#Command line option parsing
	parser = argparse.ArgumentParser(description= 'Add o=linked glycan to specified SER or THR in a PDB')
	parser.add_argument('-i','--infile', nargs='?', type=str, dest='inPDB', help = 'Path to the input PDB', required=True)
	parser.add_argument('-o','--outfile', nargs='?', type=str, dest='outPDB', help = 'Path to the output PDB', required=True)
	parser.add_argument('-r','--resid', nargs='?', type=int, dest='targetResidue', help = 'Residue number of residue to glycosylate', required=True)
	parser.add_argument('-c','--chain', nargs='?', type=str, dest='targetChain', help = 'Chain letter of Residue to glycosylate', required=True)
	parser.add_argument('-g','--glycanfile', nargs='?', type=str, dest='glycanPDB', help = 'Template of glycan o-linked to residue', required=True)
	args = parser.parse_args()
	return args

pymol.finish_launching()

# --------------------
# Start of main script
# --------------------

# Interpret command line arguments
args = parse_arguments()
targetName = args.inPDB.split('/')[-1].split('.')[0]
glycanName = 'glycan'

# Load structures - targetPDB and glycan
pymol.cmd.load(args.inPDB, targetName)
pymol.cmd.load(args.glycanPDB, glycanName)

# Define the selection used to superimpose glycosylated 
# template with the target THR or SER 
atomMask = ' and name CA+N+C+O+CB'
targetResidue = targetName + " and chain " + args.targetChain + " and resi " + str(args.targetResidue)
targetAlignSelection = targetResidue + atomMask
linkResidue = glycanName + " and resn THR"
linkResidueSelection = linkResidue + atomMask
glycanSelection = glycanName + " and (not resn THR)"

pymol.cmd.align(linkResidueSelection, targetAlignSelection)

targetChain = pymol.cmd.get_model(targetResidue, 1).atom[0].chain
pymol.cmd.remove(targetResidue)
pymol.cmd.alter(linkResidue, 'resi = ' + str(args.targetResidue))
pymol.cmd.alter(glycanName, "chain = '" +  args.targetChain + "'")

lastRes = int(pymol.cmd.get_model(targetName + " and chain A", 1).atom[-1].resi) + 1
pymol.cmd.alter(glycanSelection, 'resi = ' + str(lastRes))

pymol.cmd.create('glycan_added', glycanName + ' or ' + targetName)
pymol.cmd.sort
pymol.cmd.save(args.outPDB, 'glycan_added')


# Get out!
pymol.cmd.quit()