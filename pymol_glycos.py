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
    parser.add_argument('-r','--resid', nargs='+', type=str, dest='targetResidues', help = 'Residues to glycosylate. Format = XN; X = Chain ID, N = Residue number.', required=True)
    parser.add_argument('-g','--glycanpath', nargs='?', type=str, dest='glycanPath', help = 'Path to templates of glycan o-linked to residue', required=True)
    args = parser.parse_args()
    return args
    
def residue_list_parse(combined_list):
    split_list = []
    for res in combined_list:
        if res[0].isalpha() and res[1:].isdigit():
            split_list.append([res[0], res[1:]])
        else:
            print ("Residues must be specified in the following format = XN; X = Chain ID, N = Residue number.")
            sys.exit(1)
    return split_list
    
    
pymol.finish_launching()

# --------------------
# Start of main script
# --------------------

# Interpret command line arguments
args = parse_arguments()
# Use the file name to provide the name for the target structure object
targetStructure = args.inPDB.split('/')[-1].split('.')[0]

# Load targetPDB structure
pymol.cmd.load(args.inPDB, targetStructure)

# Define the selection used to superimpose glycosylated 
# template with the target THR or SER 
atomMask = ' and name CA+N+C+O+CB'

# Split input list of residue identifications into a list of
# lists containing separated chain ID and residue number
targetResList = residue_list_parse(args.targetResidues)

for residue in targetResList:
    
    chainID = residue[0]
    resNo = residue[1]
    strResNo = str(resNo)

    # O-links only occur on SER and THR residues 
    # Check to see that a valid residue type was specified
    targetResidue = targetStructure + " and chain " + chainID + " and resi " + strResNo    
    resType = pymol.cmd.get_model(targetResidue + " and chain " + chainID, 1).atom[0].resn
    
    if resType == 'SER' or resType == 'THR':
        
        # Load template glycan (o-linked to a SER or THR as appropriate)
        # The name of the resulting pymol pbject is stored in glycanStructure
        if resType == 'SER':
            glycanPDB = os.path.join(args.glycanPath,'ser_o-link.pdb')
        else:
            glycanPDB = os.path.join(args.glycanPath,'thr_o-link.pdb')
        glycanStructure = 'glycan'
        pymol.cmd.load(glycanPDB, glycanStructure)

        # Align the linked residue in the template with the target residue
        # using the atoms selected in atomMask
        targetAlignSelection = targetResidue + atomMask
        linkResidue = glycanStructure + " and resn " + resType
        linkResidueSelection = linkResidue + atomMask
        glycanSelection = glycanStructure + " and (not resn " + resType + ")"
        pymol.cmd.align(linkResidueSelection, targetAlignSelection)

        # Remove the target residue in the original structure
        pymol.cmd.remove(targetResidue)
        # Renumber the linked residue in the glycan link template to that of
        # the removed target residue
        pymol.cmd.alter(linkResidue, 'resi = ' + strResNo)
        # Alter the chainID of the link and glycan
        pymol.cmd.alter(glycanStructure, "chain = '" +  chainID + "'")
        # Give the glycan residues a residue number after that of the existing atoms
        # of the chain containing the target link residue
        newResNo = int(pymol.cmd.get_model(targetStructure + " and chain " + chainID, 1).atom[-1].resi) + 1    
        firstGlycanRes = int(pymol.cmd.get_model(glycanSelection,1).atom[0].resi)
        lastGlycanRes = int(pymol.cmd.get_model(glycanSelection,1).atom[-1].resi)
        
        for oldResNo in range(firstGlycanRes, lastGlycanRes + 1):
            pymol.cmd.alter(glycanSelection + ' and resi ' + str(oldResNo), 'resi = ' + str(newResNo))
            newResNo += 1

        # Create a new structure containing combining the target and the new link
        # residue and linked glycan
        pymol.cmd.create('glycan_added', glycanStructure + ' or ' + targetStructure)
        pymol.cmd.delete(glycanStructure)
        pymol.cmd.delete(targetStructure)
        # Rename the newly constructed structure to become the new target
        # for any subsequent glycan additions
        pymol.cmd.set_name('glycan_added',targetStructure)
    else:
        print "Selection '" + targetResidue + "' is not a SER or THR residue and has been ignored."

# Sort and save the glycosylated structure
pymol.cmd.sort(targetStructure)
pymol.cmd.save(args.outPDB, targetStructure)


# Get out!
pymol.cmd.quit()