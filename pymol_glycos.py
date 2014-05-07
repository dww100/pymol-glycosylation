#!/usr/bin/env python
#
#

import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

import sys, time, os
import argparse
import pymol
import yaml

def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""

    parser = argparse.ArgumentParser(description= 'Add glycan to specified residues in a PDB: SER or THR (o-linked), ASN (n-linked)')
    parser.add_argument('-i','--input', nargs='?', type=str, dest='in_pdb', 
                        help = 'Path to the input PDB', required=True)
    parser.add_argument('-o','--output', nargs='?', type=str, dest='out_pdb', 
                        help = 'Path to the output PDB', required=True)
    parser.add_argument('-t','--targetfile', nargs='?', type=str, dest='target_file',
                        help = 'YAML file describing the residues to glycosylate.', required=True)
    parser.add_argument('-g','--glycanpath', nargs='?', type=str, dest='glycan_path',
                        help = 'Path to templates of glycan o-linked to residue', default = '.')
    args = parser.parse_args()
    return args

pymol.finish_launching()

# --------------------
# Start of main script
# --------------------



# Setting to write CONECT records to PDB
pymol.cmd.set("pdb_conect_all", "on")
# Setting to read CONECT records in PDBs
pymol.cmd.set("connect_mode", 0)

# Interpret command line arguments
args = parse_arguments()
# Use the file name to provide the name for the target structure object
target_structure = args.in_pdb.split('/')[-1].split('.')[0]

# Load targetPDB structure
pymol.cmd.load(args.in_pdb, target_structure)

# Define the selection used to superimpose glycosylated
# template with the target ASN, THR or SER
atom_mask = ' and not name N+C+O'

# Read in YAML file containing the target residues
f = open(args.target_file, 'r')
target_res_list = yaml.load(f)

glycan_residues = []

for residue in target_res_list:

    chain_id = residue['chain']
    res_no = residue['resid']
    str_res_no = str(res_no)

    # O-links only occur on SER and THR residues
    # N- links at ASN residues
    # Check to see that a valid residue type was specified
    target_residue = target_structure + " and chain " + chain_id + " and resi " + str_res_no
    res_type = pymol.cmd.get_model(target_residue + " and chain " + chain_id, 1).atom[0].resn

    if res_type in ['SER', 'THR', 'ASN']:

        # Load template glycan (o-linked to a SER or THR as appropriate)
        # The name of the resulting pymol pbject is stored in glycan_structure
        if res_type == 'SER':
            glycan_pdb = os.path.join(args.glycan_path,'ser_o-link.pdb')
        elif res_type == 'THR':
            glycan_pdb = os.path.join(args.glycan_path,'thr_o-link.pdb')
        else:
            glycan_pdb = os.path.join(args.glycan_path,'n-link.pdb')
        glycan_structure = 'glycan'
        pymol.cmd.load(glycan_pdb, glycan_structure)

        # Align the linked residue in the template with the target residue
        # using the atoms selected in atom_mask
        target_align_selection = target_residue + atom_mask
        linkResidue = glycan_structure + " and resn " + res_type
        link_res_selection = linkResidue + atom_mask
        glycan_selection = glycan_structure + " and (not resn " + res_type + ")"
        pymol.cmd.align(link_res_selection, target_align_selection)

        # Remove the target residue in the original structure
        pymol.cmd.remove(target_residue)
        # Renumber the linked residue in the glycan link template to that of
        # the removed target residue
        pymol.cmd.alter(linkResidue, 'resi = ' + str_res_no)
        # Alter the chain_id of the link and glycan
        pymol.cmd.alter(glycan_structure, "chain = '" +  chain_id + "'")
        # Give the glycan residues a residue number after that of the existing atoms
        # of the chain containing the target link residue
        new_res_no = int(pymol.cmd.get_model(target_structure + " and chain " + chain_id, 1).atom[-1].resi) + 1
        first_glycan_res = int(pymol.cmd.get_model(glycan_selection,1).atom[0].resi)
        last_glycan_res = int(pymol.cmd.get_model(glycan_selection,1).atom[-1].resi)

        for old_res_no in range(first_glycan_res, last_glycan_res + 1):
            pymol.cmd.alter(glycan_selection + ' and resi ' + str(old_res_no), 'resi = ' + str(new_res_no))
            glycan_residues.append(new_res_no)
            new_res_no += 1

        # Create a new structure containing combining the target and the new link
        # residue and linked glycan
        pymol.cmd.create('glycan_added', glycan_structure + ' or ' + target_structure)
        pymol.cmd.delete(glycan_structure)
        pymol.cmd.delete(target_structure)
        # Rename the newly constructed structure to become the new target
        # for any subsequent glycan additions
        pymol.cmd.set_name('glycan_added',target_structure)
    else:
        print "Selection '" + target_residue + "' is not a ASN, SER or THR residue and has been ignored."

# Sort and save the glycosylated structure
pymol.cmd.sort(target_structure)

pymol.cmd.save(args.out_pdb, target_structure)

out = open('tmp.pdb','w')

glycan_index = []

with open(args.out_pdb, 'r') as f:
    for line in f:
        if (line[0:6] in ['ATOM  ', 'HETATM']):
            if int(line[22:26]) in glycan_residues:
                glycan_index.append(line[6:11])
            out.write(line)   
        elif line[0:6] == 'CONECT':
            if line[6:11] in glycan_index:
                out.write(line)
        else:
            out.write(line)
out.close()

os.rename('tmp.pdb',args.out_pdb)

# Get out!
pymol.cmd.quit()


