#!/usr/bin/env python
"""
Script to add 'standard' glycans to proteins.
Pymol is used to edit the structures and generate CONECT records.
Created models are suitable as input to the Glycan Reader function of
CHARMM-GUI in order to generate CHARMM PSF/PDB pairs for simulation.

Note: The method is totally naive and matches the atoms of the link residue
(except for the backbone N, C and O) in the target PDB to that in the template.
"""

# Copyright 2014 University College London

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
                        help = 'Path to YAML file describing the residues to glycosylate.', 
                        default = None)
    parser.add_argument('-e','--existingfile', nargs='?', type=str, dest='exist_file',
                        help = 'Path to YAML file describing existing glycan linkages', 
                        default = None)
    parser.add_argument('-g','--glycanpath', nargs='?', type=str, dest='glycan_path',
                        help = 'Path to templates of glycan o-linked to residue', 
                        default = __file__.rsplit(os.sep,1)[0])                       
    args = parser.parse_args()
    
    if (args.exist_file == None) and (args.target_file == None):
        print "No files containing linkages to be made were provided."
        sys.exit()
    
    return args

pymol.finish_launching()

def get_protein_link_info(sel_text):
    '''Take a selection text and if it is a ASN, SER or THR residue return 
    filename for glycan template (n-linked for the former and variation on 
    o-linked for the latter), the residue type and appropriate link atom name.

    Args:
        sel_text (str):    Pymol selection text for the protein residue 
                           candidate for link to carbohydrate
                         
    Returns:
        dict: Contains link template PDB filename, link residue type and link atom 

        Dictionary is of the form:

        {
            'template': filename of appropriate glycan template,
            'link_res': link residue type (ASN, SER or THR )
            'link_atom': selection text for the link atom 
        }

        If the passed selection text does not find a ASN, SER or THR residue
        return empty dict 
         
    '''
    
    output = {}
    
    res_type = pymol.cmd.get_model(sel_text, 1).atom[0].resn
    
    if res_type in ['SER', 'THR', 'ASN']:

        # Load template glycan (o-linked to a SER or THR as appropriate)
        # The name of the resulting pymol pbject is stored in glycan_structure
        if res_type == 'SER':
            output['template'] = 'ser_o-link.pdb'
            output['link_res'] = 'SER'
            output['link_atom'] = 'OG'
        elif res_type == 'THR':
            output['template'] = 'thr_o-link.pdb'
            output['link_res'] = 'THR'
            output['link_atom'] = 'OG1'
        else:
            output['template'] = 'n-link.pdb'
            output['link_res'] = 'ASN'
            output['link_atom'] = 'ND2'    
    
    return output

def create_gly_link(pro_res_sel, pro_atom, gly_res_sel):
    
    pro_sel = pro_res_sel + ' and name ' + pro_atom
    gly_sel = gly_res_sel + ' and name C1'
    pymol.cmd.select('b1', pro_sel)
    pymol.cmd.select('b2', gly_sel)
    pymol.cmd.bond('b2','b1')

# --------------------
# Start of main script
# --------------------


# Setting to write CONECT records to PDB
#pymol.cmd.set("pdb_conect_all", "on")
# Setting to read CONECT records in PDBs
pymol.cmd.set("connect_mode", 1)

# Interpret command line arguments
args = parse_arguments()
# Use the file name to provide the name for the target structure object
target_structure = args.in_pdb.split('/')[-1].split('.')[0]

# Load targetPDB structure
pymol.cmd.load(args.in_pdb, target_structure)

# Define the selection used to superimpose glycosylated
# template with the target ASN, THR or SER
atom_mask = ' and not name N+C+O'

if args.target_file != None:
    # Read in YAML file containing the target residues
    f = open(args.target_file, 'r')
    target_res_list = yaml.load(f)
else:
    target_res_list = []
    
if args.exist_file != None:
    # Read YAML file containing glycan links supposedly present in the original structure
    f = open(args.exist_file, 'r')
    exist_links = yaml.load(f)
    for link in exist_links:
        
        pro_sel_text = target_structure + ' and chain ' + link['pro_chain'] + ' and resi ' + str(link['pro_res'])
        link_info = get_protein_link_info(pro_sel_text)
        
        if link_info:
            # Ensure that the glycosydic linkage is included
            gly_link_res = target_structure + ' and chain ' + link['gly_chain'] + ' and resi ' + str(link['gly_res'])
            create_gly_link(pro_sel_text, link_info['link_atom'], gly_link_res)

for residue in target_res_list:

    chain_id = residue['chain']
    res_no = residue['resid']
    str_res_no = str(res_no)

    # O-links only occur on SER and THR residues
    # N-links at ASN residues
    # Check to see that a valid residue type was specified
    
    target_residue = target_structure + " and chain " + chain_id + " and resi " + str_res_no
    link_info = get_protein_link_info(target_residue)

    if link_info:
        
        target_res_selection = pymol.cmd.get_model(target_residue + " and chain " + chain_id, 1)
        res_type = target_res_selection.atom[0].resn
        seg_id = target_res_selection.atom[0].segi
            
        glycan_structure = 'glycan'
        pymol.cmd.load(os.path.join(args.glycan_path,link_info['template']), glycan_structure)

        # Align the linked residue in the template with the target residue
        # using the atoms selected in atom_mask
        target_align_selection = target_residue + atom_mask
        linkResidue = glycan_structure + " and resn " + link_info['link_res']
        link_res_selection = linkResidue + atom_mask
        glycan_selection = glycan_structure + " and (not resn " + link_info['link_res'] + ")"
        pymol.cmd.align(link_res_selection, target_align_selection)

        # Remove the link residue in the template
        pymol.cmd.remove(linkResidue)

        # Alter the chain_id of the template glycan
        pymol.cmd.alter(glycan_structure, "chain = '" +  chain_id + "'")
        pymol.cmd.alter(glycan_structure, "segi = '" +  seg_id + "'")
        
        # Give the glycan residues a residue number after that of the existing atoms
        # of the chain containing the target link residue
        new_res_no = int(pymol.cmd.get_model(target_structure + " and chain " + chain_id, 1).atom[-1].resi) + 1
        first_glycan_res = int(pymol.cmd.get_model(glycan_selection,1).atom[0].resi)
        last_glycan_res = int(pymol.cmd.get_model(glycan_selection,1).atom[-1].resi)

        for old_res_no in range(first_glycan_res, last_glycan_res + 1):
            pymol.cmd.alter(glycan_selection + ' and resi ' + str(old_res_no), 'resi = ' + str(new_res_no))
            new_res_no += 1

        # Create a new structure combining the original target structure and the added glycan
        pymol.cmd.create('glycan_added', glycan_structure + ' or ' + target_structure)
        
        # Ensure that the glycan-protein linkage is incorporated in the new model
        new_pro_link = 'glycan_added and chain ' + chain_id + ' and resi ' + str_res_no
        new_gly_link = 'glycan_added and chain ' + chain_id + ' and resi ' + str(new_res_no)
        create_gly_link(new_pro_link, link_info['link_atom'], new_gly_link)
        
        
        # Delete both old structures
        pymol.cmd.delete(glycan_structure)
        pymol.cmd.delete(target_structure)
        
        # Rename the newly constructed structure to become the new target
        # for any subsequent glycan additions
        pymol.cmd.set_name('glycan_added',target_structure)
        
    else:
        print "Selection '" + target_residue + "' is not a ASN, SER or THR residue and has been ignored."

# Sort and save the glycosylated structure
pymol.cmd.sort(target_structure)

pymol.cmd.save(args.out_pdb, target_structure + ' and not hydro')

# Get out!
pymol.cmd.quit()
