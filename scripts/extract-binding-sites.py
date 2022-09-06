#!/usr/bin/env python
import argparse
import numpy as np
from Bio.PDB import MMCIFParser
import logging
import sys
from scipy.spatial import KDTree
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("extract binding sites")

def main():
    """
    extract ligand binding site from cif format structure
    for the given ligand, consider all of its atom
    if nearest atom in RNA < 6 Ångström, the residue which the atom belongs to is consider as a binding site
    """
    parser = argparse.ArgumentParser(description='extract active site')
    parser.add_argument('--input','-i',type=str,required=True,help="Input RNA-ligand complex structure in cif format")
    parser.add_argument('--ligand-id','-l',type=str,required=True,help="ID of the ligand to consider")
    parser.add_argument('--output','-o',type=str,required=True,help="Where to save binding sites")
    parser.add_argument('-k',type=int,default=1,help="kNN to retrieve")
    parser.add_argument('--distance','-d',type=float,default=6,help="Distancec cutoff to use in Ångström")
    args = parser.parse_args()     
    
    logger.info(f"Load complex structure from {args.input}")
    pdb_id = ".".join(args.input.split("/")[-1].split(".")[:-1])
    parser = MMCIFParser()
    structure = parser.get_structure(pdb_id,args.input)
    chains = structure[0].child_dict
    logger.info(f"{len(chains)} chains are loaded .")
    atom_ids = []
    coordinates = []
    residues = {}
    ligands = []
    for chain_id in chains:
        logger.info(f"Processing chain {chain_id} ...")
        residues[chain_id] = ""
        position = 0
        for r in chains[chain_id].get_residues():
            resname = r.get_resname()
            if resname in "ACGU":
                residues[chain_id] += resname
                for atom in r.get_atoms():
                    atom_ids.append((chain_id, position, resname, atom.name))
                    coordinates.append(atom.coord)
                position += 1
            elif r.get_resname() == args.ligand_id:
                logger.info(f"Ligand {args.ligand_id} detected in chain {chain_id}.")
                ligands.append((chain_id, r))
            else:
                #logger.info(f"skip residue {resname}")
                continue
    
    if len(ligands) == 0:
        logger.error("The specified ligand is not present in input complex structure, please check .")
        sys.exit(1)
    atom_ids = np.array(atom_ids)
    coordinates = np.array(coordinates)
    logger.info("Build kd-tree of atom coordinate for efficient 1NN query ...")
    kdtree = KDTree(data=coordinates)
   
    fout = open(args.output,"w")
    for i, (chain_id, ligand) in enumerate(ligands):
        logger.info(f"Processing ligand instance {i} ...")
        ligand_atom_coordinates = []
        for atom in ligand.get_atoms():
            ligand_atom_coordinates.append(atom.coord)
        logger.info("Perform 1NN query ...")
        distances, indices = kdtree.query(np.array(ligand_atom_coordinates),k=args.k)
        # chain position
        logger.info(f"Consider Ångström < {args.distance} for binding site determination .")
        binding_sites = np.unique(atom_ids[np.unique(indices[distances<args.distance].reshape(-1))][:,[0,1]],axis=0)     
        logger.info("Saving binding sites ...")
        for chain_id in np.unique(binding_sites[:,0]):
            binding_sites_by_chain = binding_sites[binding_sites[:,0] == chain_id,:]
            active_sites = sorted(binding_sites_by_chain[:,1].astype(int))
            print(">" + chain_id + " " + str(i), file=fout)
            print(">" + chain_id + " " + str(i))
            print(residues[chain_id],file=fout)
            print(residues[chain_id])
            sequence = residues[chain_id]
            mask = np.full(len(sequence),False)
            mask[active_sites] = True
            print("".join(mask.astype(int).astype(str)),file=fout)
            print("".join(mask.astype(int).astype(str)))
    fout.close()
    logger.info("All done .")
if __name__ == "__main__":
    main()
