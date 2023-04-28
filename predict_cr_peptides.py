#       Sgourakis Lab
#   Author: Sagar Gupta
#   Date: April 28, 2023
#   Email: sagarg@sas.upenn.edu

# import required libraries
from collections import defaultdict
import csv
import math
import os
import platform
from pymol import cmd
import subprocess
from tqdm import tqdm


def get_pdb_dict(fastafile):
    '''
        Reads FASTA file
    '''
    
    fasta_file_handler = open(fastafile, "r")
    counter = 1
    raw_dict = {}
    pdbid = ""
    for line in fasta_file_handler:
        if ">" in line:
            line = line.rstrip().lstrip('>')
            pdbid = line.split('|')[0]
            allele = line.split('|')[2]
            if counter % 2 == 0:
                raw_dict[pdbid] = ""
            counter += 1
        else:
            line = line.rstrip()
            raw_dict[pdbid] = line+"_"+allele
            
    fasta_file_handler.close()
    pdb_dict = {}
    for key, value in raw_dict.items():
        if value != "":
            pdb_dict[key] = value

    return pdb_dict

def run_command(command, ignore=False):
    '''
        Helps run command line operations
    '''

    cmd = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = cmd.communicate()

    if cmd.returncode != 0:
        if ignore:
            return False
        print("--- FAIL ---\n")
        print(command)
        if output[0] != None:
            print(str(output[0], 'utf-8'))
        if output[1] != None:
            print(str(output[1], 'utf-8'))
        print("\n--- FAIL ---")
        return False
    else:
        if output[0] != None:
            output = output[0].decode('utf-8')
        return output

def distance(angle1, angle2):
    '''
        Determines the dihedral difference between two angles
        Assumes inputted angles are in degrees
    '''

    distance = 2.0 * (1.0 - math.cos(math.radians(angle1) - math.radians(angle2)))
    return distance

def make_resfile(directory, end_peptide_seq, structure_pepchain):
    '''
        Write resfile which mutates the template peptide sequence
    '''

    with open(f"{directory}/{end_peptide_seq}.resfile", "w") as resfile:

        # preserve the input rotamer (do not pack at all) - applies to everything
        resfile.write("NATRO\n")
        resfile.write("start\n")

        for pos, end_aa in enumerate(end_peptide_seq):
            resfile.write(f"{pos+1} {structure_pepchain} PIKAA {end_aa}\n")

    return f"{directory}/{end_peptide_seq}.resfile"

def main():

    '''
        THE LINE BELOW NEEDS TO BE CHANGED BASED ON THE USER
    '''
    ROSETTA_INSTALL_DIR = "/path_to_dir/rosetta/"

    os_system = platform.system()
    extension = ''
    if os_system == 'Linux':
        extension = 'linuxgccrelease'
    elif os_system == 'Darwin':
        extension = 'macosclangrelease'


    '''
    Crystal structure is WT A*24:02 (chain A) bound to PHOX2B (chain B) in complex with 10LH, heavy chain (chain H) and light chain (chain L)
    '''
    crystal_structure = f"PHOX2B_A2402_10LH.pdb"

    # pHLA-I structural database taken from https://hla3db.research.chop.edu with the 
    # following query: peptide length 9, anchor class ∆7
    # as of April 28th, 2023 there were 315 structures
    # NOTE: The manuscript uses a cutoff date of April 29th, 2022
    pdb_dict = get_pdb_dict("MHC_pdbs/MHCs.fasta")
    
    '''
    Replaces the PHOX2B peptide backbone in the A2402 and 10LH complex
    '''
    workdir1 = "template_graft/"
    os.system(f"mkdir -p {workdir1}")
    print(f"Grafting template peptide backbone via HLA alignment (n = {len(pdb_dict)}) (runtime: ~20 seconds)")
    for pdbid in tqdm(pdb_dict.keys()):

        cmd.load(crystal_structure, "10LH")
        cmd.load(f"MHC_pdbs/{pdbid}_reordered.pdb", "template")

        cmd.super("template and chain A and name ca", "10LH and chain A and name ca", cycles=0) # align HLA Ca atoms
        cmd.remove("template and chain A") # removes template HLA
        cmd.remove("10LH and chain B") # removes 10LH peptide (PHOX2B)
        cmd.save(f"{workdir1}/{pdbid}_A2402_10LH.pdb")
        cmd.delete("all")

    '''
    Replaces the new backbone sequence with PHOX2B peptide sequence and performs \extensive\ rotamer sampling
    '''
    print("Introducing hotspot residues and conducting extensive rotamer sampling (runtime: ~35 minutes)")
    workdir2 = "template_graft_hotspot/"
    os.system(f"mkdir -p {workdir2}")
    for pdbid in tqdm(pdb_dict.keys()):

        resfile_name = make_resfile(workdir2, 'AAAAARAAA', "B")  # change the peptide sequence to 'alanine' PHOX2B peptide sequence
        template_structure_new_bb = f"{workdir1}/{pdbid}_A2402_10LH.pdb"
        
        # since R6 is not aromatic, we do not increase rotamer sampling for aromatic residues
        os.system(f"{ROSETTA_INSTALL_DIR}/main/source/bin/fixbb.{extension} -database {ROSETTA_INSTALL_DIR}/main/database -s {template_structure_new_bb} -out:path:all {workdir2} -resfile {resfile_name} -suffix _PHOX2B -ex1 -ex2 -ex3 -ex4 -overwrite")
        
    '''
    Gets relevant distances/angles between R6 and D232 in the fixed backbone models
    '''
    print("Selecting template backbones with ideal target side chain geometry (runtime: ~15 seconds)")
    output_dir = "output"
    os.system(f"mkdir -p {output_dir}")
    with open(f"{output_dir}/template_graft_hotspot_dist_angle.csv", "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['pdbid', 'R6-CZ_D104-DG', "R6-CZ_R6-NH1_D104-OD2_D104-CG", "R6-CZ_R6-NH1_D104-OD2_D104-CG"])
        for pdbid in tqdm(pdb_dict.keys()):

            fixbb_model = f"{workdir2}/{pdbid}_A2402_10LH_PHOX2B_0001.pdb"
            cmd.load(fixbb_model, "struc")

            cz_cg_distance = cmd.get_distance("struc//B/ARG`6/CZ", "struc//H/ASP`104/CG")

            nh1_od2_dihedral_angle = cmd.get_dihedral("struc//B/ARG`6/CZ", "struc//B/ARG`6/NH1", "struc//H/ASP`104/OD2", "struc//H/ASP`104/CG")
            nh2_od1_dihedral_angle = cmd.get_dihedral("struc//B/ARG`6/CZ", "struc//B/ARG`6/NH2", "struc//H/ASP`104/OD1", "struc//H/ASP`104/CG")

            cmd.delete("all")

            writer.writerow([pdbid, cz_cg_distance, nh1_od2_dihedral_angle, nh2_od1_dihedral_angle])

    template_graft_hotspot_dist_angle = defaultdict(list)
    with open(f"{output_dir}/template_graft_hotspot_dist_angle.csv", "r") as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for line in reader:
            template_graft_hotspot_dist_angle[line[0]] = [float(x) for x in line[1:]]
    
    '''
    Finds fixbb models that fit the criteria for a R6/D232 interaction
    Criteria:
        1. Dihedral difference (D-score) between (CZ, NH1, and CG, OD1) + (CZ, NH2, and CG, OD2) < 3.5 and each individual D-score < 2.3
        2. Distance between CG and CZ is < 5 Å

    For reference, in 10LH:PHOX2B/A24, the distance between CG and CZ is 4.1 Å
    '''
    cmd.load(crystal_structure, "10LH")
    REF_NH1_OD2 = cmd.get_dihedral("10LH//B/ARG`6/CZ", "/10LH//B/ARG`6/NH1", "10LH//H/ASP`104/OD2", "10LH//H/ASP`104/CG")
    REF_NH2_OD1 = cmd.get_dihedral("10LH//B/ARG`6/CZ", "/10LH//B/ARG`6/NH2", "10LH//H/ASP`104/OD1", "10LH//H/ASP`104/CG")
    cmd.delete("all")

    valid_models = []
    for pdbid, metrics in template_graft_hotspot_dist_angle.items():
        cz_cg_distance = metrics[0]
        nh1_od2_dihedral_angle = metrics[1]
        nh2_od1_dihedral_angle = metrics[2]

        if cz_cg_distance > 5:
            continue

        dscore1 = distance(nh1_od2_dihedral_angle, REF_NH1_OD2)
        dscore2 = distance(nh2_od1_dihedral_angle, REF_NH2_OD1)
        combined_dscore = dscore1 + dscore2

        if combined_dscore < 4 and dscore1 < 2.3 and dscore2 < 2.3:
            valid_models.append(pdbid)
          
    print(f"There are {len(valid_models)} valid models")
    with open(f"{output_dir}/valid_models.txt", "w") as txtfile:
        for pdbid in valid_models:
            txtfile.write(f"{pdbid}\n")
    
    '''
        Generate heteroclitic peptides with target hotspot residues
    '''
    new_pep_seq = {}
    for pdbid in valid_models:
        peptide = pdb_dict[pdbid].split("_")[0] # template peptide sequence

        new_peptide = [""]*9
        for index, resi in enumerate(peptide):
            pos = index+1

            # establish A*24 anchors
            if pos == 2:
                new_peptide[index] = 'Y'
            elif pos == 9:
                new_peptide[index] = 'F'

            # establish PHOX2B hotspot residue
            elif pos == 6:
                new_peptide[index] = "R"

            # reduce positive charge around R6
            elif (pos == 7) and (resi == "R" or resi == "H" or resi == "K"):
                new_peptide[index] = "G"
            elif (pos == 4 or pos == 5) and (resi == "R" or resi == "H" or resi == "K"):
                new_peptide[index] = "I"

            else:
                new_peptide[index] = resi

        new_peptide = "".join(new_peptide)
        new_pep_seq[pdbid] = new_peptide
    
    unique_pep_seq = list(set(new_pep_seq.values()))

    print(f"There are {len(unique_pep_seq)} unique peptide sequences")

    with open(f"{output_dir}/potential_cr_pep_seq.txt", "w") as txtfile:
        for pep_seq in unique_pep_seq:
            txtfile.write(f"{pep_seq}\n")
    
        
if __name__ == "__main__":

    main()
