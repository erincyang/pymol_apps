#!/usr/bin/python3

#For running on laptop
#pymol -r /Users/audreyoneill/digs/aoneill1/KingRotation/scripts/setup_view.py -- -w ../input/ECY_I3-30_trimer.pdb -i 66W_67G_0001.pdb

#For running on workstation in lab
#pymol -r /home/aoneill1/KingRotation/scripts/setup_view.py -- -w ../input/ECY_I3-30_trimer.pdb -i 66W_67G_0001.pdb

import pymol
from pymol import cmd
import subprocess
from subprocess import call
import argparse
import os

#When running script, first model is WT trimer and second is double mutant

def view_overall():
#    obj_list = cmd.get_object_list()
#    wt = obj_list[0]
#    mut = obj_list[1]
    global wt
    global mut
    
    sele_wt = wt.split("/")[-1][:-4] + " and chain A"
    sele_mut = mut[:-4] + " and chain A"
    
    cmd.super(str(sele_mut), str(sele_wt))
    return None

def view_mut():
 #   obj_list = cmd.get_object_list()
 #   wt = obj_list[0]
 #   mut = obj_list[1]
    global wt
    global mut
    
    resi1 = str(mut).split("_")[0][:-1]
    resi2 = str(mut).split("_")[1][:-1]
    
    sele_string = "resi " + resi1 + "+" + resi2
    cmd.center(str(sele_string))
    cmd.show("sticks")
    cmd.remove("hydro")
    cmd.zoom("center", "10")
    cmd.spectrum("b", "blue_red", str(mut[:-4]))
    cmd.show("spheres", "resi " + str(resi1))
    cmd.show("spheres", "resi " + str(resi2))
    return None

#cmd.extend('view_overall', view_overall)
#cmd.extend('view_mut', view_mut)

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--pdb_to_read', type=str, default="", help="PDB file to get per residue energies from")
parser.add_argument('-w', '--pdb_of_wt', type=str, default="", help="PDB file of WT")
args = parser.parse_args()

wt = args.pdb_of_wt
cmd.load(str(wt))

#For running on laptop
subprocess.call(["/Users/audreyoneill/digs/aoneill1/KingRotation/scripts/e_to_b-value.sh", str(args.pdb_to_read).strip()])    

#For running on lab workstation
#subprocess.call(["/home/aoneill1/KingRotation/scripts/e_to_b-value.sh", str(args.pdb_to_read).strip()])

#Name of new PDB file generated from bash script    
mut = str(args.pdb_to_read).strip()[:-4] + "_energies.pdb"

#Load now PDB file with per residue energies in B values column
cmd.load(str(mut))
   
#Remove new PDB file that contained per residue energies in B values column 
os.remove(mut)
    
cmd.extend('view_overall', view_overall)
cmd.extend('view_mut', view_mut)

view_overall()
view_mut()
