# -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# previous line and header box below doc string from Will

'''
Pymol script to superimpose residues and other molecules and generate constraint files

Use it in pymol:

	type 'run PATH_TO/enzdes_pymol_gen.py'
	or put that in your .pymolrc


	enter editing mode and click on three atoms


	type 'name_for_res_1 = res(1,"pk3","pk1","pk2")'
		or just 'res2 = res(2) which means 'res2 = res(2,"pk1","pk2","pk3")'


	this makes a res object, a tiny residue class that only considers three atoms

		THE ORDER of the pk atoms given when making the res object defines atoms 1,2,3 in cst calculations and superimpositions
		the order you click doesn't matter unless you use the default args like res77 = res(77)

	
	Do it again with a second three atoms


	Now either:
		type 'cst(res1,res2,'path_to/your_file.cst') to generate cst file template with measurements in place

		this doesn't do anything about symmetric equilivency

		file arguement is optional, and the path starts in your home directory (or where-ever your pymol dumps fetched pdbs)

	Or:
		type 'suppose(res1,res2)' to put res1 on top of res2

	To repeat for another cst you must first run:
	cst_reset()
	then you must repeat the whole process for the next residue pair (even if one of the residues is the same; sorry)

Contact Harley at pylesharley@gmail.com with questions  

'''

# Adapted from Will Sheffler's header from pymol_util
##################################
import sys,os,inspect
newpath = os.path.dirname(inspect.getfile(inspect.currentframe())) # script directory
if not newpath in sys.path:
	sys.path.append(newpath)

import string, math, re, gzip, itertools, glob, sets
from random import randrange
from math import sqrt
# import xyzMath as xyz
# from xyzMath import Ux,Uy,Uz,Imat
# from functools import partial
import cProfile

try:
	import pymol
	from pymol import cmd,CmdException,cgo,util
	from chempy import cpv
	def inpymol(): return True
except ImportError as e:
	print "can't load pymol, Mocking it for testing/doc"
	from minimock import Mock
	pymol = Mock("pymol")
	cmd = Mock("cpv")
	cmd = Mock("cmd")
	cgo = Mock("cgo")
	cgo = Mock("util")
	cmd.extend = (lambda x,y: x)
	def inpymol(): return False

# End of elements from Will's pymol_utils.py
############################################

# importing modules
from pymol import stored
import numpy as np
import subprocess
import time
import os

OneLetter = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

class res:
	""" three atom residue constraints for generating cst files """
	def __init__(self, id_number=0, atom1 = "pk1", atom2 = "pk2", atom3 = "pk3"):

		self.a1_name = 'residue_%d_atom_1'%id_number
		self.a2_name = 'residue_%d_atom_2'%id_number
		self.a3_name = 'residue_%d_atom_3'%id_number
		print atom1,'->',self.a1_name
		print atom2,'->',self.a2_name
		print atom3,'->',self.a3_name

		# pymol cmd function does the rest
		cmd.set_name(atom1,self.a1_name)
		cmd.set_name(atom2,self.a2_name)
		cmd.set_name(atom3,self.a3_name)

		self.atom1 = cmd.get_model(self.a1_name, 1).atom[0]
		self.atom2 = cmd.get_model(self.a2_name, 1).atom[0]
		self.atom3 = cmd.get_model(self.a3_name, 1).atom[0]

def suppose(res1,res2):
	''' Aligns residue 1 to residue 2. Only considers each residue's three defined atoms '''
	cmd.pair_fit(res1.a1_name,res2.a1_name,res1.a2_name,res2.a2_name,res1.a3_name,res2.a3_name)

def res_reset():
	cmd.delete('res*atom*')

def cst_reset():
	cmd.delete('angle_*')
	cmd.delete('distance_*')
	cmd.delete('torsion_*')
	cmd.delete('res*atom*')

def resset(): # no underscores 
	res_reset()

def recst():
 	cst_reset()

def cst(res1,res2,outfile=0):
	''' generates enzdes style constraint '''	
	cst_header = '''CST::BEGIN
  TEMPLATE::   ATOM_MAP: 1 atom_name: %s %s %s
  TEMPLATE::   ATOM_MAP: 1 residue3: %s

  TEMPLATE::   ATOM_MAP: 2 atom_name: %s %s %s ,
  TEMPLATE::   ATOM_MAP: 2 residue1:  %s'''%(	res1.atom1.name,res1.atom2.name,res1.atom3.name, 
  												res1.atom1.resn,
												res2.atom1.name,res2.atom2.name,res2.atom3.name, 
												OneLetter[res2.atom1.resn])
	
	distanceAB = cmd.distance('distance_AB', res1.a1_name, res2.a1_name)
	cst_distance = '''
  CONSTRAINT:: distanceAB:   %.2f   0.20   80.0  0     0'''%(distanceAB)

	angleA = cmd.angle('angle_A',	res1.a2_name, res1.a1_name, res2.a1_name)
	angleB = cmd.angle('angle_B',	res2.a2_name, res2.a1_name, res1.a1_name)
	cst_angles = '''
  CONSTRAINT::    angle_A:   %.2f  10.0   10.0  360   0
  CONSTRAINT::    angle_B:   %.2f  10.0   10.0  360   0'''%(angleA,angleB)

	torsionA = cmd.dihedral('torsion_A',	res1.a3_name, res1.a2_name, res1.a1_name, res2.a1_name)
	torsionAB = cmd.dihedral('torsion_AB',	res1.a2_name, res1.a1_name, res2.a1_name, res2.a2_name)
	torsionB = cmd.dihedral('torsion_B',	res2.a3_name, res2.a2_name, res2.a1_name, res1.a1_name)
	
	if torsionA < 0: torsionA += 360.0
	if torsionAB < 0: torsionAB += 360.0
	if torsionB < 0: torsionB += 360.0
	
	cst_torsions = '''
  CONSTRAINT::  torsion_A:   %.2f   30.0   10.0  360   0
  CONSTRAINT:: torsion_AB:   %.2f   30.0   10.0  360   0
  CONSTRAINT::  torsion_B:   %.2f   30.0   10.0  360   0'''%(torsionA, torsionAB, torsionB)

	full_cst = '''%s\n%s%s%s
  ALGORITHM_INFO:: match
   MAX_DUNBRACK_ENERGY 4.0 
   IGNORE_UPSTREAM_PROTON_CHI
  ALGORITHM_INFO::END

CST::END''' %(cst_header, cst_distance, cst_angles, cst_torsions)
	
	# prints cst in pymol window
	print full_cst

	# prints cst to path/file in filename
	if outfile:
		with open(outfile,'a') as openoutfile:
			print>>openoutfile, full_cst
