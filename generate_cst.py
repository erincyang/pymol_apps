#! /usr/bin/env python

InfoString = ''' 
This script is to generate ATOMPAIR constraints for Rosetta,
 by default for all downstreams and oxygens
                    within 3 angstroms 
                    on non-neighbor residues 
                    within an input pose. 

'''

# uncomment just next line and copy block in multiline string for ipython mode
# '''

from multiprocessing import Process
from scipy import spatial
from Bio import PDB
import numpy as np
import subprocess
import argparse
import sys
import os
import re

if '-h' not in sys.argv:
  import solenoid_tools
  import rosetta
  rosetta.init(extra_options = "-mute basic -mute core -mute protocols")

# '''

# sys.argv.extend(['-pdbs', '1EZG.pdb', '-out', './' ])

ThreeToOne = {'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','MET':'M','PRO':'P','PHE':'F','TRP':'W','SER':'S','THR':'T','ASN':'N','GLN':'Q','TYR':'Y','CYS':'C','CYD':'C','LYS':'K','ARG':'R','HIS':'H','ASP':'D','GLU':'E','STO':'*','UNK':'U'}
ChainAlphabetIndices = {'A':1, 'B':2, 'C':3, 'D':4, 'E':5, 'F':6, 'G':7, 'H':8, 'I':9, 'J':10, 'K':11, 'L':12, 'M':13, 'N':14, 'O':15, 'P':16, 'Q':17, 'R':18, 'S':19, 'T':20, 'U':21, 'V':22, 'W':23, 'X':24, 'Y':25, 'Z':26 }

BackboneAtomList = ['C', 'CA', 'N', 'O']

def pymol_commands(Pdb, Repeat, ReportedRepeatCount):
  return 'fetch %s\tselect rep%d, resi %s'%( Pdb, ReportedRepeatCount, '+'.join([str(Res) for Res in Repeat]) )

def get_pose_constraints(Pose, MaxDist, MinPositionSeperation, UpstreamGrep, DownstreamGrep, Weight, NeedHydrogen=True):

    # for making full atom kd tree
    ResAtmCoordLists = []
    # for translating from kd tree index to ( residue, atom ) coord
    ResAtmRecordLists = []

    # loop through all residue numbers
    for Res in range(1, Pose.n_residue() + 1):
      # remade for each residue
      AtmRecordList = []
      AtmCoordList = []
      # loop through residue's atom numbers
      for Atm in range(1, Pose.residue(Res).natoms() + 1):
        # add (residue, atom) coord to residue's list
        AtmRecordList.append((Res, Atm))
        # add atom xyz coord to residue's list
        AtmCoordList.append( np.array(list(Pose.residue(Res).atom(Atm).xyz())) )
      
      # add residue's lists to respective global lists
      ResAtmCoordLists.extend(AtmCoordList)
      ResAtmRecordLists.extend(AtmRecordList)

    ResidueAtomArray = np.array( ResAtmCoordLists )
    ResidueAtomKDTree = spatial.KDTree( ResidueAtomArray )

    ResidueAtomNeighbors = ResidueAtomKDTree.query_ball_point( ResidueAtomArray, MaxDist )
    # ResidueAtomNearNeighbors = ResidueAtomKDTree.query_ball_point( ResidueAtomArray, 2.0 )
    ResidueAtomHydrogens = ResidueAtomKDTree.query_ball_point( ResidueAtomArray, 1.1 )

    # holds constraints before printing
    AllConstraints = [] 
    # holds sorted cst
    AllBackboneBackboneCst = []
    AllBackboneSidechainCst = []
    AllSidechainSidechainCst = []

    # All contacts are from upstream to downstream residues to avoid double counting
    Upstream = []
    for UpIndex, UpXyzCoords in enumerate(ResAtmCoordLists):
      UpRes, UpAtm = ResAtmRecordLists[UpIndex]

      # # loop through residues storing info on oxygens
      # for UpRes in range( 1, Pose.n_residue() + 1 ):
      #   # loop through atoms
      #   for UpAtm in range( 1, Pose.residue(UpRes).natoms() + 1 ):
      UpName = Pose.residue(UpRes).atom_name(UpAtm).replace(' ', '')

      # skip virtual residues
      if Pose.residue(UpRes).is_virtual(UpAtm):
        continue

      #                                this guy 
      #                                 /
      # checks upstream name           V
      if re.match(UpstreamGrep, UpName ): 
        # print '\n'*2
        # print 'UpRes, UpName', UpRes, UpName

        # get neighbors of upstream residues
        NeighborsOfUpstream = ResidueAtomNeighbors[UpIndex]
        
        # prep for loop
        Downstreams = []

        Constraints = []
        BackboneBackboneCst = []
        BackboneSidechainCst = []
        SidechainSidechainCst = []

        # ArbitrayOrderOfAtomNames = {}
        for DownIndex in NeighborsOfUpstream:
          # name presumes downstream, checks with if imediately below
          DownRes, DownAtm = ResAtmRecordLists[DownIndex]

          # checks that downstream residue is dowstream of upstream and passes min primary sequence spacing
          if DownRes - UpRes >= MinPositionSeperation:
            DownName = Pose.residue(DownRes).atom_name(DownAtm).replace(' ', '')
            
            # skip if same atom
            if UpRes == DownRes:
              if UpName == DownName:
                continue

            # skip virtual residues
            if Pose.residue(DownRes).is_virtual(DownAtm):
              continue

            # checks downstream name
            if re.match( DownstreamGrep, DownName ):
              # print 'DownRes, DownName', DownRes, DownName

              PotentialUpstreamHydrogens = ResidueAtomHydrogens[UpIndex]
              UpstreamHydrogens = []
              # print 'PotentialUpstreamHydrogens', PotentialUpstreamHydrogens
              for UpH_I in PotentialUpstreamHydrogens:
                UpH_Res, UpH_Atm = ResAtmRecordLists[UpH_I]
                UpH_Name  = Pose.residue(UpH_Res).atom_name(UpH_Atm).replace(' ', '')
                # print 'UpH_Name', UpH_Name
                if 'H' in UpH_Name:
                  UpstreamHydrogens.append((UpH_Res, UpH_Atm, UpH_Name))
                # print 'UpstreamHydrogens', UpstreamHydrogens

              PotentialDownstreamHydrogens = ResidueAtomHydrogens[DownIndex]
              DownstreamHydrogens = []
              # print 'PotentialDownstreamHydrogens', PotentialDownstreamHydrogens
              for DownH_I in PotentialDownstreamHydrogens:
                DownH_Res, DownH_Atm = ResAtmRecordLists[DownH_I]
                DownH_Name = Pose.residue(DownH_Res).atom_name(DownH_Atm).replace(' ', '')
                # print 'DownH_Name', DownH_Name
                if 'H' in DownH_Name:
                  DownstreamHydrogens.append((DownH_Res, DownH_Atm, DownH_Name))
                # print 'DownstreamHydrogens', DownstreamHydrogens

              # check their is at least one hydrogen in system before adding constraint
              if len(UpstreamHydrogens) or len(DownstreamHydrogens) or NeedHydrogen == False:

                # if/elif/else seperates cst into groups for
                # (BBBB) just backbone-backbone interactions from 
                # (BBSC) just backbone-sidechain interactions from
                # (SCSC) just sidechain-sidechain interactinos
                # 
                
                if UpName in BackboneAtomList and DownName in BackboneAtomList:
                  BBBB = 1
                  BBSC = 0
                  SCSC = 0
                elif not UpName in BackboneAtomList and not DownName in BackboneAtomList:
                  BBBB = 0
                  BBSC = 0
                  SCSC = 1
                else:
                  BBBB = 0
                  BBSC = 1
                  SCSC = 0

                # print 'UpName', UpName
                # print 'DownName', DownName
                # print 

                # print 'found downstream neighbor %s'%DownName
                DownXyzCoords = np.array( list(Pose.residue(DownRes).atom(DownAtm).xyz()) )
                # print 'DownRes, DownName', DownRes, DownName
                # print 'DownXyzCoords', DownXyzCoords

                # ## Get neighbors for angles and torsions to use with AtomPairs

                SelectUpNeighbors = []
                # iterates through upstream atom neighbors for references for angle
                for UpNeighborIndex in NeighborsOfUpstream:
                  UpNeighborRes, UpNeighborAtm = ResAtmRecordLists[UpNeighborIndex]
                  UpNeighborName = Pose.residue(UpNeighborRes).atom_name(UpNeighborAtm).replace(' ', '')

                  # keep looking if neighbor is hyrdogen
                  if 'H' in UpNeighborName:
                    continue                

                  # skip virtual residues
                  if Pose.residue(UpNeighborRes).is_virtual(UpNeighborAtm):
                    continue

                  # keep looking if neighbor is self
                  if UpNeighborName == UpName and UpNeighborRes == UpRes:
                    continue
                  # keep looking if neighbor is downstream residue again
                  if UpNeighborName == DownName and UpNeighborRes == DownRes:
                    continue
                  UpNeighborCoords = ResAtmCoordLists[UpNeighborIndex]
                  DistanceToNeighbor = solenoid_tools.vector_magnitude( UpXyzCoords - UpNeighborCoords )
                  SelectUpNeighbors.append( (DistanceToNeighbor, UpNeighborName, UpNeighborRes, UpNeighborCoords) )

                # sort by distance to atom, nearest first
                SelectUpNeighbors.sort()                
                UpNeighbor1Tuple = SelectUpNeighbors[0]
                UpNeighbor2Tuple = SelectUpNeighbors[1]
                # print '\n'*2
                # print 'UpRes, UpName', UpRes, UpName
                # print 'UpstreamHydrogens', UpstreamHydrogens
                # print 'SelectUpNeighbors', SelectUpNeighbors

                 # get neighbors of upstream residues
                NeighborsOfDownstream = ResidueAtomNeighbors[DownIndex]
                SelectDownNeighbors = []
                # iterates through upstream atom neighbors for references for angle
                for DownNeighborIndex in NeighborsOfDownstream:
                  DownNeighborRes, DownNeighborAtm = ResAtmRecordLists[DownNeighborIndex]
                  DownNeighborName = Pose.residue(DownNeighborRes).atom_name(DownNeighborAtm).replace(' ', '')

                  # keep looking if neighbor is hyrdogen
                  if 'H' in DownNeighborName:
                    continue                

                  # skip virtual residues
                  if Pose.residue(DownNeighborRes).is_virtual(DownNeighborAtm):
                    continue

                  # keep looking if neighbor is self
                  if DownNeighborName == DownName and DownNeighborRes == DownRes:
                    continue
                  # keep looking if neighbor is upstream residue
                  if DownNeighborName == UpName and DownNeighborRes == UpRes:
                    continue

                  DownNeighborCoords = ResAtmCoordLists[DownNeighborIndex]
                  DistanceToNeighbor = solenoid_tools.vector_magnitude( DownXyzCoords - DownNeighborCoords )
                  SelectDownNeighbors.append( (DistanceToNeighbor, DownNeighborName, DownNeighborRes, DownNeighborCoords) )

                # sort by distance to atom, nearest first
                SelectDownNeighbors.sort()
                DownNeighbor1Tuple = SelectDownNeighbors[0]
                DownNeighbor2Tuple = SelectDownNeighbors[1]
                # print 'DownRes, DownName', DownRes, DownName
                # print 'DownstreamHydrogens', DownstreamHydrogens
                # print 'SelectDownNeighbors', SelectDownNeighbors

                Distance = solenoid_tools.vector_magnitude(DownXyzCoords - UpXyzCoords)
                DistanceCst = 'AtomPair %s %d %s %d SCALARWEIGHTEDFUNC %f HARMONIC %.2f 1.0' %( UpName, UpRes, DownName, DownRes, Weight, Distance )

                # Use Biopython for angle and dihedral calculations
                # here 'Vec' means PDB.Vector of atom's xyz coord
                UpstreamVec = PDB.Vector(UpXyzCoords)
                DownstreamVec = PDB.Vector(DownXyzCoords)
                
                UpNeighbor1Vec = PDB.Vector(UpNeighbor1Tuple[3])
                UpNeighbor2Vec = PDB.Vector(UpNeighbor2Tuple[3])
                DownNeighbor1Vec = PDB.Vector(DownNeighbor1Tuple[3])
                DownNeighbor2Vec = PDB.Vector(DownNeighbor2Tuple[3])

                Angle1 = PDB.calc_angle(UpNeighbor1Vec, UpstreamVec, DownstreamVec)
                AngleCst1 = 'Angle %s %d %s %d %s %d SCALARWEIGHTEDFUNC %f CIRCULARHARMONIC %.2f 0.5' %( UpNeighbor1Tuple[1], UpNeighbor1Tuple[2], UpName, UpRes, DownName, DownRes, Weight, Angle1 )
                Angle2 = PDB.calc_angle(UpstreamVec, DownstreamVec, DownNeighbor1Vec)
                AngleCst2 = 'Angle %s %d %s %d %s %d SCALARWEIGHTEDFUNC %f CIRCULARHARMONIC %.2f 0.5' %( UpName, UpRes, DownName, DownRes, DownNeighbor1Tuple[1], DownNeighbor1Tuple[2], Weight, Angle2 )

                Torsion1 = PDB.calc_dihedral(UpNeighbor2Vec, UpNeighbor1Vec, UpstreamVec, DownstreamVec)
                TorsionCst1 = 'Dihedral %s %d %s %d %s %d %s %d SCALARWEIGHTEDFUNC %f CIRCULARHARMONIC %.2f 0.5' %( UpNeighbor2Tuple[1], UpNeighbor2Tuple[2], UpNeighbor1Tuple[1], UpNeighbor1Tuple[2], UpName, UpRes, DownName, DownRes, Weight, Torsion1 )
                Torsion2 = PDB.calc_dihedral(UpNeighbor1Vec, UpstreamVec, DownstreamVec, DownNeighbor1Vec)
                TorsionCst2 = 'Dihedral %s %d %s %d %s %d %s %d SCALARWEIGHTEDFUNC %f CIRCULARHARMONIC %.2f 0.5' %( UpNeighbor1Tuple[1], UpNeighbor1Tuple[2], UpName, UpRes, DownName, DownRes, DownNeighbor1Tuple[1], DownNeighbor1Tuple[2], Weight, Torsion2 )
                Torsion3 = PDB.calc_dihedral(UpstreamVec, DownstreamVec, DownNeighbor1Vec, DownNeighbor2Vec)
                TorsionCst3 = 'Dihedral %s %d %s %d %s %d %s %d SCALARWEIGHTEDFUNC %f CIRCULARHARMONIC %.2f 0.5' %( UpName, UpRes, DownName, DownRes, DownNeighbor1Tuple[1], DownNeighbor1Tuple[2], DownNeighbor2Tuple[1], DownNeighbor2Tuple[2], Weight, Torsion3 )

                # adds constraint to running lists of constraints
                Constraints.extend( [DistanceCst, AngleCst1, AngleCst2, TorsionCst1, TorsionCst2, TorsionCst3] )
                if BBBB: BackboneBackboneCst.extend( [DistanceCst, AngleCst1, AngleCst2, TorsionCst1, TorsionCst2, TorsionCst3] )
                if BBSC: BackboneSidechainCst.extend( [DistanceCst, AngleCst1, AngleCst2, TorsionCst1, TorsionCst2, TorsionCst3] )
                if SCSC: SidechainSidechainCst.extend( [DistanceCst, AngleCst1, AngleCst2, TorsionCst1, TorsionCst2, TorsionCst3] )

              # else:
              #   print 'No hydrogen!'
              #   sys.exit()

        AllConstraints.extend(Constraints)
        AllBackboneBackboneCst.extend(BackboneBackboneCst)
        AllBackboneSidechainCst.extend(BackboneSidechainCst)
        AllSidechainSidechainCst.extend(SidechainSidechainCst)

    SortedConstraints = (AllBackboneBackboneCst, AllBackboneSidechainCst, AllSidechainSidechainCst)

    return AllConstraints, SortedConstraints


# sys.argv = [sys.argv[0]]+['-pdbs', '1rwr_Relax.pdb']

def main(argv=None):
  if argv != None:                                                             
    sys.argv =[ sys.argv[0] ]+[ arg for arg in argv ]
    
  ArgParser = argparse.ArgumentParser(description=' generate_cst.py arguments ( -help ) %s'%InfoString)
  # Required arguments:
  ArgParser.add_argument('-pdbs', type=str, nargs='+', help=' input pdbs ', required=True)
  # Optional arguments:
  ArgParser.add_argument('-out', type=str, help=' output directory ', default='./')
  ArgParser.add_argument('-max_dist', type=float, default=3.4, help=' distance between the oxygens and downstreams ')
  ArgParser.add_argument('-min_seq_sep', type=int, default=3, help=' minimum seperation in primary sequece ')
  ArgParser.add_argument('-upstream_atom', type=str, default='[ON]\w?\d?', help=' grep for upstream atoms ')
  ArgParser.add_argument('-downstream_atom', type=str, default='[ON]\w?\d?', help=' grep for downstream atoms ')
  ArgParser.add_argument('-num_repeats', type=int, default=5, help=' number of repeats to extrapolate contacts for ')
  ArgParser.add_argument('-weight',  type=float, default=1.0,  help=' weighting for constraints ')
  ArgParser.add_argument('-renumber_pose', type=bool, default=True, help='True|False renumber pdb residues ' )
  
  ArgParser.add_argument('-disulfide', type=bool, default=True, help='True|False also include disulfide constraints ' )  

  Args = ArgParser.parse_args()
  
  # if len(Args.pdbs[0]) == 1:
  #   Args.pdbs = [''.join(Args.pdbs)]

  if Args.out [-1] != '/':
    Args.out = Args.out + '/'

  import rosetta
  rosetta.init(extra_options = "-mute basic -mute core -mute protocols")

  ReportedRepeatCount = 0
  TotalPdbs = len(Args.pdbs)

  for iPdb, Pdb in enumerate(Args.pdbs):
    print ' Working with %s; %d of %d total pdbs '%(Pdb, iPdb+1, TotalPdbs)
    # Starting rosetta  
    Pose = rosetta.pose_from_pdb(Pdb)
    OutputPdb = Args.out+Pdb

    # Sets pdb info so residues in dumped pdbs are same as index 
    Pose.pdb_info(rosetta.core.pose.PDBInfo( Pose ))
    if Args.renumber_pose:
      rosetta.dump_pdb(Pose, OutputPdb)
    else:
      rosetta.dump_pdb(Pose, OutputPdb.replace('.pdb', '_renumbered.pdb'))

    AllConstraints, SortedConstraints = get_pose_constraints(Pose, Args.max_dist, Args.min_seq_sep, Args.upstream_atom, Args.downstream_atom, Args.weight, True)
    
    if Args.disulfide:
      DisulfAllConstraints, DisulfSortedConstraints = get_pose_constraints(Pose, 3.5, 2, 'SG', 'SG', Args.weight, False)
      AllConstraints.extend(DisulfAllConstraints)

    # print AllConstraints
    # print SortedConstraints
    # print 
    # print
    # print DisulfAllConstraints
    # print DisulfSortedConstraints
    # sys.exit()

    CstName = OutputPdb.replace('.pdb', '_All.cst')
    with open(CstName, 'w') as CstFile:
      print>>CstFile, '\n'.join(AllConstraints) 

    BackboneBackboneCst, BackboneSidechainCst, SidechainSidechainCst = SortedConstraints

    CstName = OutputPdb.replace('.pdb', '_BBBB.cst')
    with open(CstName, 'w') as CstFile:
      print>>CstFile, '\n'.join(BackboneBackboneCst) 
    
    CstName = OutputPdb.replace('.pdb', '_BBSC.cst')
    with open(CstName, 'w') as CstFile:
      print>>CstFile, '\n'.join(BackboneSidechainCst) 
    
    CstName = OutputPdb.replace('.pdb', '_SCSC.cst')
    with open(CstName, 'w') as CstFile:
      print>>CstFile, '\n'.join(SidechainSidechainCst) 
    
    CstName = OutputPdb.replace('.pdb', '_Disulf.cst')
    with open(CstName, 'w') as CstFile:
      print>>CstFile, '\n'.join(DisulfAllConstraints)

if __name__ == "__main__":
  sys.exit(main())


