#command example: pymol -c I32-10_arm_loop.pdb I32-10_stem_loop.pdb -qr ~/scripts/pymol/random_dock_optimization.py --chain_a I32-10_stem_loop --chain_b I32-10_arm_loop --geometry x20_ih

import sys
sys.path.append("/home/akhmelin/scripts/pymol")
import argparse
from pymol import cmd
import numpy as np
from scipy.spatial.distance import cdist
import os
#from os import path
from pymol_util import com
import math
from flatten_obj import flatten_obj

#chain_a corresponds to the trimer
#chain_b corresponds to the dimer
#geometry as defined by Reigun, e.g. x20_ih

def euler_angles_degrees(zvector, xvector):
#Return the three euler angles (in radians) that describe this HomogeneousTransform as the series
#of a Z axis rotation by the angle phi (returned in position 1 of the output vector), followed by
#an X axis rotation by the angle theta (returned in position 3 of the output vector), followed by another
#Z axis rotation by the angle psi (returned in position 2 of the output vector).
    zvector = zvector/np.linalg.norm(zvector)
    xvector = xvector/np.linalg.norm(xvector)
    yvector = np.cross(zvector, xvector)
    if zvector[2] >= 1:
        ang1 = - math.atan2(yvector[0],xvector[0])
        ang2 = 0
        ang3 = 0
    elif zvector[2] <= -1:
        ang1 = - math.atan2(yvector[0],xvector[0])
        ang2 = 0
        ang3 = math.pi	
    else:
        pos_sin_theta = math.sqrt(1 - zvector[2]*zvector[2])
        ang3 = math.asin(pos_sin_theta)
        if zvector[2] < 0:
            ang3 = math.pi - ang3
        ang1 = math.atan2(zvector[0], -zvector[1])
        ang2 = math.atan2(xvector[2], yvector[2])
    ang1 = ang1 * 180/math.pi
    ang2 = ang2 * 180/math.pi
    ang3 = ang3 * 180/math.pi
    return ang1, ang2, ang3

def over_dist(chain_a_i, chain_b, ab_con):
#calculate overlap distance between linkers in the bundle and arm
    aM = np.empty((0,3), dtype=float)
    aD = np.empty((0,3), dtype=float)
    for resi in [23,47,71]:
        cmd.do('select M{resi}a, {chain_a_i} and resi {resi} and name CA'.format(chain_a_i=chain_a_i, resi=resi))
        com_aM = com("M"+str(resi)+"a")
        aM = np.append(aM, [[com_aM.x, com_aM.y, com_aM.z]], axis = 0)
    for resi in [24,48,72]:
        cmd.do('select D{resi}a, {chain_a_i} and resi {resi} and name CA'.format(chain_a_i=chain_a_i, resi=resi))
        com_aD = com("D"+str(resi)+"a")
        aD = np.append(aD, [[com_aD.x, com_aD.y, com_aD.z]], axis = 0)
    bMs = np.empty((0,3), dtype=float)
    bDs = np.empty((0,3), dtype=float)
    for n in ab_con:
        cmd.do('select Mb, {chain_b}_{n} and resi 1+175 and name CA near_to 10 of {chain_a_i}'.format(chain_a_i=chain_a_i, chain_b=chain_b, n=n))
        cmd.do('select Db, {chain_b}_{n} and resi 2+176 and name CA near_to 10 of {chain_a_i}'.format(chain_a_i=chain_a_i, chain_b=chain_b, n=n))
        com_bM = com("Mb")
        com_bD = com("Db")
        bMs = np.append(bMs, [[com_bM.x,com_bM.y,com_bM.z]],axis=0)
        bDs = np.append(bDs, [[com_bD.x,com_bD.y,com_bD.z]],axis=0)
    Mdist = np.average(np.amin(cdist(aM,bMs)))
    Ddist = np.average(np.amin(cdist(aD,bDs)))
    return bMs, bDs, Mdist, Ddist

def rot_angle(chain_a, i, vector, bMs, bDs, step):
#calculate rotation angle for a particular bundle
    rot_aM = np.empty((0,3), dtype=float)
    rot_aD = np.empty((0,3), dtype=float)
    rot_Mdist = {}
    rot_Ddist = {}
    for ang in np.arange(-60,60,step):
        cmd.rotate(vector, ang, chain_a+"_"+str(i))
        for resi in [23,47,71]:
            rot_aM = np.append(rot_aM, [[com("M"+str(resi)+"a").x, com("M"+str(resi)+"a").y, com("M"+str(resi)+"a").z]], axis = 0)
        for resi in [24,48,72]:
            rot_aD = np.append(rot_aD, [[com("D"+str(resi)+"a").x, com("D"+str(resi)+"a").y, com("D"+str(resi)+"a").z]], axis = 0)
        rot_Mdist[ang] = np.average(np.amin(cdist(rot_aM,bMs)))
        rot_Ddist[ang] = np.average(np.amin(cdist(rot_aD,bDs)))
        cmd.rotate(vector, -ang, chain_a+"_"+str(i))
    min_Mdist = min(rot_Mdist.values())
    min_Ddist = min(rot_Ddist.values())
    angleM = min([key for key in rot_Mdist if rot_Mdist[key] == min_Mdist])
    angleD = min([key for key in rot_Ddist if rot_Ddist[key] == min_Ddist])
    print("Rotation angles found for ", chain_a+"_"+str(i), ". angleM = ", angleM, "and angleD = ", angleD)
    if angleM == angleD:
        angle = angleM
    else:
        angle = np.average(np.array([angleM, angleD]))
    return angle

def main(chain_a, chain_b, geometry):
#define total number of chains
    os.mkdir(geometry)

    num_a = int(geometry[1:3])
    num_b = int(3/2 * num_a)

#import vectors defining each subunit
    azvectors = np.genfromtxt(f"/home/akhmelin/I32-10/mathematical_models/vector_files/{geometry}_ztrimer.csv", delimiter=" ")
    bzvectors = np.genfromtxt(f"/home/akhmelin/I32-10/mathematical_models/vector_files/{geometry}_zdimer.csv", delimiter=" ")
    axvectors = np.genfromtxt(f"/home/akhmelin/I32-10/mathematical_models/vector_files/{geometry}_xtrimer.csv", delimiter=" ")
    bxvectors = np.genfromtxt(f"/home/akhmelin/I32-10/mathematical_models/vector_files/{geometry}_xdimer.csv", delimiter=" ")
    ab_con = np.genfromtxt(f"/home/akhmelin/I32-10/mathematical_models/vector_files/{geometry}_modcon.csv", delimiter=" ", dtype=int)	
#print(ab_con)	

#generate copies of each chain
    for i in range(num_a):
        cmd.do('create {chain_a}_{i}, {chain_a}'.format(chain_a=chain_a, i=i))
    cmd.do('delete {chain_a}'.format(chain_a=chain_a))
	
    for j in range(num_b):
        cmd.do('create {chain_b}_{j}, {chain_b}'.format(chain_b=chain_b, j=j))
    cmd.do('delete {chain_b}'.format(chain_b=chain_b))

#define selection for B carbons without the ovelap residues from each original chain
    cmd.do('select c_trimer, not resi 21+22+23+24+45+46+47+48+69+70+71+72 and {chain_a}_0 and name CB'.format(chain_a=chain_a))
    cmd.do('select c_dimer, not resi 1+2+3+4+175+176+177+178 and {chain_b}_0 and name CB'.format(chain_b=chain_b))
    cmd.do('select chains, all'.format())

#translate by r in [10:30]nm in 0.2 angstrong steps		
    rotation = 1
    overlap_dist={}
    for r in np.arange(50, 301, 0.1):
		
        if r == 50:
            if geometry == "x20_ih":
#in the first step, align to new axis and then translate
                print("Forming ", geometry, " symmetry")
                for i in range(num_a):
                    a_ang1, a_ang2, a_ang3 = euler_angles_degrees(azvectors[i], axvectors[i])
                    cmd.rotate([0,0,1], a_ang2, chain_a+"_"+str(i))
                    cmd.rotate([1,0,0], a_ang3, chain_a+"_"+str(i))
                    cmd.rotate([0,0,1], a_ang1, chain_a+"_"+str(i))
                    cmd.translate(np.array(r*azvectors[i]/np.linalg.norm(azvectors[i])).tolist(), chain_a+"_"+str(i))
                for j in range(num_b):
                    b_ang1, b_ang2, b_ang3 = euler_angles_degrees(bzvectors[j], bxvectors[j])
                    cmd.rotate([0,0,1], b_ang2, chain_b+"_"+str(j))
                    cmd.rotate([1,0,0], b_ang3, chain_b+"_"+str(j))
                    cmd.rotate([0,0,1], b_ang1, chain_b+"_"+str(j))
                    r2 = r + 0.23456134 #correction of dimer vs trimer translation
                    cmd.translate(np.array(r2*bzvectors[j]/np.linalg.norm(bzvectors[j])).tolist(), chain_b+"_"+str(j))
	    
            else:
#just as above but translate according to real zz vector with no correction of translational DOF
                print("Forming ", geometry, " symmetry")
                for i in range(num_a):
                    a_ang1, a_ang2, a_ang3 = euler_angles_degrees(azvectors[i], axvectors[i])
                    cmd.rotate([0,0,1], a_ang2, chain_a+"_"+str(i))
                    cmd.rotate([1,0,0], a_ang3, chain_a+"_"+str(i))
                    cmd.rotate([0,0,1], a_ang1, chain_a+"_"+str(i))
                    cmd.translate(np.array(r*azvectors[i]).tolist(), chain_a+"_"+str(i))
                for j in range(num_b):
                    b_ang1, b_ang2, b_ang3 = euler_angles_degrees(bzvectors[j], bxvectors[j])
                    cmd.rotate([0,0,1], b_ang2, chain_b+"_"+str(j))
                    cmd.rotate([1,0,0], b_ang3, chain_b+"_"+str(j))
                    cmd.rotate([0,0,1], b_ang1, chain_b+"_"+str(j))
                    cmd.translate(np.array(r*bzvectors[j]).tolist(), chain_b+"_"+str(j))

#in following steps, just translate by small increments
        else:
            for i in range(num_a):
                cmd.translate(np.array(0.1*azvectors[i]/np.linalg.norm(azvectors[i])).tolist(), chain_a+"_"+str(i))
            for j in range(num_b):
                cmd.translate(np.array(0.1*bzvectors[j]/np.linalg.norm(bzvectors[j])).tolist(), chain_b+"_"+str(j))

#get xyz coordinates for all carbon atoms in each subunit and calculate pairwise distances
#if minimum distance is below 3 ang, continue the loop
#if minimum distance equla or above 3 ang, rotate bundle to minimize distance between two
#overlaping residues of loop in bundle and arm
#after initial rotation, continue translation to minimize distance between overlapping residues
#when minimized, find translational distance for which overlapping distance is minimum,
#minimize all bundle rotations and then stop!
        xyz_ctrimer = cmd.get_coords("c_trimer", 1)
        xyz_cdimer = cmd.get_coords("c_dimer", 1)
        distmap = cdist(xyz_ctrimer, xyz_cdimer)
        min_distance = np.amin(distmap)
		
        if min_distance < 2.0:
            print("Clash found (min c-c distance = ",min_distance,") found for r=",r)
            continue
		
        elif 2.0 <= min_distance <= 4.0:
            print("No clashes (min c-c distance = ",min_distance,") found for r=",r)
            if rotation > 0:
                print("Look for appropriate bundle rotation")
                bMs, bDs, Mdist, Ddist = over_dist(chain_a+"_0", chain_b, ab_con[0])
                print("Distance between overlapping residues is: ", Mdist, " and ", Ddist)
                if np.average(np.array([Mdist, Ddist])) > 2:
                    step = 1
                if 1.5 < np.average(np.array([Mdist, Ddist])) <= 2:
                    step = 0.5
                if np.average(np.array([Mdist, Ddist])) <= 1.5:
                    step =0.2
                vector = azvectors[0].tolist()
                angle = rot_angle(chain_a, 0, vector, bMs, bDs, step)
                print("Starting rotation by angle = ", angle)
                cmd.rotate(vector, angle, chain_a+"_0")
                bMs, bDs, Mdist, Ddist = over_dist(chain_a+"_0", chain_b, ab_con[0])
                print("Distance between overlapping residues after rotation ",rotation," is: ", Mdist, " and ", Ddist)
                overlap_dist[r] = np.average(np.array([Mdist, Ddist]))
                print(overlap_dist)
                flatten_obj(name=geometry+'_'+str(r), selection="chains", state=0, rename=0, chain_map=geometry+'_map_'+str(r))
                cmd.save(geometry+"/"+geometry+"_"+str(r)+".pse", selection=geometry+"_"+str(r))
                print("Rotation cycle ", rotation, "finished")
                rotation = rotation + 1
                sampl_trans=list(overlap_dist.keys())
                sampl_over=list(overlap_dist.values())
                sampl_over=['%.2f' % elem for elem in sampl_over]
                if len(sampl_over)<6:
                    continue
                else:
                    if sampl_over[-1] < sampl_over[-6]:
                        continue
                    else:
                        print("Minimum overlap found!")
                        rf = sampl_trans[-1]
                        rmin = sampl_trans[sampl_over.index(min(sampl_over))]
                        for i in range(num_a):
                            cmd.translate(np.array(-rf*azvectors[i]).tolist(), chain_a+"_"+str(i))
                            cmd.translate(np.array(rmin*azvectors[i]).tolist(), chain_a+"_"+str(i))
                            print("Look for appropriate bundle ", i, " rotation")
                            bMs, bDs, Mdist, Ddist = over_dist(chain_a+"_"+str(i), chain_b, ab_con[i])
                            print("Distance between overlapping residues is: ", Mdist, " and ", Ddist)
                            vector = azvectors[i].tolist()
                            angle = rot_angle(chain_a, i, vector, bMs, bDs, 0.2)
                            print("Starting rotation by angle = ", angle)
                            cmd.rotate(vector, angle, chain_a+"_"+str(i))
                            bMs, bDs, Mdist, Ddist = over_dist(chain_a+"_"+str(i), chain_b, ab_con[i])
                            print("Distance between overlapping residues after rotation ",rotation," is: ", Mdist, " and ", Ddist)
                        for j in range(num_b):
                            cmd.translate(np.array(-rf*bzvectors[j]).tolist(), chain_b+"_"+str(j))
                            cmd.translate(np.array(rmin*bzvectors[j]).tolist(), chain_b+"_"+str(j))
                        flatten_obj(name=geometry+'_minimized_'+str(rmin), selection="chains", state=0, rename=0, chain_map=geometry+'_map_minimzed'+str(rmin))
                        cmd.save(geometry+"/"+geometry+"_minimized_"+str(rmin)+".pse", selection=geometry+"_minimized_"+str(rmin))
                        print("Job finished!")
                        break

            else:
                print("No rotation executed")
                bMs, bDs, Mdist, Ddist = over_dist(chain_a+"_0", chain_b, ab_con[0])
                print("Distance between overlapping residues is: ", Mdist, " and ", Ddist)
                overlap_dist[r] = np.average(np.array([Mdist, Ddist]))
                print(overlap_dist)
                flatten_obj(name=geometry+'_'+str(r), selection="chains", state=0, rename=0, chain_map=geometry+'_map_'+str(r))
                cmd.save(geometry+"/"+geometry+"_"+str(r)+".pse", selection=geometry+"_"+str(r))
                continue

        elif min_distance > 4.0:			
            print("No clashes (min c-c distance = ",min_distance,") found for r=",r)
            flatten_obj(name=geometry+'_'+str(r), selection="chains", state=0, rename=0, chain_map=geometry+'_map_'+str(r))
            cmd.save(geometry+"/"+geometry+"_"+str(r)+".pse", selection=geometry+"_"+str(r))
            print("Translation terminated")
            print("Looking for optimal configuration")
            sampl_trans=list(overlap_dist.keys())
            sampl_over=list(overlap_dist.values())
            sampl_over=['%.2f' % elem for elem in sampl_over]
            rf = sampl_trans[-1]
            rmin = sampl_trans[sampl_over.index(min(sampl_over))]
            for i in range(num_a):
                cmd.translate(np.array(-rf*azvectors[i]).tolist(), chain_a+"_"+str(i))
                cmd.translate(np.array(rmin*azvectors[i]).tolist(), chain_a+"_"+str(i))
                print("Look for appropriate bundle ", i, " rotation")
                bMs, bDs, Mdist, Ddist = over_dist(chain_a+"_"+str(i), chain_b, ab_con[i])
                print("Distance between overlapping residues is: ", Mdist, " and ", Ddist)
                vector = azvectors[i].tolist()
                angle = rot_angle(chain_a, i, vector, bMs, bDs, 0.2)
                print("Starting rotation by angle = ", angle)
                cmd.rotate(vector, angle, chain_a+"_"+str(i))
                bMs, bDs, Mdist, Ddist = over_dist(chain_a+"_"+str(i), chain_b, ab_con[i])
                print("Distance between overlapping residues after rotation ",rotation," is: ", Mdist, " and ", Ddist)
            for j in range(num_b):
                cmd.translate(np.array(-rf*bzvectors[j]).tolist(), chain_b+"_"+str(j))
                cmd.translate(np.array(rmin*bzvectors[j]).tolist(), chain_b+"_"+str(j))
            flatten_obj(name=geometry+'_minimized_'+str(rmin), selection="chains", state=0, rename=0, chain_map=geometry+'_map_minimzed'+str(rmin))
            cmd.save(geometry+"/"+geometry+"_minimized_"+str(rmin)+".pse", selection=geometry+"_minimized_"+str(rmin))
            print("Job finished!")
            break

def parse_args():

    parser = argparse.ArgumentParser(description="Find midpoint and respective xyz axis from connectivity and xyz file")
    parser.add_argument('--chain_a', type=str, required=True, help="trimer building block")
    parser.add_argument('--chain_b', type=str, required=True, help="dimer building block")
    parser.add_argument('--geometry', type=str, required=True, help="geometry to reach")
    args = parser.parse_args()
    chain_a = args.chain_a
    chain_b = args.chain_b
    geometry = args.geometry
    return chain_a, chain_b, geometry

if __name__=="__main__":
    main(*parse_args())

#cmd.extend('dock', dock)

