#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 14:13:18 2020

Parametrically generate xyz for arbitrarily tiled polyhedra, eventually make symdef
Currently takes the penton edge length, the in-plane edge length of the hexon (not
the 556-566 edge length) and h number where T=h**2 and writes out an xyz file with
C atoms placed on each penton vertex (556) and a P atom at each icosahedral 3-fold axis

@author: Quinton
"""
#TODO: GET THIS WORKING SO THAT IT CAN ACTUALLY RUN MAKEICOS
from xyzMath import Vec,Xform,Mat,rotation_around, alignaroundaxis
from pymol_util import trans, rot, xform
from symgen import makeicos, makesym
import numpy as np
from pymol import cmd, CmdException
from math import atan2

 
def is_even(number):
    if (number % 2) == 0:
        return True
    else:
        return False

def range_skip_nth(start, stop, nth, first):
    #create a list from start to stop, skipping every nth number starting with the first
    #for example, range_skip_nth(-5,2,3,0) will generate the list [-4,-3,-1,0,2,] skipping
    #-5, -2, and 1.
    out_list = []
    for i,n in enumerate(range(start,stop+1)):
        if not ((i - first) % nth) == 0:
            out_list.append(n)
    return out_list

def skew(x):
    return np.array([[0, -x[2], x[1]],
                     [x[2], 0, -x[0]],
                     [-x[1], x[0], 0]])

def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return np.sqrt(dotproduct(v, v))

def vec_angle(v1, v2):
  return np.arccos(dotproduct(v1, v2) / (length(v1) * length(v2)))

def isParallel(v1,v2):
    try:
        theta = np.round(vec_angle(v1, v2), 4) #theta needs to be reasonably accurate. Too accurate and this returns False too readily.
        if (theta == 0 or theta == 180):
            return True
    except ZeroDivisionError:
        return True

    return False

def describe_particle(h):
    #given an h number, assuming k=0
    #return the number of vertices 556, 556_2fold, 566, 666
    if h < 1:
        return
    c_num = 0
    b_num = 0

    for h in np.arange(1,h + 1):
        c_num = ( h - 1 ) + c_num
        if h > 2:
            b_num = ( h - 3 ) + b_num
    if h == 1:
        aab = 20
    else:
        aab = 60
    abb = ( h - 2 ) * 30 * 2
    ccc = c_num * 20
    bbb = b_num * 20
    
        
    if abb < 0:
        abb = 0 

    return aab, abb, ccc, bbb


def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

def rotate_about_vec(angle, vec, rad=True):
    w_matrix = skew((vec / np.linalg.norm(vec)).reshape(3))
    if not rad:
        angle = angle * ( np.pi / 180 )    
    rotation_matrix = np.eye(3) + np.sin(angle)*w_matrix + (2*np.sin( angle / 2)**2)*w_matrix.dot(w_matrix)
    return rotation_matrix

def normal_vec(plane):
    
    vecs_in_plane = [ plane[0] - plane[n] for n in range(1,len(plane))]
    cross_prods = [ np.cross(vecs_in_plane[0], vecs_in_plane[n]) for n in range(1,len(vecs_in_plane))]

    if len(plane) < 3:
        print("Two points does not a plane make")
        return None
    elif len(plane) > 3:
        all_in_plane = True#check that all of the points are in the plane
        for norm_vec in cross_prods: #if num chains is 3 then this will error.
            if not isParallel(cross_prods[0], norm_vec):
                all_in_plane = False
        if not all_in_plane:
            print ("All points are not in the same plane")
            return None
    else:
        return ( cross_prods[0] / np.linalg.norm(cross_prods[0])).reshape(3)

def translate(to_translate, vec, mag):
    #takes an object to translate, a vector for direction, and a magnitude
    #returns the object appropriately translated along vec
    
    #first scale the vector appropriately
    to_translate
    dX, dY, dZ = (mag / np.linalg.norm(vec)) * np.array(vec)
    translate_matrix = np.array([[1, 0, 0, dX],
                               [0, 1, 0, dY],
                               [0, 0, 1, dZ],
                               [0, 0, 0, 1]])
    
    for i,point in enumerate(to_translate):
        point = np.append(point, 1)
        to_translate[i] = np.delete(translate_matrix.dot(point.T).T, -1)
    return to_translate

def renumber(selection='all', start=1, startsele=None, quiet=1):
    '''
DESCRIPTION

    Set residue numbering (resi) based on connectivity.

ARGUMENTS

    selection = string: atom selection to renumber {default: all}

    start = integer: counting start {default: 1}

    startsele = string: residue to start counting from {default: first in
    selection}
    '''
    start, quiet = int(start), int(quiet)

    model = cmd.get_model(selection)
    cmd.iterate(selection, 'next(atom_it).model = model',
            space={'atom_it': iter(model.atom), 'next': next})
    if startsele is not None:
        startidx = cmd.index('first (' + startsele + ')')[0]
        for atom in model.atom:
            if (atom.model, atom.index) == startidx:
                startatom = atom
                break
        else:
            print(' Error: startsele not in selection')
            raise CmdException
    else:
        startatom = model.atom[0]
    for atom in model.atom:
        atom.adjacent = []
        atom.visited = False
    for bond in model.bond:
        atoms = [model.atom[i] for i in bond.index]
        atoms[0].adjacent.append(atoms[1])
        atoms[1].adjacent.append(atoms[0])
    minmax = [start, start]

    def traverse(atom, resi):
        atom.resi = resi
        atom.visited = True
        for other in atom.adjacent:
            if other.visited:
                continue
            if (atom.name, other.name) in [('C', 'N'), ("O3'", 'P')]:
                minmax[1] = resi + 1
                traverse(other, resi + 1)
            elif (atom.name, other.name) in [('N', 'C'), ('P', "O3'")]:
                minmax[0] = resi - 1
                traverse(other, resi - 1)
            elif (atom.name, other.name) not in [('SG', 'SG')]:
                traverse(other, resi)
    traverse(startatom, start)
    
    cmd.alter(selection, 'resi = next(atom_it).resi',
            space={'atom_it': iter(model.atom), 'next': next})
    if not quiet:
        print(' Renumber: range (%d to %d)' % tuple(minmax))

def get_sequence(sele):
    aa1 = list("ACDEFGHIKLMNPQRSTVWY")
    aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
    aa321 = dict(zip(aa3,aa1))
    sequence = [ aa321[i.resn] for i in cmd.get_model(sele + " and n. ca").atom ]
    return sequence

#TODO: Clean this up, it's overly complicated
def renumber_across_chains(sele, make_single_chain=True):
    chains = {}
    for i,c in enumerate(cmd.get_chains(sele)):
        seq = get_sequence(f"{sele} and chain {c}")
        seq_len = len(seq)
        chains[c] = {"seq" : seq, "length" : seq_len}
    
    start = 0
    for c in chains.keys():
        chains[c]["start"] = start
        chains[c]["end"] = start + chains[c]["length"]
        start = start + chains[c]["length"] + 1
    
    for c in chains.keys():
        print(f"RENUMBERING: Chain {c} in {sele} start resi {chains[c]['start']}")
        renumber(f"{sele} and chain {c}", quiet=False)
        cmd.alter(f"{sele} and chain {c}", f"resi=str(int(resi) + {chains[c]['start']})")

    if make_single_chain==True:
        cmd.alter(sele, "chain='A'")
    

"""
Building the polyhedron from an underlying icosahedron
""" 
    
phi = (1 + np.sqrt(5)) / 2
dodec_centers = np.array([[0,1,phi],
                 [0,-1,phi],
                 [0,-1,-phi],
                 [0,1,-phi],
                 [phi,0,1],
                 [-phi,0,1],
                 [-phi,0,-1],
                 [phi,0,-1],
                 [1,phi,0],
                 [-1,phi,0],
                 [-1,-phi,0],
                 [1,-phi,0]])

icos_centers = np.array([[1,1,1],
                [-1,1,1],
                [-1,-1,1],
                [-1,-1,-1],
                [1,-1,-1],
                [1,-1,1],
                [1,1,-1],
                [-1,1,-1],
                [0,phi,1/phi],
                [0,-phi,1/phi],
                [0,-phi,-1/phi],
                [0,phi,-1/phi],
                [1/phi,0,phi],
                [-1/phi,0,phi],
                [-1/phi,0,-phi],
                [1/phi,0,-phi],
                [phi,1/phi,0],
                [-phi,1/phi,0],
                [-phi,-1/phi,0],
                [phi,-1/phi,0]])


#p_comp1 = cmd.centerofmass("chain A")
#p_comp2 = cmd.centerofmass("chain C")
#a_p = np.linalg.norm(p_comp1 - p_comp2) 
#are not flat...



def makegoldberg(penton, trimer, h_number, het_resis, name="GOLD", n=1, align=False):
    #penton is a selection of all five pseudosymmetric trimers in the penton, where each trimer is labeled as one chain
    #trimer is a selection of all three chains in the homotrimer component
    #h_number determines the desired T_number by T=h_number**2
    #name is the name used for the full cages
    #n is the number of ASUs to display
    #het_resis is a list of touples with the start and stop resis for the "A" and "B" chains in the heterotrimer.
    
    tmp_penton = "TMP_penton"
    tmp_trimer = "TMP_trimer"
    cmd.delete(tmp_penton)
    cmd.delete(tmp_trimer)
    cmd.delete("tmp*")
    cmd.delete("GOLD*")
    cmd.create(tmp_penton, penton)
    cmd.create(tmp_trimer, trimer)
    PENTONCHAINS="ACEGI"
    TRIMERCHAINS="BDF"
    for i,c in enumerate(cmd.get_chains(tmp_penton)):
        print(f"{tmp_penton} and chain {c} chain='{PENTONCHAINS[i]}'")
        cmd.alter(f"{tmp_penton} and chain {c}", f"chain='{PENTONCHAINS[i]}'")
    seq_len = 0    
    for i,c in enumerate(cmd.get_chains(tmp_trimer)):
        seq = get_sequence(f"{tmp_trimer} and chain {c}")
        #rename chain
        print(f"Trimer sequence is {len(seq)} AA long, adding {seq_len}")
        cmd.alter(f"{tmp_trimer} and chain {c}", f"resi=str(int(resi) + {seq_len})")
        cmd.alter(f"{tmp_trimer} and chain {c}", f"chain='{TRIMERCHAINS[0]}'")
        
        seq_len = seq_len + len(seq)

    #Store an extra copy of the tmp_penton A chain for later use
    cmd.create("TMP_penton_A", f"{tmp_penton} and chain A")
    
    
    #get hexon and penton lengths from the design model
    pent_center = Vec(cmd.centerofmass(tmp_penton))
    #TODO: This should detect nearest trimer to penton, as of now relies on the selected trimer to be in contact with chain A.
    p1vec = Vec(cmd.centerofmass(f"{tmp_penton} and chain A"))
    p2vec = Vec(cmd.centerofmass(f"{tmp_penton} and chain C"))
    t1vec = Vec(cmd.centerofmass(tmp_trimer))
    oop_vec = p1vec - t1vec
    h_vec = (t1vec.dot(oop_vec) / t1vec.length()**2) * t1vec
    d_vec = -(oop_vec - h_vec)
    h = h_vec.length()
    d = d_vec.length()
    
    #get the rotation matrix to align the sheet to the y axis
    rot_mat = rotation_matrix_from_vectors(t1vec, np.array([0,0,1]))

    #Check that the trimer is "inside" or "outside" the cage    
    if h_vec.dot(t1vec) < 0:        
        h = -h
        
    a_hx = d
    a_p = (p1vec - p2vec).length()
    print(f"The penton is {pent_center.length()} A from the center.")
    print(f"The penton A chain is {p1vec.length()} A from the center and the penton distance is {a_p}")
    print(a_hx)
        
    """
    #Find the circumradius of the base icosahedron given the T-number
    #given the edge length of the penton a_p (distance between adjacent 566 vertices)
    #and the edge length of the hexon a_hx (distance between 566 vertex and adjacent 666 vertex)
    #and the h number from T = h**2 + hk + k**2 and k=0 always, then the circumradius Cr depends on a_p, a_hx, and h by
    
    edge is the edge length of the underlying icosahedron, ie. the edge length of the 3-fold facet, plus the length of a
    line drawn along the two-fold symmetry axis from the edge of the facet to the C5 symmetry axis
    the edge length of the facet depends on the length of side of the hexon by:
    
        eq. 1   sqrt(3) * a_hx * (h-1)
        
    and the length of the penton edge by
    
        eq. 2   2 * a_p * (3 / np.sqrt(6 + 2*np.sqrt(5))))
        
        Imagine the underlying icosahedron has edge lengths that are l longer than the facets defined by eq.1. 
        l depends on the penton edge length a_p by
        
        eq. 3   l = l' / sin(30)  = 3/sqrt(3) * l'
        
        and
        
        eq. 4   l'**2 = (1/2*a_p)**2 + h**2
        
        where
        
        eq. 5   h = (1/2)*a_p tan(theta / 2)
        
        and theta = the icosahedral dihedral angle. Thus
        
        eq. 6   h = (1/2)* a_p * sqrt((3 - sqrt(5)) / (3 + sqrt(5)))
        
        Substituting equation 4 into 3 and 6 into 4, and simplifying, we get eq. 2
    
    """
    edge = (np.sqrt(3) * a_hx * (h_number-1)) + (2 * a_p * (3 / np.sqrt(6 + 2*np.sqrt(5))))
     
    phi = ( 1 + np.sqrt(5) ) / 2
    
    """
    The circumradius of an icosahedron, R, is:
        
        eq.1a    R = (edge / 2) * np.sqrt(phi * np.sqrt(5)) 
        
        the inradius is:
        
        eq.1b    Ci = (phi**2 * edge) / (2*np.sqrt(3) )
        
        where edge is the edge length of the icosahedron, discussed above, and phi is the golden ration = ( 1 + np.sqrt(5) ) / 2
       
        For a goldberg particle this radius needs to be corrected dependent on the size of the penton.
        
        eq.2    Cr = R - Cr_p
        
        From the trigonometric identeties of half-angles we derive the following:
            
            eq.3a    sin(theta) = Ci/R
            eq.3b    cos(theta) = np.sqrt(R**2 - Ci**2) / R
            
        And therefore
        
            eq.4    Cr_p = rp * cos(theta) / sin(theta)
            
            where:
                
            eq.5    rp = np.sqrt(2 / (5 - np.sqrt(5))) * a_p
            
            Substituting eq. 5, and 3a and 3b into 4 and 4 into 2 yields 
        
            eq.6    Cr = R - (( np.sqrt(R**2 - Ci**2) / Ci ) * r_p)
            
    """
    Ci = (phi**2 * edge) / (2*np.sqrt(3) )
    R = ( edge / 2 ) * np.sqrt(phi * np.sqrt(5))
    r_p = np.sqrt(2 / (5 - np.sqrt(5))) * a_p
    Cr = R - (( np.sqrt(R**2 - Ci**2) / Ci ) * r_p)
    
    #Translate the penton along it's axis some distance 
    #Cr is the final radius, so translate Cr less it's current radius.
    print(f"Cage radius {Cr}")
    print(f"{type(Cr)}")
    tr_pent_dist = Cr - pent_center.length()
    pent_unit = pent_center.normalized()
    
    print(tr_pent_dist)
    trans(tmp_penton, pent_unit * float(tr_pent_dist))
    
    #   Use the hexagonal grid notation to fill out one asu. Given all possible values of 
    #   i*q, j*r, k*s, the ASU lattice contains all trimers where i > 0 and j,k <= 0.
    #   meeting certain requirements for the values of i, j, and k
    #1. Create a list of unit vectors that point to every three-fold position in the facet
    #   That is, 
    
    #Define basis vectors for 2D array.
    #Transform d_vec into plane normal to t1_vec and including origin
    n1 = t1vec.normalized()
    q = d_vec - ( ( d_vec.dot(n1) ) * n1)
    #q = d_vec
    r = rotation_around(t1vec, 2*np.pi/3) * q
    s = rotation_around(t1vec, -2*np.pi/3) * q

    hexons = []
    hex_coords = []
    hexon_type = []
    a = -( h_number - 2 )
    b = int((h_number-1) / 2)
    for i in range( a, b + 1 ):
        for j in range( a, b + 1 ):
            for k in range( a, b + 1 ):
                sum_ijk = sum([i,j,k])
                abs_sum = sum([abs(i), abs(j), abs(k)])
                if ((i <= 0 and (j >= 0 and k >=0)) or (j <= 0 and (i >= 0 and k >=0)) or (k <= 0 and (i >= 0 and j >=0))) and (i*j*k == 0):
                    if (abs_sum <= abs(a)) and (sum_ijk >= a) and (sum_ijk <= b) and (sum_ijk != a + 1):
                            #if not (abs_sum == 1 and sum_ijk == -1):
                            #if (i < 0 and (abs(i) >= j and abs(i) >= k) and (abs(i + 1) != j) and (abs(i + 1) != k)) or (j < 0 and (abs(j) >= i and abs(j) >= k)  and (abs(j + 1) != i) and (abs(j + 1) != k)) or (k < 0 and (abs(k) >= i and abs(k) >= j) and (abs(k + 1) != i) and (abs(k + 1) != j)) :
                            #also, allowed points, like [-7,0,2] are off the edge of the penton, so this needs to be handled differently
                            #The edge of the facet is for example [a:0,0,0:b]
                            
                            #if i,j,k in this list then min(j,k) == 0, else == 0
                            zero_start_list = range_skip_nth(a,0,3,1)
                            
                            #the maximum allowed j,k given i <= 0 is defined by int( ( a - (i - 1)) / 2 )
                            max_k_i = int( ( (i + 1) - a) / 2 )
                            max_j_i = int( ( (i + 1) - a) / 2 )
                            
                            max_i_j = int( ( (j + 1) - a) / 2 )
                            max_k_j = int( ( (j + 1) - a) / 2 )
                            
                            max_i_k = int( ( (k + 1) - a) / 2 )
                            max_j_k = int( ( (k + 1) - a) / 2 )
                            
                            
                            if i in zero_start_list:
                                if i in zero_start_list[::2]:
                                    allowed_k_i = range_skip_nth(0,max_k_i,3,1)
                                    allowed_j_i = range_skip_nth(0,max_k_i,3,1)
                                else:
                                    allowed_k_i = range_skip_nth(0,max_k_i,3,2)
                                    allowed_j_i = range_skip_nth(0,max_k_i,3,2)
                            else:
                                allowed_k_i = range_skip_nth(1,max_k_i,3,2)
                                allowed_j_i = range_skip_nth(1,max_j_i,3,2)
                            if j in zero_start_list:
                                if j in zero_start_list[::2]:
                                    allowed_i_j = range_skip_nth(0,max_i_j,3,1)
                                    allowed_k_j = range_skip_nth(0,max_k_j,3,1)
                                else:
                                    allowed_i_j = range_skip_nth(0,max_i_j,3,2)
                                    allowed_k_j = range_skip_nth(0,max_k_j,3,2)
                            else:
                                allowed_k_j = range_skip_nth(1,max_k_j,3,2)
                                allowed_i_j = range_skip_nth(1,max_i_j,3,2)
                            if k in zero_start_list:
                                if k in zero_start_list[::2]:
                                    allowed_i_k = range_skip_nth(0,max_i_k,3,1)
                                    allowed_j_k = range_skip_nth(0,max_j_k,3,1)
                                else:
                                    allowed_i_k = range_skip_nth(0,max_i_k,3,2)
                                    allowed_j_k = range_skip_nth(0,max_j_k,3,2)
                            else:
                                allowed_i_k = range_skip_nth(1,max_i_k,3,2)
                                allowed_j_k = range_skip_nth(1,max_j_k,3,2)
                            if (i <= 0 and ( ((k in allowed_k_i) and j==0) or ((j in allowed_j_i) and k==0))) or (j <= 0 and ( ((i in allowed_i_j) and k==0) or ((k in allowed_k_j) and i==0))) or (k <= 0 and ( ((i in allowed_i_k) and j==0) or ((j in allowed_j_k) and i==0))):
                                
                                #TODO: This vector needs to be modifide to get the in-plane vector
                                #
                                hexons.append(i*q + j*r + k*s)
                                hex_coords.append([i,j,k])
                                neg_list = [z for z in range(a+1,1)][::2]
                                #For the edge of the facet, a < a, j, i < b + 1
                                if (i+j+k) in [z for z in range(a,b+1)][::3]: #All chain C homotrimer positions
                                    print("Adding trimer of type CCC to facet")
                                    hexon_type.append("CCC")
                                elif ((i in neg_list) and (j == neg_list.index(i) + 1) and k==0) or ((i in neg_list) and (k == neg_list.index(i) + 1) and j==0) or ((j in neg_list) and (i == neg_list.index(j) + 1) and k==0) or ((j in neg_list) and (k == neg_list.index(j) + 1) and i==0) or ((k in neg_list) and (i == neg_list.index(k) + 1) and j==0) or ((k in neg_list) and (j == neg_list.index(k) + 1) and i==0):  #All chain ABB heterotrimer positions    
                                    print("Adding trimer of type ABB to Edge")    
                                    hexon_type.append("ABB")
                                else: #Otherwise the position is a BBB homotrimer
                                    print("Adding trimer of type BBB to Edge") 
                                    hexon_type.append("BBB")

    #Create a copy of the penton trimer and split it up into a true heterotrimer
    #Use this for making AAB, ABB, and BBB
    c = ["AAB_A", "AAB_B"]
    for i, resis in enumerate(het_resis):
        cmd.create(c[i], f"penton and chain A and resi {resis[0]}-{resis[1]}")
        trans(c[i], -p1vec)
    #MAKE BBB Component
    cmd.create("TMP_B_1", "AAB_B")
    cmd.create("TMP_B_2", "AAB_B")
    rot("TMP_B_1", t1vec, 120)#2*np.pi/3)
    rot("TMP_B_2", t1vec, -120)#-2*np.pi/3)
    cmd.alter("TMP_B_1", "chain='B'")
    cmd.alter("TMP_B_2", "chain='C'")
    cmd.create("TMP_BBB", "AAB_B or TMP_B_1 or TMP_B_2")
    cmd.delete("TMP_B_1")
    cmd.delete("TMP_B_2")
    
    #MAKE ABB component from one copy of AAB_A and two copies of AAB_B
    cmd.create("TMP_B_1", "AAB_B")
    #Measure the rotation of the trimer around the axis.

    avec = Vec(cmd.centerofmass("AAB_A"))
    bvec = Vec(cmd.centerofmass("AAB_B"))
    ab_angle = atan2((bvec.cross(avec)).dot(n1), bvec.dot(avec))        
        
    if ab_angle > 0:
        rot("AAB_B", t1vec, -120)
    else:
        rot("AAB_B", t1vec, 120)
    #chain A and B are comp "B" and chain C is comp "A"
    cmd.alter("TMP_B_1", "chain='B'")
    cmd.alter("AAB_A", "chain='C'")
    
    cmd.create("TMP_ABB", "AAB_B or AAB_A or TMP_B_1")
    cmd.delete("AAB_B")
    cmd.delete("AAB_A")
    cmd.delete("TMP_B_1")
        
    #else:
    #    #TODO: Use the het-trimer as a stand in for homotrimer and single component.
    #    #TODO: Translate het to origin in 3-fold plane
    
    #Move the homotrimer to the origin
    trans(tmp_trimer, -t1vec)
    
    #Iterate through hexon positions, make copies of appropriate component, and translate onto position.
    b_met = False
    alignment_obj = Xform(Mat(Vec(0.000000, 0.000000, 0.000000), Vec(0.000000, 0.000000, 0.000000), Vec(0.000000, -0.000000, 0.000000)), Vec(0.000000, 0.000000, 0.000000))
    basis_vector = Vec(x=0, y=0, z=0)
    for i, hex_vec in enumerate(hexons):
        
        ident = hexon_type[i]
        
        if ident == "CCC":
            if (hex_coords[i][0] <= 0) & (hex_coords[i][1] >= 0) & (hex_coords[i][2] >= 0) & (not b_met):
                if (sum(hex_coords[i]) == b):
                    b_met = True
                hex_trans_vec = hex_vec + float(Ci -  h)*n1
                #print(f"Hex plane distance is {Ci} and trimer distance is {Ci - h}")
                #print(f"Translating along vector by {n1.length()}")
                if (all(coord == 0  for coord in hex_coords[i])):
                    print("Trimer at center of facet")
                    cmd.create("TMP", f"{trimer}")
                    chains = cmd.get_chains("TMP")
                    for c,chain in enumerate(chains):
                        cmd.alter(f"TMP and chain {chain}", f"chain='{TRIMERCHAINS[c]}'")
                    cmd.create(f"tmp_trimer_{i}", f"TMP and chain {TRIMERCHAINS[0]}")
                    cmd.delete("TMP")
                    trans(f"tmp_trimer_{i}", -t1vec)
                else:
                    cmd.create(f"tmp_trimer_{i}",tmp_trimer)
                trans(f"tmp_trimer_{i}", hex_trans_vec)
                print(f"Creating tmp_trimer_{i} at coord {hex_coords[i]}")
                
                if i == 0:
                    #check that the hex_trans_vec in plane with the tri_vec, pent_vec plane
                    tmp_penton_vec = Vec(cmd.centerofmass(f"{tmp_penton} and chain A"))
                    pt_align_vec = tmp_penton_vec - t1vec
                    alignment_obj = alignaroundaxis(t1vec, hex_trans_vec, pt_align_vec)
                    basis_vector = hex_vec
                    print(f"alignment {alignment_obj}")
                if align:
                    xform(f"tmp_trimer_{i}", alignment_obj)
                    
                if n > 1:
                    renumber_across_chains(f"tmp_trimer_{i}")
                    makeicos(f"tmp_trimer_{i}", n=n, name=f"{name}_trimer_{i}")
                    cmd.delete(f"tmp_trimer_{i}")
                    
        elif ident == "ABB":
            if (hex_coords[i][0] <= 0) & (hex_coords[i][1] >= 0) & (hex_coords[i][2] >= 0) & (not b_met):
                if (sum(hex_coords[i]) == b):
                    b_met = True
                if (hex_coords[i][1] > 0):
                    if ab_angle < 0:
                        rotation=-1 #Clockwise rotation necessary to align with facet
                    else:
                        rotation=0
                elif (hex_coords[i][2] >= 0):
                    if ab_angle >  0:
                         #Counter-clockwise rotation necessary to align with facet
                        rotation=1
                    else:
                        rotation=0
                else:
                    print("WARNING: ROTATION COULD NOT BE DETERMINED")
                
                tmp_name = f"tmp_ABB_trimer_{i}"
                cmd.create(tmp_name,"TMP_ABB")
                rot(tmp_name, t1vec, rotation*120)
                #trans(tmp_name, -p1vec)
                hex_trans_vec = hex_vec + float(Ci)*n1
                trans(tmp_name, hex_trans_vec)
                print(f"Creating tmp_ABB_trimer_{i} at coord {hex_coords[i]} with rotation {rotation}")
                if align:
                    #get angle between basis_vector and hex_vector
                    #rot_angle = basis_vector.angle(hex_vec)
                    #align_vec = rotation_around(t1vec, rot_angle) * pt_align_vec
                    #alignment_obj = alignaroundaxis(t1vec, hex_trans_vec, align_vec)
                    xform(f"tmp_ABB_trimer_{i}", alignment_obj)
                if n > 1:
                    renumber_across_chains(f"tmp_ABB_trimer_{i}")
                    makeicos(f"tmp_ABB_trimer_{i}", n=n, name=f"{name}_AAB_trimer_{i}")
                    cmd.delete(f"tmp_ABB_trimer_{i}")
        else:
            if (hex_coords[i][0] <= 0) & (hex_coords[i][1] >= 0) & (hex_coords[i][2] >= 0):
                    
                if (all(coord == 0  for coord in hex_coords[i])):
                    print("BBB trimer at center of facet")
                    tmp_name = f"tmp_BBB_trimer_{i}"
                    cmd.create(tmp_name,"TMP_BBB and chain A")
                    hex_trans_vec = hex_vec + float(Ci)*n1
                    trans(tmp_name, hex_trans_vec)
                    print(f"Creating tmp_BBB_trimer_{i} at coord {hex_coords[i]}")
                    if align:
                        #if hex_coords[i][1] + hex_coords[i][2] == 0:
                        #    rot_angle = 0
                        #else:
                        #    rot_angle = basis_vector.angle(hex_vec)
                        #align_vec = rotation_around(t1vec, rot_angle) * pt_align_vec
                        #alignment_obj = alignaroundaxis(t1vec, hex_trans_vec, align_vec)
                        xform(f"tmp_BBB_trimer_{i}", alignment_obj)
                    if n > 1:
                        renumber_across_chains(f"tmp_BBB_trimer_{i}")
                        makeicos(f"tmp_BBB_trimer_{i}", n=n, name=f"{name}_BBB_trimer_{i}")
                        cmd.delete(f"tmp_BBB_trimer_{i}")                
                
                elif (hex_coords[i][0] < 0) & (hex_coords[i][1] >= 0) & (hex_coords[i][2] >= 0) & (not b_met):
                    tmp_name = f"tmp_BBB_trimer_{i}"
                    cmd.create(tmp_name,"TMP_BBB")
                    hex_trans_vec = hex_vec + float(Ci)*n1
                    trans(tmp_name, hex_trans_vec)
                    print(f"Creating tmp_BBB_trimer_{i} at coord {hex_coords[i]}")
                    if align:
                        #rot_angle = basis_vector.angle(hex_vec)
                        #align_vec = rotation_around(t1vec, rot_angle) * pt_align_vec
                        #alignment_obj = alignaroundaxis(t1vec, hex_trans_vec, align_vec)
                        xform(f"tmp_BBB_trimer_{i}", alignment_obj)
                    if n > 1:
                        renumber_across_chains(f"tmp_BBB_trimer_{i}")
                        makeicos(f"tmp_BBB_trimer_{i}", n=n, name=f"{name}_BBB_trimer_{i}")
                        cmd.delete(f"tmp_BBB_trimer_{i}")
                        
                elif (hex_coords[i][0] == 0) & (hex_coords[i][1] == 0) & (hex_coords[i][2] >= 0):
                    tmp_name = f"tmp_BBB_trimer_{i}"
                    cmd.create(tmp_name,"TMP_BBB")
                    hex_trans_vec = hex_vec + float(Ci)*n1
                    trans(tmp_name, hex_trans_vec)
                    print(f"Creating tmp_BBB_trimer_{i} at coord {hex_coords[i]}")
                    if align:
                        #rot_angle = basis_vector.angle(hex_vec)
                        #align_vec = rotation_around(t1vec, rot_angle) * pt_align_vec
                        #alignment_obj = alignaroundaxis(t1vec, hex_trans_vec, align_vec)
                        xform(f"tmp_BBB_trimer_{i}", alignment_obj)
                    if n > 1:
                        renumber_across_chains(f"tmp_BBB_trimer_{i}")
                        makeicos(f"tmp_BBB_trimer_{i}", n=n, name=f"{name}_BBB_trimer_{i}")
                        cmd.delete(f"tmp_BBB_trimer_{i}")
    #cmd.create(tmp_trimer, trimer)
    #translate Ci length.
    #cmd.delete("TMP*")
    if n > 1:
        renumber_across_chains(f"{tmp_penton} and chain A")
        makeicos(f"{tmp_penton} and chain A", n=n, name=f"{name}_penton")
        cmd.delete(tmp_penton)
    cmd.delete(tmp_trimer)
    cmd.delete("TMP_penton_A")
    cmd.delete("TMP_BBB")
    cmd.delete("TMP_ABB")
