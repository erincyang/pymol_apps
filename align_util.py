#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 18:27:08 2020

@author: Quinton
@edited by: Erin Yang on 2020.08.20
usage in pymol: align_sele_to_vec("cage",cmd.centerofmass("sele"),[0,0,1])
"""

from pymol import cmd
import numpy as np

def trans(sel, v):
    if xyz.isvec(v):
        cmd.translate([v.x, v.y, v.z], sel, 0, 0)
    elif xyz.isnum(v):
        cmd.translate([v, v, v], sel, 0, 0)
    else:
        raise NotImplementedError

def rot(sel, axis, ang, cen=xyz.V0):
    if not xyz.isvec(axis):
        raise NotImplementedError
    if not xyz.isnum(ang):
        raise NotImplementedError
    if not xyz.isvec(cen):
        raise NotImplementedError
    if cen is None:
        cen = com(sel)
    # if abs(axis.x) < 0.00001: axis.x = 0.0
    # if abs(axis.y) < 0.00001: axis.y = 0.0
    # if abs(axis.z) < 0.00001: axis.z = 0.0[[]]
    # cmd.rotate([round(axis.x,5), round(axis.y,5), round(axis.z,5)], ang, sel, 0, 0, None, [round(cen.x,5), round(cen.y,5), round(cen.z,5) ])
    axis = "[ %9.6f, %9.6f, %9.6f ]" % (axis.x, axis.y, axis.z)
    cen = "[ %9.6f, %9.6f, %9.6f ]" % (cen.x, cen.y, cen.z)
    cmd.rotate(axis, ang, sel, 0, 0, None, cen)


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

def align_sele_to_vec(sele, sele_vec, align_vec):
    vec1 = sele_vec
    vec2 = align_vec

    rot_mat = rotation_matrix_from_vectors(vec1, vec2)

    translate_mat1 = np.append(rot_mat[0], vec2[0])
    translate_mat2 = np.append(rot_mat[1], vec2[1])
    translate_mat3 = np.append(rot_mat[2], vec2[2])

    translate_mat_A = np.append(np.append(translate_mat1, np.append(translate_mat2, translate_mat3)), [0,0,0,1])

    #translate penton onto cage
    cmd.transform_object(sele, translate_mat_A)

def vec_angle(v1, v2):
    angle = np.arccos(dotproduct(v1, v2) / (length(v1) * length(v2))) 
    print( angle ) 
    return angle

def rotate_cage(cage, align_vec1, align_vec2, rot_axis):
	c2_1 = f'chain F and {cage} + chain T and {cage}'
	align_sele_to_vec(cage, cmd.centerofmass(c2_1), align_vec1)

	trans(cage, -Vec(cmd.centerofmass(cage)))

	c2_2 = f'chain B and {cage} + chain H and {cage}'
	ang_rad = vec_angle(cmd.centerofmass(c2_2), align_vec2)
	ang_deg = 180/np.pi*ang_rad

	rot(cage, rot_axis, -ang_deg)

def dist_to_plane(plane, point):
    norm_vec = normal_vec(plane)
    X = point - plane[1]
    D = norm_vec.dot(X)
    return D

def coplanar(plane):
    #check that all points are in a plane
    vecs_in_plane = [ plane[0] - plane[n] for n in range(1,len(plane))]
    cross_prods = [ np.cross(vecs_in_plane[0], vecs_in_plane[n]) for n in range(1,len(vecs_in_plane))]
    
    if len(plane) < 3:
        print("Two points does not a plane make")
        return None
    elif len(plane) == 3:
        print ("All points are in the same plane because only three points were provided")
        return True
    else:
        all_in_plane = True#check that all of the points are in the plane
        for norm_vec in cross_prods: #if num chains is 3 then this will error.
            if not isParallel(cross_prods[0], norm_vec):
                all_in_plane = False
        if not all_in_plane:
            print ("All points are not in the same plane")
            return False

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

def rotate_about_vec(angle, vec, rad=True):
    w_matrix = skew((vec / np.linalg.norm(vec)).reshape(3))
    if not rad:
        angle = angle * ( np.pi / 180 )    
    rotation_matrix = np.eye(3) + np.sin(angle)*w_matrix + (2*np.sin( angle / 2)**2)*w_matrix.dot(w_matrix)
    return rotation_matrix

def normal_vec(plane):
    vecs_in_plane = [plane[0] - plane[n] for n in range(1, len(plane))]
    cross_prods = [np.cross(vecs_in_plane[0], vecs_in_plane[n]) for n in range(1, len(vecs_in_plane))]

    if len(plane) < 3:
        print("Two points does not a plane make")
        return None
    elif len(plane) > 3:
        all_in_plane = True  # check that all of the points are in the plane
        for norm_vec in cross_prods:  # if num chains is 3 then this will error.
            if not isParallel(cross_prods[0], norm_vec):
                all_in_plane = False
        if not all_in_plane:
            print("All points are not in the same plane")
            return None
    else:
        return (cross_prods[0] / np.linalg.norm(cross_prods[0])).reshape(3)

'''
def rotate_about_vec(angle, vec, rad=True):
    w_matrix = skew((vec / np.linalg.norm(vec)).reshape(3))
    if not rad:
        angle = angle * ( np.pi / 180 )    
    rotation_matrix = np.eye(3) + np.sin(angle)*w_matrix + (2*np.sin( angle / 2)**2)*w_matrix.dot(w_matrix)
    return rotation_matrix

def project_vec_to_plane(normal, vec):
    norm_mag = np.sqrt(sum([x**2 for x in normal]))
    projection = vec - ( ( dotproduct(vec, normal) / norm_mag**2 ) * normal )
    return projection

def rotate_cage(sele, angle, sele_vec, align_vec):
	rotation_matrix = rotate_about_vec(angle, align_vec, rad=True)

	translate_mat1 = np.append(rotation_matrix[0], sele_vec[0])
	translate_mat2 = np.append(rotation_matrix[1], sele_vec[1])
	translate_mat3 = np.append(rotation_matrix[2], sele_vec[2])

	translate_mat = np.append(np.append(translate_mat1, np.append(translate_mat2, translate_mat3)), [0,0,0,1])

	cmd.transform_object(sele, translate_mat)
'''