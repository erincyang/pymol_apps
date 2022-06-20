#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 11:03:48 2021

@author: quintond
"""
import argparse
import numpy as np

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

def dist_to_plane(plane, point):
    norm_vec = normal_vec(plane)
    X = point - plane[1]
    D = norm_vec.dot(X)
    return D