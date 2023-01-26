from pymol.cgo import *
from pymol import cmd
from pymol.vfont import plain

w = 0.2 # cylinder width 
l = 150.0 # cylinder length
h = 1.0 # cone hight
d = w * 1.618 # cone base diameter
 
obj = [CYLINDER, 0.0, 0.0, 0.0,   l, 0.0, 0.0, w, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
       CYLINDER, 0.0, 0.0, 0.0, 0.0,   l, 0.0, w, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
       CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   l, w, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
       CONE,   l, 0.0, 0.0, h+l, 0.0, 0.0, d, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 
       CONE, 0.0,   l, 0.0, 0.0, h+l, 0.0, d, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 
       CONE, 0.0, 0.0,   l, 0.0, 0.0, h+l, d, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0]

cyl_text(obj,plain,[-5.,-5.,-1],'Origin',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
cyl_text(obj,plain,[10.,0.,0.],'X',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
cyl_text(obj,plain,[0.,10.,0.],'Y',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
cyl_text(obj,plain,[0.,0.,10.],'Z',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])

cmd.load_cgo(obj, 'axes')
