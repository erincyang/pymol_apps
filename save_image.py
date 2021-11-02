# Usage:
# save_image(imgname)
#

from pymol import cmd

def save_image( imgname="img" ):
	cmd.do( "ray 1600, 1200" )
	cmd.do( "png /home/yhsia/Desktop/pymol_images/%s.png" %(imgname) )
cmd.extend( "save_image", save_image );

"""
#rotate
rotate axis=z, angle=15, camera=0, origin=[0,0,0], object=objname

#translate
translate vector=[0,0,15], camera=0, object=objname
"""
