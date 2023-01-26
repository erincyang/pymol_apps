# Usage:
# drawline (objname, color, x1, y1, z1, x2, y2, z2)
#

from pymol import cmd

def drawline( objname="line", color="white", x1=0, y1=0, z1=0, x2=0, y2=0, z2=0 ):
	posString1="[ %s, %s, %s ]" %(x1, y1, z1)
	posString2="[ %s, %s, %s ]" %(x2, y2, z2)
	
	cmd.pseudoatom ("aa", pos=posString1 )
	cmd.pseudoatom ("bb", pos=posString2 )
	
	cmd.distance ("%s" %(objname), "aa", "bb" )
	#cmd.distance ("%s" %(objname), "/aa/////1", "/bb/////1" )
	
	#cmd.set("dash_gap", 0)
	#cmd.set("dash_radius", 0.55)
	#cmd.set("dash_round_ends", 0)
	cmd.set("dash_color", "%s" %(color), "%s" %(objname) )
	cmd.hide("labels", "%s" %(objname))

	cmd.delete("aa")
	cmd.delete("bb")
cmd.extend( "drawline", drawline );

def drawD32( objname="line", color="white", x1=0, y1=0, z1=0, x2=0, y2=0, z2=0 ):
	drawline ( "z_axis", "black", 0,0,75, 0,0,-75 )
	drawline ( "x_axis_1", "black", 0,0,0, 75,0,0 )
	drawline ( "x_axis_2", "black", 0,0,0, 75,0,0 )
	cmd.rotate ( axis="z", angle=120, camera=0, object="x_axis_2", origin=[0,0,0] )
	drawline ( "x_axis_3", "black", 0,0,0, 75,0,0 )
	cmd.rotate ( axis="z", angle=-120, camera=0, object="x_axis_3", origin=[0,0,0]  )

	cgoCircle ( 0,0,0, r=50, w=2, cr=0, cg=0, cb=0 )
cmd.extend( "drawD32", drawD32 );


#%s = str
#%i = int