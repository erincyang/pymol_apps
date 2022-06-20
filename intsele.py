# Usage:
# intsele (distance, bb1, bb2)
#

from pymol import cmd

def intsele( distance=3, bb1="A", bb2="B" ):
	newobj = "%sA_interface" %(distance)
	seleA = "chain %s around %s and chain %s" %(bb1, distance, bb2)
	seleB = "chain %s around %s and chain %s" %(bb2, distance, bb1)
	cmd.select( newobj, "(%s)+(%s)" %(seleA, seleB) )
cmd.extend( "intsele", intsele );

#%s = str
#%i = int