# Usage:
# marshmallow (objname, color, resol)
#

from pymol import cmd

"""
def marshmallow(objname="marsh", color="white", resol="12" ):
	cmd.set("surface_quality", 1)
	cmd.alter("%s" %(objname), "b=50")
	cmd.alter("%s" %(objname), "q=1")
	cmd.set("gaussian_resolution", "%s" %(resol))
	
	cmd.map_new("map_%s" %(objname), "gaussian", 1, "%s and name ca"  %(objname), 12)
	cmd.isosurface("surf_%s" %(objname), "map_%s" %(objname))
	
	cmd.color("%s" %(color), "surf_%s" %(objname))
	cmd.show("surface", "surf*")
cmd.extend( "marshmallow", marshmallow );
"""

def marshmallow(objname="marsh", color="white", resol="12" ):
	cmd.set("surface_quality", 1)
	cmd.alter("%s" %(objname), "b=50")
	cmd.alter("%s" %(objname), "q=1")
	cmd.set("gaussian_resolution", "%s" %(resol))
	cmd.do("remove (hydro)")
	
	chains = cmd.get_chains("%s" %(objname))
	numchains = len(chains)
	
	for x in range(0, numchains):
		chainname= 'chain ' + chains[x-1]
		mapname= 'map_' + objname + '_' + chains[x-1]
		surfacename= 'surf_' + objname + '_' + chains[x-1]
		selection= objname + ' and ' + chainname + ' and name ca'
		cmd.map_new( mapname, "gaussian", 1, selection, "%s" %(resol) )
		cmd.isosurface( surfacename, mapname )
		x+=1
		
	cmd.color("%s" %(color), "surf_%s*" %(objname))
	cmd.show("surface", "surf*")
cmd.extend( "marshmallow", marshmallow );

def marsh_chains():
	cmd.alter( "chain !", "chain='A'" )
	cmd.alter( "chain #", "chain='B'" )
	cmd.alter( "chain $", "chain='C'" )
	cmd.alter( "chain %", "chain='D'" )
	cmd.alter( "chain &", "chain='E'" )
	cmd.alter( "chain -", "chain='F'" )
	cmd.alter( "chain .", "chain='G'" )
	cmd.alter( "chain <", "chain='H'" )
	cmd.alter( "chain =", "chain='I'" )
	cmd.alter( "chain >", "chain='J'" )
	cmd.alter( "chain ?", "chain='K'" )
	cmd.alter( "chain @", "chain='L'" )
	cmd.alter( "chain ]", "chain='M'" )
	cmd.alter( "chain _", "chain='N'" )
	cmd.alter( "chain {", "chain='O'" )
	cmd.alter( "chain |", "chain='P'" )
	cmd.alter( "chain }", "chain='Q'" )
	cmd.alter( "chain ~", "chain='R'" )
	cmd.alter( "chain 1", "chain='S'" )
	cmd.alter( "chain 2", "chain='T'" )
	cmd.alter( "chain 3", "chain='U'" )
	cmd.alter( "chain 4", "chain='V'" )
	cmd.alter( "chain 5", "chain='W'" )
	cmd.alter( "chain 6", "chain='X'" )
	cmd.alter( "chain 7", "chain='Y'" )
	cmd.alter( "chain 8", "chain='Z'" )
	cmd.alter( "chain 9", "chain='A'" )
	cmd.alter( "chain 0", "chain='B'" )
cmd.extend( "marsh_chains", marsh_chains );

#%s = str
#%i = int