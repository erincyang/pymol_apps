from pymol import cmd

def clrcage( geom="I53" ):
    #parse arche
    geom_list = list(geom)
    arche=geom_list[0]
    comp1=geom_list[1]
    comp2=geom_list[2]
    
    #color
    if ( comp1 == "5" ) or ( comp2 == "5" ):
        cmd.color("firebrick","c5")
    if ( comp1 == "3" ) or ( comp2 == "3" ):
        cmd.color("marine","c3")
    if ( arche == "T" ) and ( comp2 == "3" ):
        cmd.color("teal","c3b")
    if ( comp1 == "2" ) or ( comp2 == "2" ):
        cmd.color("forest","c2")
cmd.extend( "clrcage", clrcage );

def gencage( geom="I53",n="20" ):
    #parse arche
    geom_list = list(geom)
    arche=geom_list[0]
    comp1=geom_list[1]
    comp2=geom_list[2]
    
    #select chains
    cmd.select( "chainA", "chain A and vis" )
    cmd.select( "chainB", "chain B and vis" )
    
    #make cage
    if arche == "T":
        maketet(sel="chainA",name="c{0}".format(comp1),n=int(n))
        if comp2 != None:
            if comp2 == "3":
                maketet(sel="chainB",name="c3b".format(comp2),n=int(n))
            else:
                maketet(sel="chainB",name="c{0}".format(comp2),n=int(n))
    elif arche == "O":
        makeoct(sel="chainA",name="c{0}".format(comp1),n=int(n))
        if comp2 != None:
            makeoct(sel="chainB",name="c{0}".format(comp2),n=int(n))
    elif arche == "I":
        makeicos(sel="chainA",name="c{0}".format(comp1),n=int(n))
        if comp2 != None:
            makeicos(sel="chainB",name="c{0}".format(comp2),n=int(n))
    
    #color comps
    clrcage( geom=geom )
cmd.extend( "gencage", gencage );

def seqres(sel1="vis", offset=1):
    """
        # Written by Spencer Bliven
        # based on http://pymolwiki.org/index.php/Zero_residues
        PURPOSE: renumbers the selected residues sequentially, regardless of gaps
        USAGE: sequential_residues protName    # first residue is 1
        USAGE: sequential_residues protName, 5 # first residue is 5
        EXAMPLE: sequential_residues *
    """
    offset = int(offset)

    # A counter from offset up
    stored.offset = int(offset) - 1
    stored.curr_res = None

    cmd.alter(
        sel1,
        """
if stored.curr_res != int(resi):
    stored.offset=stored.offset + 1
    stored.curr_res=int(resi)
    resi=stored.offset
else:
    resi=stored.offset
""",
    )
    cmd.sort()

# let pymol know about the function
cmd.extend("seqres", seqres)