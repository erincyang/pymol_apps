from pymol import cmd

def meshview(sel="sele"):
	cmd.do( "remove chain_*" )
	for i,c in enumerate(cmd.get_chains(sel)):
		cmd.do( "create chain_%s, %s and chain %s" %(c, sel, c) )
		cmd.do( "hide everything, chain_%s" %(c) )
		cmd.do( "show mesh, chain_%s" %(c) )
cmd.extend( "meshview", meshview );
