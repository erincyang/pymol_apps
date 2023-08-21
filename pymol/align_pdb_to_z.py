import sys, os
#sys.path.append("/home/erinyang/pymol_apps/pymol")
from optparse import OptionParser
from pymol import cmd

def parse_args():
	#### Command-line flags ####

	parser = OptionParser(description="Align symmetric pdb to Z axis")
	parser.add_option("-i", dest="pdb", default="", help="Input pdb") 
	parser.add_option("-l", dest="pdb_list", default="", help="A file with a list of input pdbs. Optionally, when getting from the database can add _biounitnumber to specify specifc biounits for each pdbid") 
	parser.add_option("-o", dest="output_dir", default="", help="Where should the output files be placed?") 
	parser.add_option("-p", dest="pymol_scripts_dir", default="/home/erinyang/pymol_apps/pymol/", help="Location of Will's util.py pymol module") 
	parser.add_option("-z", dest="alignz", action="store_true", default=False, help="Do you want to extract the shortest chain and align radial symmetry axis with the z-axis?")
	parser.add_option("-e", dest="expected_symmetry", default=None, help="Optional. If you have an expected symmetry, example C3, and the determined symmetry does not match, then put a warning in the output pdb.")

	(options, args) = parser.parse_args()

	return options, args

def clean_pdbname(pdb):
	## replace and strip invalid characters 
	pdb_clean = pdb.split('/')[-1][:-6]

	## replace ' ' with ''
	pdb_clean = pdb_clean.replace(' ', '')

	## replace '=' and ',' with '_'
	pdb_clean = pdb_clean.replace('=', '_')
	pdb_clean = pdb_clean.replace(',', '_')
	print(pdb_clean)
	return pdb_clean

def align_pdb_to_z(pdb, expected_symmetry):
	cmd.do("delete all")
	cmd.do(f"load {pdb}")
	
	pdb_clean = clean_pdbname(pdb)
	cmd.do(f'set_name {pdb_clean}, base_model')
	chains = cmd.get_chains("base_model")

	if len(chains) == 1: 
		total_fasta_length = len(cmd.get_model("base_model").get_residues())
		chain_length = total_fasta_length / int(options.expected_symmetry)
		chains=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','Q','X','Y','Z']
		for n in range(int(options.expected_symmetry)):
			cmd.alter("base_model and resi " + str(n*chain_length+1) + ":" + str((n+1)*chain_length), "chain= '"+str(chains[n]+"'"))
	
	cmd.do('create align_model, base_model')

	cmd.do("run %s" %(options.pymol_scripts_dir + "sym_util.py"))
	cmd.do("run %s" %(options.pymol_scripts_dir + "xyzMath.py"))
	cmd.do(f"aligncx('align_model',{expected_symmetry})")

	## Check direction of N terminus

	cmd.do("trans('align_model', -Vec(com(align_model)))") 
	z_coord = cmd.get_coords('align_model and resi 1 and name CA')[0][2]
	print(z_coord)
	
	if z_coord < 0:
		cmd.do(f"rot('align_model', Ux, 180)")
	
	cmd.do("create chA, chain A and align_model")
	cmd.do(f"save {options.output_dir}/{pdb_clean}_asu.pdb, chA")

#### Begin code ####
if __name__=="__main__":

	options, args = parse_args()

	pdb_list = []
	if options.pdb_list != "":
		pdb_list_file = open(options.pdb_list,"r")
		pdb_list = pdb_list_file.readlines()
		pdb_list_file.close()
	else:
		pdb_list.append(options.pdb)

	for pdb in pdb_list: 
		align_pdb_to_z(pdb, options.expected_symmetry)
