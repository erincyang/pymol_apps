#!/usr/bin/python

__author__ = "Jacob Bale"
__version__ = "1.0.4"
__email__ = "balej@uw.edu"

#### Script for cleaning up a pdb and aligning cyclic complexes with the z axis (if desired) ####
#### CAUTION: For use with crystal structures and homo-oligomeric protein assemblies only.  Logic has not been tested with NMR or EM structures.
#### Ex. Command lines:

# From biounits -> To user-specified path: single pdbid, single biounit: 
	#pymol -ck path_to_script/auto_prep_symmetric_pdb.py -- -i 1qah -o ~/foo/ -b 2 -g -z -e C3

# From biounits => To user-specified path: single pdbid, all biounits: 
	#pymol -ck path_to_script/auto_prep_symmetric_pdb.py -- -i 1qah -o ~/foo/ -g -z -e C3

# From biounits -> To database: single pdbid, single biounit: 
	#pymol -ck path_to_script/auto_prep_symmetric_pdb.py -- -i 1qah -b 2 -a -g -z -e C3

# From biounits => To database: single pdbid, all biounits: 
	#pymol -ck path_to_script/auto_prep_symmetric_pdb.py -- -i 1qah -a -g -z -e C3

# From somewhere else -> To user-specified path: single pdb: 
	#pymol -ck path_to_script/auto_prep_symmetric_pdb.py -- -i ~/foo/1qah.pdb3.gz -o ~/foo/ -z -e C3
	#or
	#pymol -ck path_to_script/auto_prep_symmetric_pdb.py -- -i ~/foo/1qah_3.pdb -o ~/foo/ -z -e C3

# From somewhere else -> To user-specified path: multiple pdbs: 
	#pymol -ck path_to_script/auto_prep_symmetric_pdb.py -- -l YOURPATHSTOPDBS.list -o ~/foo/ -z -e C3

# From somewhere else -> To database: single pdb: 
	#pymol -ck path_to_script/auto_prep_symmetric_pdb.py -- -i ~/foo/1qah.pdb3.gz -a -z -e C3

# From somewhere else -> To database: single pdb: 
	#pymol -ck path_to_script/auto_prep_symmetric_pdb.py -- -i ~/foo/1qah.pdb3.gz -a -z -e C3

# From biounits -> To user-specified path: multiple pdbids, specific biounit numbers specified with _#: 
	#pymol -ck path_to_script/auto_prep_symmetric_pdb.py -- -l PISA_HOMODIMERS.list -o ~/foo/ -g -z -e C2 

# From biounits -> To user-specified path: multiple pdbids, same biounit number for all: 
	#pymol -ck path_to_script/auto_prep_symmetric_pdb.py -- -l PISA_HOMODIMERS.list -o ~/foo/ -b 3 -g -z -e C2 

# From biounits -> To User-specified path: multiple pdbids, all biounits: 
	#pymol -ck path_to_script/auto_prep_symmetric_pdb.py -- -l PISA_HOMODIMERS.list -o ~/foo/ -g -z -e C2 

# From biounits -> To database: multiple pdbids, specific biounit numbers specified with _#: 
	#pymol -ck path_to_script/auto_prep_symmetric_pdb.py -- -l PISA_HOMODIMERS.list -a -g -z -e C2 

# From biounits -> To database: multiple pdbids, same biounit number for all: 
	#pymol -ck path_to_script/auto_prep_symmetric_pdb.py -- -l PISA_HOMODIMERS.list -b 3 -a -g -z -e C2 

# From biounits -> To database: multiple pdbids, all biounits: 
	#pymol -ck path_to_script/auto_prep_symmetric_pdb.py -- -l PISA_HOMODIMERS.list -a -g -z -e C2 

import string, sys, os, inspect, datetime, subprocess, re
newpath = os.path.dirname(inspect.getfile(inspect.currentframe())) # script directory
if not newpath in sys.path: sys.path.append(newpath)
from optparse import OptionParser
try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3
import numpy as np

#### Command-line flags ####

parser = OptionParser()
#### Possible Inputfiles ####
parser.add_option("-i", dest="pdb", default="", help="Input pdb") 
parser.add_option("-l", dest="pdb_list", default="", help="A file with a list of input pdbs. Optionally, when getting from the database can add _biounitnumber to specify specifc biounits for each pdbid") 
parser.add_option("-o", dest="output_dir", default="", help="Where should the output files be placed?") 
parser.add_option("-b", dest="biounit_number", default="", help="Optional Argument If For single_biounit_mode: Which biounit should be used? If not given then will try all available.") 
parser.add_option("-p", dest="pymol_scripts_dir", default="%s" %(os.getenv("HOME")+"/pymol_apps/pymol/"), help="Location of Will's util.py pymol module") 
parser.add_option("-a", dest="add_to_database", action="store_true", default=False, help="Do you want to add this to the oligomer database?")
parser.add_option("-g", dest="get_from_database", action="store_true", default=False, help="Do you want to get the structure from the shared biounits database?")
parser.add_option("-z", dest="alignz", action="store_true", default=False, help="Do you want to extract the shortest chain and align radial symmetry axis with the z-axis?")
parser.add_option("-v", dest="view_in_pymol", action="store_true", default=False, help="If executed without the -c option this will enable you to view the results in pymol")
parser.add_option("-1", dest="max_distance_error", default=0.3, help="Threshold for the error in the distance for an acceptable alignement")
parser.add_option("-2", dest="max_angle_error", default=0.3, help="Threshold for the error in the angle for an acceptable alignement")
parser.add_option("-e", dest="expected_symmetry", default=None, help="Optional. If you have an expected symmetry, example C3, and the determined symmetry does not match, then put a warning in the output pdb.")
(options, args) = parser.parse_args()

pdb_list = []
if options.pdb_list != "":
	pdb_list_file = open(options.pdb_list,"r")
	pdb_list = pdb_list_file.readlines()
	pdb_list_file.close()
else:
	pdb_list.append(options.pdb)

single_biounit_mode = False
biounit_input_pattern = re.compile('(\w{4})\.pdb([1-9])')

database_array = None
new_rows = []
database_headers = []
refids = []
if options.add_to_database:
	database_file = open(os.getenv("HOME")+"/databases/symmetric_scaffolds/symmetric_scaffolds.dat","r")
	database_string = database_file.read()
	database_file.close()

	database_array =  np.genfromtxt(StringIO(database_string), names=True, dtype=None, delimiter="\t")	
	database_headers = list(database_array.dtype.names)
	refids = database_array['refid']

# For the sake of speed, this assumes that all fields are already present in the database (ie, no new columns need to be added).
def add_entry(database_headers,refids,my_dict):
	new_row = []
	version = 1
	while (my_dict['refid']+"."+str(version) in refids):
		version+=1
	my_dict['refid'] = my_dict['refid']+"."+str(version)
	for header in database_headers:
		try:
			new_row.append(str(my_dict[header]))
		except KeyError:
			new_row.append("-")
	return new_row

for pdb in pdb_list:

	another_biounit = True

	input_pdb_path, input_pdb_name, pdbid_bunit, output_pdb_path, output_pdbname = "", "", "", "", ""
	log_string1, log_string2, log_string3, warning_string1 = "", "", "", ""

	biounit_number=1
	if options.biounit_number != "":
		biounit_number = int(options.biounit_number)
		single_biounit_mode = True

	if options.get_from_database:
		if len(pdb.split("_")) != 1:
			single_biounit_mode = True
			biounit_number = int(pdb.split("_")[1].strip())
			input_pdbname = pdb.split("_")[0].strip()+".pdb"+str(biounit_number)+".gz"
			pdbid_bunit = pdb.split("_")[0].strip()+"_"+str(biounit_number)
			output_pdbname = pdbid_bunit + ".pdb"
		else:
			input_pdbname = pdb.strip() + ".pdb" + str(biounit_number) + ".gz"
			pdbid_bunit = pdb.strip() + "_" + str(biounit_number)
			output_pdbname = pdbid_bunit + ".pdb"
		input_pdb_path = os.getenv("HOME")+"/databases/biounits/"+input_pdbname[1:3]+"/"	
	elif biounit_input_pattern.match(os.path.basename(pdb)):
		match_obj = biounit_input_pattern.search(pdb)
		if os.path.dirname(pdb) != "":
			input_pdb_path = os.path.dirname(pdb) + "/"
		input_pdbname = os.path.basename(pdb)
		output_pdbname = match_obj.group(1)+"_"+ match_obj.group(2)+".pdb"
		pdbid_bunit = match_obj.group(1) + "_" + match_obj.group(2)
		single_biounit_mode = True

	else:
		if os.path.dirname(pdb) != "":
			input_pdb_path = os.path.dirname(pdb) + "/"
		input_pdbname = os.path.basename(pdb)
		output_pdbname = input_pdbname.replace(".gz","")
		pdbid_bunit = pdb.strip() + "_" + str(biounit_number)
		single_biounit_mode = True

	if options.add_to_database:
		output_pdb_path = os.getenv("HOME")+"/databases/symmetric_scaffolds/"+input_pdbname[1:3]+"/"
		if (not os.path.exists(output_pdb_path)):
			print("Making new output directory: " + output_pdb_path)
			subprocess.check_output("mkdir -p " + output_pdb_path,stderr=subprocess.STDOUT,shell=True)
		if os.path.exists(output_pdb_path+output_pdbname):
			subprocess.check_output("rm " + output_pdb_path + output_pdbname,stderr=subprocess.STDOUT,shell=True) 
		if os.path.exists(input_pdb_path + input_pdbname):
			subprocess.check_output("cp " + input_pdb_path + input_pdbname + " " + output_pdb_path + output_pdbname + ".gz",stderr=subprocess.STDOUT,shell=True) 
			subprocess.check_output("gunzip " + output_pdb_path + output_pdbname + ".gz",stderr=subprocess.STDOUT,shell=True) 
			log_string1 = "PROCESSING   InputPDB: " + input_pdb_path + input_pdbname
			subprocess.check_output("cp " + output_pdb_path + output_pdbname + " " + output_pdb_path + output_pdbname.replace(".pdb","") + "_ref.pdb",stderr=subprocess.STDOUT,shell=True) 
		else:
			print(input_pdbname + " reports failure: File could not be found.")
			continue
	else:
		if  options.output_dir != "":
			output_pdb_path = options.output_dir.rstrip("/") + "/"
		if (pdb[-3:] == ".gz") or options.get_from_database:
			subprocess.check_output("cp " + input_pdb_path + input_pdbname + " " + output_pdb_path + output_pdbname + ".gz",stderr=subprocess.STDOUT,shell=True) 
			if os.path.exists(output_pdb_path+output_pdbname):
				subprocess.check_output("rm " + output_pdb_path + output_pdbname,stderr=subprocess.STDOUT,shell=True) 
			subprocess.check_output("gunzip " + output_pdb_path + output_pdbname + ".gz",stderr=subprocess.STDOUT,shell=True) 
			log_string1 = "PROCESSING   InputPDB: " + input_pdb_path + input_pdbname 
			subprocess.check_output("cp " + output_pdb_path + output_pdbname + " " + output_pdb_path + output_pdbname.replace(".pdb","") + "_ref.pdb",stderr=subprocess.STDOUT,shell=True) 
		else:
			subprocess.check_output("cp " + input_pdb_path + input_pdbname + " " + output_pdb_path + output_pdbname,stderr=subprocess.STDOUT,shell=True) 
			subprocess.check_output("cp " + input_pdb_path + input_pdbname + " " + output_pdb_path + output_pdbname.replace(".pdb","") + "_ref.pdb",stderr=subprocess.STDOUT,shell=True) 
			log_string1 = "PROCESSING   InputPDB: " + input_pdb_path + input_pdbname 

	while another_biounit: 

		entry_dict = {}
		entry_dict['refid'] = pdbid_bunit
		entry_dict['pdbid'] = pdbid_bunit.split("_")[0]
		entry_dict['bunit'] = pdbid_bunit.split("_")[1]

		# Open the input pdb.
		input_pdbfile = open(output_pdb_path + output_pdbname, "r")
		input_pdbfile_list = input_pdbfile.readlines()
		input_pdbfile.close()

		# In the pdb for alignment: replace selenomethionines with methionines and remove alternate conformations other than A (do we need to adjust bond lengths? is there a way to do this properly in Pymol? Removing alt confs in pymol doesn't seem to remove them from behind the scenes...)

		header_list = []
		newline_list = []
		header = True
		selenomet = False
		altscconfs = False
		ex_chains = False
		ex_models = False
		construct_length = 0
		for line in input_pdbfile_list:
			newline = ""
			if line[0:4] == "ATOM" or (line[0:6] == "HETATM" and line[17:20] == "MSE"):
				# Replace selenomets.
				if line[0:6] == "HETATM" and line[17:20] == "MSE":
					if line[12:14] == "SE":
							newline = "ATOM" + line[4:11] + "  SD  MET " + line[21:76] + "SD" + line[78:].rstrip()
					else:
							newline = "ATOM" + line[4:16] + " MET " + line[21:].rstrip()			
					selenomet = True
					newline_list.append(newline)
				# Skip alternative conformations.
				elif (line[16] != " "):
					altscconfs = True
					if line[16] != "A":
						continue
					else:
						newline = line[:16] + " " + line[17:]
						newline_list.append(newline)
				else:
					newline = line.rstrip()
					newline_list.append(newline)
				header = False
			elif (line[:6] == "HETATM") or (line[:3] == "TER") or (line[:3] == "END") or (line[:5] == "MODEL") or (line[:6] == "ENDMDL"):
				newline = line.rstrip()
				header = False
				newline_list.append(newline)
			elif header:
				header_list.append(line.rstrip())
				if line[:6] == "SEQRES":
					construct_length = int(line[13:17])

		newlines = "\n".join(newline_list)
		header_string = "\n".join(header_list)
		# Save the new pdb with selenomethionines replaced and alt sc confs removed.
		output_pdbfile = open(output_pdb_path + output_pdbname, "w")
		output_pdbfile.write(header_string + "\n" + newlines)
		output_pdbfile.close()

		# Keep track of missing or unrecognized residues. If only one chain per model then convert models to chains (indicates assembly composed of copies of the asymmetric unit: most likely cyclic symmetry), if multiple models with more than one chain per model then only keep the first model for evaluation/alignment and make note in header (assuming one is only dealing with homooligomers, this indicates higher symmetry which we are currently not able to assess/align automatically), and pick out the chain with the fewest atoms and reorder to make that chain A. Save this as the "clean" pdb, which should be easy to align with design models for visual comparisons to the native with ligands/waters/etc intact, but no need to split_states or make copies of chain A to align with the design model.

		#Create an initial reference pdb object and an initial alignment pdb object (remove hetams from pdb for alignment).

		object1 = output_pdbname.replace(".pdb","")
		cmd.do("load %s" %(output_pdb_path + output_pdbname))
		cmd.do("load %s" %(output_pdb_path + object1 + "_ref.pdb"))
		cmd.do("set_name %s, base_model" %(object1+"_ref"))
		chains = cmd.get_chains(object1)
		states = cmd.count_states(object1)
		#If only one chain per model then convert models to chains (indicates assembly composed of copies of the asymmetric unit: most likely cyclic symmetry). Keep the first 6 models/chains.
		model_warning = ""
		chain_warning = ""
		chain_model_warning = ""
		single_chain_per_model = False
		my_misc_dict = { 'tmp_list' : [] , 'align_resi_set' : set() }
		input_output_chain_array = []
		if states > 1 and len(chains) == 1:
			single_chain_per_model = True
			cmd.do("split_states %s" %(object1))
			cmd.do("split_states base_model")
			merge_list = []
			base_merge_list = []
			for state in range(1,states+1):
				if state > 6:
					model_warning = "\nPROCESSING   Model_Warning: "+str(states)+" models in the biological unit with 1 chain each. Only the first 6 were considered. Possible higher order symmetry exists than detected."   
					ex_models = True
					#try:
					#	entry_dict['notes']+=" "+model_warning[14:]
					#except KeyError:
					#	entry_dict['notes']=model_warning[14:] 
					break
				chain_ord = 64 + state
				chain_id = chr(chain_ord)
				input_output_chain_array.append((chains[0],chain_id))
				state_suffix = "_%04d" %(state)
				cmd.do("set_name %s, chain_%s" %(object1+state_suffix, chain_id))
				cmd.do("set_name %s, base_chain_%s" %("base_model"+state_suffix, chain_id))
				cmd.do("alter chain_%s, chain='%s'" %(chain_id, chain_id)) 
				cmd.do("alter base_chain_%s, chain='%s'" %(chain_id, chain_id)) 
				merge_list.append("chain_" + chain_id )
				base_merge_list.append("base_chain_" + chain_id )
			merge_string = " + ".join(merge_list)
			base_merge_string = " + ".join(base_merge_list)
			cmd.do("create base_ref_model, (%s)" %(base_merge_string))
			cmd.do("create align_model, (%s)" %(merge_string))
			cmd.do("remove (align_model and hetatm)")
			cmd.do("remove (align_model and resn a+g+c+t+u)")

		#If multiple models with more than one chain per model then convert to one model with A x B chains (where A is the number of models and B is the number of chains per model).
		elif (states > 1) and (len(chains) > 1):
			chain_model_warning = "\nPROCESSING   Chain_Model_Warning: "+str(states)+" models with "+str(len(chains)) +" chains each in the biological unit. Possible higher order symmetry exists than detected."   
			ex_models = True
			merge_list = []
			base_merge_list = []
			#try:
			#	entry_dict['notes']+=" "+chain_model_warning[14:]
			#except KeyError:
			#	entry_dict['notes']=chain_model_warning[14:] 
			cmd.do("split_states %s" %(object1))
			cmd.do("split_states base_model")
			failure = False
			#Catch potential cases with different numbers of chains per model.  This code is not equiped to handle such cases.
			for state in range(1,states+1):
				state_suffix = "_%04d" %(state)
				if len(chains) != len(cmd.get_chains("%s" %(object1+state_suffix))):
					print(input_pdbname + " reports failure: The number of chains in the different states do not match.")
					biounit_number+=1
					if single_biounit_mode:
						another_biounit = False
					elif not os.path.exists(input_pdb_path+pdb.strip()+".pdb"+str(biounit_number)+".gz"):
						another_biounit = False
					else:
						input_pdbname = pdb.strip() + ".pdb" + str(biounit_number) + ".gz"
						pdbid_bunit = pdb.strip() + "_" + str(biounit_number)
						output_pdbname = pdbid_bunit + ".pdb"
						if os.path.exists(output_pdb_path+output_pdbname):
							subprocess.check_output("rm " + output_pdb_path + output_pdbname,stderr=subprocess.STDOUT,shell=True) 
						subprocess.check_output("cp " + input_pdb_path + input_pdbname + " " + output_pdb_path + output_pdbname + ".gz",stderr=subprocess.STDOUT,shell=True) 
						subprocess.check_output("gunzip " + output_pdb_path + output_pdbname + ".gz",stderr=subprocess.STDOUT,shell=True) 
						log_string1 = "PROCESSING   InputPDB: " + input_pdb_path + input_pdbname
						subprocess.check_output("cp " + output_pdb_path + output_pdbname + " " + output_pdb_path + output_pdbname.replace(".pdb","") + "_ref.pdb",stderr=subprocess.STDOUT,shell=True) 
					cmd.reinitialize()
					failure = True
					break
			if failure:
				continue

			#Convert A models with B chains into 1 model with AxB chains.
			too_many_chains = False
			chain_ord = 64
			for state in range(1,states+1):
				state_suffix = "_%04d" %(state)
				state_chains = cmd.get_chains("%s" %(object1+state_suffix))
				for ch in state_chains:
					chain_ord+=1
					chain_id = chr(chain_ord)
					if chain_ord-64 > 6:
						chain_warning = "\nPROCESSING   Chain_Warning: "+str(len(chains))+" chains per model in the biological unit. Only the first 6 were considered. Possible higher order symmetry exists than detected."   
						ex_chains = True
						#try:
						#	entry_dict['notes']+=" "+chain_warning[14:]
						#except KeyError:
						#	entry_dict['notes']=chain_warning[14:] 
						too_many_chains = True
						break
					input_output_chain_array.append((state_chains[0],chain_id)) #Does this make sense for these cases with multiple chains and multiple models?
					cmd.do("create chain_%s, %s and chain %s" %(chain_id, object1+state_suffix, ch))
					cmd.do("create base_chain_%s, %s and chain %s" %(chain_id, "base_model"+state_suffix, ch))
					cmd.do("alter chain_%s, chain='%s'" %(chain_id, chain_id)) 
					cmd.do("alter base_chain_%s, chain='%s'" %(chain_id, chain_id)) 
					merge_list.append("chain_" + chain_id )
					base_merge_list.append("base_chain_" + chain_id )
				if too_many_chains:
					break
			merge_string = " + ".join(merge_list)
			base_merge_string = " + ".join(base_merge_list)
			cmd.do("create base_ref_model, (%s)" %(base_merge_string))
			cmd.do("create align_model, (%s)" %(merge_string))
			cmd.do("remove (align_model and hetatm)")
			cmd.do("remove (align_model and resn a+g+c+t+u)")

		#If one model and multiple chains then take first 6 chains.
		elif states == 1 and len(chains) > 1:
			cmd.do("set_name base_model, base_ref_model")
			cmd.do("set_name %s, align_model" %(object1))
			cmd.do("remove (align_model and hetatm)")
			cmd.do("remove (align_model and resn a+g+c+t+u)")
			chain_count = 0
			for ch in chains:
				chain_count+=1
				if chain_count > 6:
					chain_warning = "\nPROCESSING   Chain_Warning: "+str(len(chains))+" chains in the biological unit. Only the first 6 were considered. Possible higher order symmetry exists than detected."   
					ex_chains = True
					#try:
					#	entry_dict['notes']+=" "+chain_warning[14:]
					#except KeyError:
					#	entry_dict['notes']=chain_warning[14:] 
					cmd.do("remove chain %s" %(ch))
		else:
			print(input_pdbname + " reports failure: Only 1 chain 1 model found.")
			biounit_number+=1
			if single_biounit_mode:
				another_biounit = False
			elif not os.path.exists(input_pdb_path+pdb.strip()+".pdb"+str(biounit_number)+".gz"):
				another_biounit = False
			else:
				input_pdbname = pdb.strip() + ".pdb" + str(biounit_number) + ".gz"
				pdbid_bunit = pdb.strip() + "_" + str(biounit_number)
				output_pdbname = pdbid_bunit + ".pdb"
				if os.path.exists(output_pdb_path+output_pdbname):
					subprocess.check_output("rm " + output_pdb_path + output_pdbname,stderr=subprocess.STDOUT,shell=True) 
				subprocess.check_output("cp " + input_pdb_path + input_pdbname + " " + output_pdb_path + output_pdbname + ".gz",stderr=subprocess.STDOUT,shell=True) 
				subprocess.check_output("gunzip " + output_pdb_path + output_pdbname + ".gz",stderr=subprocess.STDOUT,shell=True) 
				log_string1 = "PROCESSING   InputPDB: " + input_pdb_path + input_pdbname
				subprocess.check_output("cp " + output_pdb_path + output_pdbname + " " + output_pdb_path + output_pdbname.replace(".pdb","") + "_ref.pdb",stderr=subprocess.STDOUT,shell=True) 
			cmd.reinitialize()
			continue
		# Save the current objects and reinitialize to get rid of extraneous structures/info.
		cmd.do("save %s, base_ref_model" %(output_pdb_path + object1 + "_ref.pdb"))
		cmd.do("save %s, align_model" %(output_pdb_path + object1 +"_align_model.pdb"))
		cmd.reinitialize()
		cmd.do("load %s" %(output_pdb_path + object1 + "_align_model.pdb"))
		cmd.do("set_name %s, align_model" %(object1 +"_align_model"))
		cmd.do("create target_model, align_model")

		# Record rmsds after refinement using super for each pairwise combination of chains.
		algn = None
		rmsd_array = []
		chains = cmd.get_chains("align_model")
		conserved = None
		complete_consensus = True
		for index_1 in range(0,len(chains)):
			other_chains = chains[:]
			chain_1 = other_chains.pop(index_1)
			if index_1 == 0:
				cmd.iterate("(align_model and chain %s and name ca)" %(chain_1),"tmp_list.append((resi,resn))",space=my_misc_dict)
			index_1_nresi = len(my_misc_dict['tmp_list'])
			my_misc_dict['tmp_list'] = [] 
			for index_2 in range(0,len(other_chains)):
				if index_2 == 0:
					rmsd_array.append([float(cmd.super("(align_model and chain %s)" %(chains[index_1]), "(target_model and chain %s)" %(other_chains[index_2]), object="algn")[0])])
				else:
					rmsd_array[index_1].append(float(cmd.super("(align_model and chain %s)" %(chains[index_1]), "(target_model and chain %s)" %(other_chains[index_2]), object="algn")[0]))

				# Get the residues that had atoms matched in each superposition and determine the intersection. Call this the "consensus sequence". Only need to do with one target_chain.
				if index_1 == 0:
					cmd.iterate("(target_model and chain %s and algn)" %(other_chains[index_2]),"align_resi_set.add(resi)",space=my_misc_dict)
					if index_2 == 0:
						conserved = my_misc_dict['align_resi_set'].copy()
						len1 = len(conserved)
						if len1 != index_1_nresi:
							complete_consensus = False
					else:
						conserved = set.intersection(conserved,my_misc_dict['align_resi_set'])
						if len1 != len(conserved):
							complete_consensus = False 
					my_misc_dict['align_resi_set'].clear()

		# Select the chain with the lowest average RMSD when aligned to all other chains.
		chain_index = 0
		selected_chain_index = 0
		min_average_rmsd = None
		for target_chain_list in rmsd_array:
			total_rmsd = 0
			for source_chain_rmsd in target_chain_list:
				total_rmsd+=source_chain_rmsd
			if chain_index == 0:
				min_average_rmsd = round(total_rmsd/len(target_chain_list),3)
			elif round(total_rmsd/len(target_chain_list),3) < min_average_rmsd:
				min_average_rmsd = round(total_rmsd/len(target_chain_list),3)
				selected_chain_index = chain_index
			chain_index+=1
		selected_chain = chains[selected_chain_index]

		if options.add_to_database:
			entry_dict['source_chain'] = selected_chain
			entry_dict['ave_intchain_rmsd'] = "%.3f" % min_average_rmsd
			if complete_consensus:
				entry_dict['comp_consen'] = True
			else:
				entry_dict['comp_consen'] = False

		# Generated a pymol selection string for the "consensus sequence"
		conserved_list = list(conserved)
		conserved_int_list = []
		for resi in conserved_list:
			conserved_int_list.append(int(re.sub(r'[a-zA-Z]','',"%s" % resi))) #Account for insertion codes.
		conserved_int_list.sort()
		previous_resi = None
		resi_string = "\nPROCESSING   Consensus:"
		if len(conserved_int_list) > 0:
			previous_resi = conserved_int_list[0]
			endcapped = False
			resi_string+=" select consensus, "+object1+"_ref and resi "+ str(previous_resi)
			for resi in conserved_int_list:
				endcapped = False
				if resi-previous_resi > 1:
					resi_string+="-"+str(previous_resi)+"+"+str(resi)
					endcapped = True
				previous_resi = resi
			if not endcapped:
				resi_string+="-"+str(conserved_int_list[-1])
		else:
			resi_string+=" No Consensus Positions Detected! Possibly Due To Inconsistent Chain Numbering In The Input PDB."

		# Reinitialize to get rid of extraneous structures/info and load back in the two models in progress.
		cmd.reinitialize()
		cmd.do("load %s" %(output_pdb_path + object1 + "_align_model.pdb"))
		cmd.do("set_name %s, align_model" %(object1 +"_align_model"))
		cmd.do("load %s" %(output_pdb_path + object1 + "_ref.pdb"))
		cmd.do("set_name %s, base_ref_model" %(object1 + "_ref"))

		# Reorder/rename chains to make the selected chain be chain A.

		chain_count = 0
		new_first_chain = chains.pop(selected_chain_index)
		chains.insert(0,new_first_chain)
		merge_list = []
		base_merge_list = []
		for ch in chains:
			chain_count+=1
			chain_ord = 64 + chain_count
			chain_id = chr(chain_ord)
			cmd.do("create chain_%s, align_model and chain %s" %(chain_id,ch))
			cmd.do("create base_chain_%s, base_ref_model and chain %s" %(chain_id,ch))
			cmd.do("alter chain_%s, chain='%s'" %(chain_id, chain_id)) 
			cmd.do("alter base_chain_%s, chain='%s'" %(chain_id, chain_id)) 
			if not single_chain_per_model:
				input_output_chain_array.append((ch,chain_id))
			merge_list.append("chain_%s" %(chain_id))
			base_merge_list.append("base_chain_%s" %(chain_id))

		merge_string = " + ".join(merge_list)
		base_merge_string = " + ".join(base_merge_list)
		merge_list = []
		cmd.do("create align_model, (%s)" %(merge_string))
		cmd.do("create base_ref_model, (%s)" %(base_merge_string))

		# Save the base_model as a reference structure and add in the header information plus min chain, min_chain_length, number atoms fewer than longest chain, chain mapping.
		cmd.do("save %s, base_ref_model" %(output_pdb_path + object1 + "_ref.pdb"))

		inputfile = open(output_pdb_path + object1 + "_ref.pdb","r")
		inputfile_list = inputfile.readlines()
		inputfile.close()
		subprocess.check_output("rm " + output_pdb_path + object1 + "_ref.pdb",stderr=subprocess.STDOUT,shell=True)

		header = True
		chain_list = []
		for line in inputfile_list:
			if (line[:4] == "ATOM") or (line[:6] == "HETATM"):
				header = False
				chain_list.append(line.rstrip())  
			elif not header:
				chain_list.append(line.rstrip())  

		chain_string = "\n".join(chain_list)

		chain_mapping_list = ["PROCESSING   Chain_Mapping:  InputPDBChain  NewPDBChain"]
		for ch in input_output_chain_array:
			chain_mapping_list.append("PROCESSING   Chain_Mapping:  " + ch[0].center(14) + ch[1].center(14))
		chain_mapping_string = "\n".join(chain_mapping_list)

		log_string2 = "PROCESSING   PDB_Parsing_Summary: A new reference version of " + input_pdb_path + input_pdbname + " was generated in which the chain with the lowest average RMSD\nPROCESSING   PDB_Parsing_Summary: to all other chains (" + str(min_average_rmsd) + "), chain " + selected_chain + ", was set to be chain A and all other ids renamed alphabetically after A.\nPROCESSING   PDB_Parsing_Summary: The new pdb was saved as " + output_pdb_path + object1 + "_ref.pdb" + chain_warning + model_warning + chain_model_warning 

		log_string2 = log_string2 + "\n" + chain_mapping_string

		now = datetime.datetime.now()	

		comment_string = "PROCESSING   Date: " + str(now.month) + "/" + str(now.day) + "/" + str(now.year) + "\n" + log_string1 + "\n" + log_string2
		clean_pdb_outfile = open(output_pdb_path + object1 + "_ref.pdb","w")
		clean_pdb_outfile.write(comment_string + "\n" + header_string + "\n" + chain_string)
		clean_pdb_outfile.close() 

		if selenomet and altscconfs:
			comment_string+="\nPROCESSING   PDB_Parsing_Summary: Selenomethionines were replaced with MET, alternative side-chains other than A were removed, and hetams removed."
		elif selenomet and not altscconfs:
			comment_string+="\nPROCESSING   PDB_Parsing_Summary: Selenomethionines were replaced with MET and hetams removed."
		elif not selenomet and altscconfs:
			comment_string+="\nPROCESSING   PDB_Parsing_Summary: Alternative side-chains other than A and hetams were removed."
		else:
			pass

		# Remove residues that have missing mainchain atoms.

		missing_mc_resi_list = []
		atoms = cmd.get_model("chain_A and name N+CA+C+O").atom
		mainchain_count = 0
		prev_resi = atoms[0].resi
		for a in atoms:
			if a.resi != prev_resi:
				if mainchain_count < 4:
					missing_mc_resi_list.append(prev_resi)
				mainchain_count=0
			prev_resi = a.resi
			mainchain_count+=1
		for resi in missing_mc_resi_list:
			cmd.do("remove chain_A and resi %s" %(resi))
		cmd.do("save %s, chain_A" %(output_pdb_path + object1 +"_chain_A_align_model.pdb"))
		cmd.do("save %s, align_model" %(output_pdb_path + object1 +"_align_model.pdb"))
		cmd.reinitialize()

		# Renumber and map missing residues. Easiest to do manually.

		input_pdbfile = open(output_pdb_path + object1 +"_chain_A_align_model.pdb", "r")
		input_pdbfile_list = input_pdbfile.readlines()
		input_pdbfile.close()

		rescount = 0
		atomcount = 0
		previous_resnum = None
		previous_resi = None
		missing_resi_list = []
		newline_list = []
		for line in input_pdbfile_list:
			newline = ""
			if line[:4] == "ATOM":
				atomcount+=1
				if previous_resnum == None:
					rescount = 1
					if int(line[22:26]) > 1:
						missing_resi_list.append(("ND",line[22:26].strip(),str(0),str(1)))
				elif line[22:27] != previous_resi: #Account for insertion codes.
					rescount+=1
					if (int(line[22:26]) - previous_resnum) > 1: #May not be correct in all cases.
						missing_resi_list.append((str(previous_resnum),line[22:26].strip(),str(rescount-1),str(rescount)))
				newline = "ATOM " + str(atomcount).rjust(6) + line[11:22] + str(rescount).rjust(4) + line[26:] 
				previous_resnum = int(line[22:26])
				previous_resi = line[22:27] #Account for insertion codes.
				newline_list.append(newline.rstrip())
			elif line[:3] == "TER" or line[:3] == "END":
				newline_list.append(line[:3])
		if construct_length - previous_resnum > 0:
			missing_resi_list.append((str(previous_resnum),str(construct_length),str(rescount),"NA"))

		if options.add_to_database:
			entry_dict['chain_len'] = rescount
			entry_dict['construct_len'] = construct_length 

		newlines = "\n".join(newline_list)

		# Save the new renumbered pdb.
		output_pdbfile = open(output_pdb_path + object1 +"_chain_A_align_model.pdb", "w")
		output_pdbfile.write(newlines)
		output_pdbfile.close()

		# Generate the missing resi string for the output file headers.

		missing_resi_header = "PROCESSING   Missing_or_Unrecognized_Residues:   RefPDB    AlignedPDB"
		missing_resi_strings_list = [missing_resi_header]
		missing_resi_string = ""
		n_missing_resis = 0
		if missing_resi_list != []:
			for i in missing_resi_list:
				missing_resi_strings_list.append("PROCESSING   Missing_or_Unrecognized_Residues:  "+ str(i[0]+"-"+i[1]).center(10) +" " + str(i[2]+"-"+i[3]).center(10))
		n_missing_resis = (construct_length - rescount)

		if options.add_to_database:
			entry_dict['n_missing_resis'] = n_missing_resis

		if missing_resi_list != []:
			missing_resi_string = "\n"+"\n".join(missing_resi_strings_list)
			comment_string+=missing_resi_string

		# Iterate through chain A-F: Use pymol to align copies of the selected chain onto all others and Will Sheffler's cx alignment script to attempt to align each according to the current number of subunits being considered.
		# Note: Currently only considering Cx symmetries, but the actually assembly could be higher order.  

		if options.alignz:
			results_list = []
			tmp = None
			cmd.do("run %s" %(options.pymol_scripts_dir + "sym_util.py"))
			cmd.do("run %s" %(options.pymol_scripts_dir + "xyzMath.py"))
			for nsubs in range(2,len(chains)+1):
				merge_list = []
				cmd.do("load %s" %(output_pdb_path + object1 + "_align_model.pdb"))
				cmd.do("load %s" %(output_pdb_path + object1 + "_chain_A_align_model.pdb"))
				cmd.do("set_name %s, align_model" %(object1 + "_align_model"))
				cmd.do("set_name %s, chain_A" %(object1 + "_chain_A_align_model"))
				for i in range(0,nsubs):
					chain_ord = 65 + i
					chain_id = chr(chain_ord)
					if i != 0:
						cmd.do("create chain_%s, chain_A" %(chain_id))
						cmd.do("alter %s, chain='%s'" %("chain_"+chain_id, chain_id)) 
						cmd.do("super %s, (align_model and chain %s)" %("chain_"+chain_id, chain_id))
					merge_list.append("chain_" + chain_id )

				merge_string = " + ".join(merge_list)
				cmd.do("create merged, (%s)" %(merge_string))

				# Run Will's pymol scripts to align the assembly along the z-axis and save the new transformed pdbs. Run twice because seems to be starting point dependent and will sometimes misalign the first time.

				merged = "merged"
				if nsubs == 2:
					cmd.do("tmp = aligncx(merged,2)")
					cmd.do("tmp = aligncx(merged,2)")
				elif nsubs == 3:
					cmd.do("tmp = aligncx(merged,3)")
					cmd.do("tmp = aligncx(merged,3)")
				elif nsubs == 4:
					cmd.do("tmp = aligncx(merged,4)")
					cmd.do("tmp = aligncx(merged,4)")
				elif nsubs == 5:
					cmd.do("tmp = aligncx(merged,5)")
					cmd.do("tmp = aligncx(merged,5)")
				else:
					cmd.do("tmp = aligncx(merged,6)")
					cmd.do("tmp = aligncx(merged,6)")
				results_list.append(tmp)
				cmd.reinitialize()

			best_sym = None
			# Pick the result that gave acceptable errors with the highest number of chains.  A C4 assembly could for instance happen to also align well as C2 if the right 2 chains were analyzed, but if it aligns well as C4,
			# then this higher symmetry would typically be the correct choice. Caution: It could be a D2, however, which is not distinguished well from C4 by Will's alignment method.
			# Could we look for orthoganol axes?
			for index in range(0,len(results_list)):
				if ((float(results_list[index][2]) < float(options.max_distance_error)) and (float(results_list[index][3]) < float(options.max_angle_error))):  
					best_sym = index+2
			if best_sym != None:

				if options.add_to_database:
					entry_dict['sym'] = "C" + str(best_sym)

				merge_list = []
				cmd.do("load %s" %(output_pdb_path + object1 + "_align_model.pdb"))
				cmd.do("load %s" %(output_pdb_path + object1 + "_chain_A_align_model.pdb"))
				cmd.do("set_name %s, align_model" %(object1 + "_align_model"))
				cmd.do("set_name %s, chain_A" %(object1 + "_chain_A_align_model"))
				super_string = ""
				for i in range(0,best_sym):
					chain_ord = 65 + i
					chain_id = chr(chain_ord)
					if i != 0:
						cmd.do("create chain_%s, chain_A" %(chain_id))
						cmd.do("alter %s, chain='%s'" %("chain_"+chain_id, chain_id)) 
						super_string+="\nPROCESSING   Superposition_Report: chain %s" %(chain_id) + " RMSD = " + str(cmd.super("chain_%s" %(chain_id), "(align_model and chain %s)" %(chain_id), object=algn)[0])
					merge_list.append("chain_" + chain_id )

				merge_string = " + ".join(merge_list)
				cmd.do("create merged, (%s)" %(merge_string))

				merged = "merged"
				if best_sym == 2:
					cmd.do("tmp = aligncx(merged,2)")
					cmd.do("tmp = aligncx(merged,2)")
				elif best_sym == 3:
					cmd.do("tmp = aligncx(merged,3)")
					cmd.do("tmp = aligncx(merged,3)")
				elif best_sym == 4:
					cmd.do("tmp = aligncx(merged,4)")
					cmd.do("tmp = aligncx(merged,4)")
				elif best_sym == 5:
					cmd.do("tmp = aligncx(merged,5)")
					cmd.do("tmp = aligncx(merged,5)")
				else:
					cmd.do("tmp = aligncx(merged,6)")
					cmd.do("tmp = aligncx(merged,6)")
				m = cmd.get_view(0)
				ttt = [m[0], m[1], m[2], 0.0, m[3], m[4], m[5], 0.0, m[6], m[7], m[8], 0.0, 0.0,   0.0,  0.0, 1.0]
				cmd.transform_object("merged",ttt)
				cmd.do("create aligned_chainA, (merged and chain A)")
				cmd.do("save %s, merged" %(output_pdb_path + "C" + str(best_sym) + "_" + output_pdbname[:-4] +"_full_tmp.pdb"))
				cmd.do("save %s, aligned_chainA" %(output_pdb_path + "C" + str(best_sym) + "_" + output_pdbname[:-4] + "_tmp.pdb"))

				if ((best_sym == len(chains)) and ((options.expected_symmetry == None) or ("C" + str(best_sym) == options.expected_symmetry))) :
					log_string3 = "PROCESSING   Alignment_Summary: Chain A was copied, aligned to each of the other " + str(best_sym-1) + " chains, merged into one object and aligned with the z axis. \nPROCESSING   Alignment_Summary: The new aligned chain A was saved as " + output_pdb_path + "C" + str(best_sym) + "_" + output_pdbname +"\nPROCESSING   Construct_Chain_Length: " + str(construct_length) + "\nPROCESSING   Chain_Length: " + str(rescount) 
				elif ((best_sym != len(chains)) and (options.expected_symmetry != None) and ("C" + str(best_sym) != options.expected_symmetry)):
					log_string3 = "PROCESSING   Alignment_Summary: Chain A was copied, aligned to each of the next " + str(best_sym-1) + " chains, merged into one object and aligned with the z axis. \nPROCESSING   Alignment_Summary: The new aligned chain a was saved as " + output_pdb_path + "c" + str(best_sym) + "_" + output_pdbname + "\nPROCESSING   Construct_Chain_Length:" + str(construct_length) + "\nPROCESSING   Chain_Length: " + str(rescount) + "\nPROCESSING   Expected_Symmetry_Warning: The determined symmetry, C" + str(best_sym) +" , does not match with the expected symmetry, " + options.expected_symmetry + ".\nPROCESSING   Extra_Chains_Warning: More chains in " + output_pdb_path + object1 + "_ref.pdb than in the aligned complex. The actual assembly may be higher order."
					ex_chains = True
					#try:
					#	entry_dict['notes']+=" Expected_Symmetry_Warning: The determined symmetry, C" + str(best_sym) +" , does not match with the expected symmetry, " + options.expected_symmetry + " Extra_Chains_Warning: More chains in " + output_pdb_path + object1 + "_ref.pdb than in the aligned complex. The actual assembly may be higher order."
					#except KeyError:
					#	entry_dict['notes']="Expected_Symmetry_Warning: The determined symmetry, C" + str(best_sym) +" , does not match with the expected symmetry, " + options.expected_symmetry + " Extra_Chains_Warning: More chains in " + output_pdb_path + object1 + "_ref.pdb than in the aligned complex. The actual assembly may be higher order." 
				elif ((best_sym == len(chains)) and (options.expected_symmetry != None) and ("C" + str(best_sym) != options.expected_symmetry)):
					log_string3 = "PROCESSING   Alignment_Summary: Chain A was copied, aligned to each of the next " + str(best_sym-1) + " chains, merged into one object and aligned with the z axis.\nPROCESSING   Alignment_Summary: The new aligned chain a was saved as " + output_pdb_path + "c" + str(best_sym) + "_" + output_pdbname + "\nPROCESSING   Construct_Chain_Length:" + str(construct_length)  + "\nPROCESSING   Chain_Length: " + str(rescount) + "\nPROCESSING   Expected_Symmetry_Warning: The determined symmetry, C" + str(best_sym) +" , does not match with the expected symmetry, " + options.expected_symmetry + "."
					#try:
					#	entry_dict['notes']+=" Expected_Symmetry_Warning: The determined symmetry, C" + str(best_sym) +" , does not match with the expected symmetry, " + options.expected_symmetry
					#except KeyError:
					#	entry_dict['notes']="Expected_Symmetry_Warning: The determined symmetry, C" + str(best_sym) +" , does not match with the expected symmetry, " + options.expected_symmetry

				else:
					log_string3 = "PROCESSING   Alignment_Summary: Chain A was copied, aligned to each of the next " + str(best_sym-1) + " chains, merged into one object and aligned with the z axis.\nPROCESSING   Alignment_Summary: The new aligned chain a was saved as " + output_pdb_path + "c" + str(best_sym) + "_" + output_pdbname + "\nPROCESSING   Construct_Chain_Length:" + str(construct_length)  + "\nPROCESSING   Chain_Length: " + str(rescount) + "\nPROCESSING   Extra_Chains_Warning: more chains in " + output_pdb_path + object1 + "_ref.pdb than in the aligned complex. The actual assembly may be higher order."
					ex_chains = True
					#try:
					#	entry_dict['notes']+=" Extra_Chains_Warning: More chains in " + output_pdb_path + object1 + "_ref.pdb than in the aligned complex. The actual assembly may be higher order."
					#except KeyError:
					#	entry_dict['notes']="Extra_Chains_Warning: More chains in " + output_pdb_path + object1 + "_ref.pdb than in the aligned complex. The actual assembly may be higher order."

				alignment_report_string = "\nPROCESSING   Axis_Alignment_Report: Axis = " + str(tmp[0])
				alignment_report_string+= "\nPROCESSING   Axis_Alignment_Report: Center = " + str(tmp[1])
				alignment_report_string+= "\nPROCESSING   Axis_Alignment_Report: Distance_Error = " + str(tmp[2])
				alignment_report_string+= "\nPROCESSING   Axis_Alignment_Report: Angle_Error = " + str(tmp[3])

				if options.add_to_database:
					entry_dict['dist_err'] = "%.3f" % round(tmp[2],3)
					entry_dict['ang_err'] = "%.3f" % round(tmp[3],3) 
					entry_dict['ex_chains'] = str(ex_chains)
					entry_dict['ex_models'] = str(ex_models)

				log_string3 = "\n"+log_string3 + super_string + alignment_report_string
				if not complete_consensus:
					log_string3+="\nPROCESSING   Complete_Consensus: False"+resi_string
				else:
					log_string3+="\nPROCESSING   Complete_Consensus: True"

				# Add Header Info, comment_string, and new log_string3 to the pdb. Remove the temporary files that were missing the header info.

				comment_string+=log_string3

				full_tmp_file = open(output_pdb_path + "C" + str(best_sym) + "_" + output_pdbname[:-4] +"_full_tmp.pdb","r")
				full_tmp_string = full_tmp_file.read()
				full_tmp_file.close()

				A_tmp_file = open(output_pdb_path + "C" + str(best_sym) + "_" + output_pdbname[:-4] + "_tmp.pdb","r")
				A_tmp_string = A_tmp_file.read()
				A_tmp_file.close()

				outfile_full_aligned_pdb = open(output_pdb_path + "C" + str(best_sym) + "_" + output_pdbname[:-4] + "_full.pdb", "w")
				outfile_full_aligned_pdb.write(comment_string+"\n"+header_string+"\n"+full_tmp_string)
				outfile_full_aligned_pdb.close()

				outfile_A_aligned_pdb = open(output_pdb_path + "C" + str(best_sym) + "_" + output_pdbname, "w")
				outfile_A_aligned_pdb.write(comment_string+"\n"+header_string+"\n"+A_tmp_string)
				outfile_A_aligned_pdb.close()

				subprocess.check_output("rm " + output_pdb_path + "C" + str(best_sym) + "_" + output_pdbname[:-4] + "_full_tmp.pdb",stderr=subprocess.STDOUT,shell=True)
				subprocess.check_output("rm " + output_pdb_path + "C" + str(best_sym) + "_" + output_pdbname[:-4] + "_tmp.pdb",stderr=subprocess.STDOUT,shell=True)
				if not options.view_in_pymol:
					subprocess.check_output("rm " + output_pdb_path + object1 + "_align_model.pdb",stderr=subprocess.STDOUT,shell=True)
					subprocess.check_output("rm " + output_pdb_path + object1 + "_chain_A_align_model.pdb",stderr=subprocess.STDOUT,shell=True)
					subprocess.check_output("rm " + output_pdb_path + object1 + ".pdb",stderr=subprocess.STDOUT,shell=True)
					print(input_pdbname + " reports success: Aligned to Z axis with C" + str(best_sym) + " symmetry.")
					if options.add_to_database:
						new_rows.append(add_entry(database_headers=database_headers,refids=refids,my_dict=entry_dict))
				elif options.view_in_pymol and (len(pdb_list) == 1):
					subprocess.check_output("rm " + output_pdb_path + object1 + "_align_model.pdb",stderr=subprocess.STDOUT,shell=True)
					subprocess.check_output("rm " + output_pdb_path + object1 + "_chain_A_align_model.pdb",stderr=subprocess.STDOUT,shell=True)
					subprocess.check_output("rm " + output_pdb_path + object1 + ".pdb",stderr=subprocess.STDOUT,shell=True)
					print(input_pdbname + " reports success: Aligned to Z axis with C" + str(best_sym) + " symmetry.")
					if options.add_to_database:
						new_rows.append(add_entry(database_headers=database_headers,refids=refids,my_dict=entry_dict))
					#cmd.do("run %s/axes.py" %(options.pymol_scripts_dir))
			else:
				subprocess.check_output("rm " + output_pdb_path + object1 + "_align_model.pdb",stderr=subprocess.STDOUT,shell=True)
				subprocess.check_output("rm " + output_pdb_path + object1 + "_chain_A_align_model.pdb",stderr=subprocess.STDOUT,shell=True)
				subprocess.check_output("rm " + output_pdb_path + object1 + ".pdb",stderr=subprocess.STDOUT,shell=True)
				print(input_pdbname + " reports failure: Not found to possess C2, C3, C4, C5, or C6 symmetry.")

		biounit_number+=1
		if single_biounit_mode:
			another_biounit = False
		elif not os.path.exists(input_pdb_path+pdb.strip()+".pdb"+str(biounit_number)+".gz"):
			another_biounit = False
		else:
			input_pdbname = pdb.strip() + ".pdb" + str(biounit_number) + ".gz"
			pdbid_bunit = pdb.strip() + "_" + str(biounit_number)
			output_pdbname = pdbid_bunit + ".pdb"
			if os.path.exists(output_pdb_path+output_pdbname):
				subprocess.check_output("rm " + output_pdb_path + output_pdbname,stderr=subprocess.STDOUT,shell=True) 
			subprocess.check_output("cp " + input_pdb_path + input_pdbname + " " + output_pdb_path + output_pdbname + ".gz",stderr=subprocess.STDOUT,shell=True) 
			subprocess.check_output("gunzip " + output_pdb_path + output_pdbname + ".gz",stderr=subprocess.STDOUT,shell=True) 
			log_string1 = "PROCESSING   InputPDB: " + input_pdb_path + input_pdbname
			subprocess.check_output("cp " + output_pdb_path + output_pdbname + " " + output_pdb_path + output_pdbname.replace(".pdb","") + "_ref.pdb",stderr=subprocess.STDOUT,shell=True) 
		cmd.reinitialize()

if options.add_to_database:
	data_string = ""
	new_array = [database_headers]
	formatted_string_list = []
	try:
		for row in database_array:
			new_array.append(row.tolist())
	except TypeError:
		new_array.append(database_array.tolist())
	for row in new_rows:
		new_array.append(row)
	for row in new_array:
		formatted_data_list = []
		for i in range(0,len(database_headers)):
			formatted_data_list.append(str(row[i]).rjust(len(database_headers[i].strip())))
		row_string = "\t".join(formatted_data_list)
		formatted_string_list.append(row_string)

	database_file = open(os.getenv("HOME")+"/databases/symmetric_scaffolds/symmetric_scaffolds.dat","w")
	database_file.write("\n".join(formatted_string_list))
	database_file.close()

cmd.quit() 
