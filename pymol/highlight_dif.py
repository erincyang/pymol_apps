from pymol import cmd

def comp_obj_seq():
    #assume first model is ref model
    obj_list = cmd.get_object_list()

    ref = obj_list[0]
    sele_string = ref + " and chain A+B"
    ref_fasta = cmd.get_fastastr(sele_string).split("\n", 1)[1].replace("\n","")
    for obj in obj_list[1:]:
        sele_string = obj + " and chain A+B"
        obj_fasta = cmd.get_fastastr(sele_string).split("\n", 1)[1].replace("\n","")
        if len(ref_fasta) != len(obj_fasta):
            print( "sequence lengths are not the same!")
        resi_list = [i for i in xrange(len(ref_fasta)) if ref_fasta[i] != obj_fasta[i]]
    for resi in resi_list:
        cmd.show("sticks", "resi " + str(resi+1))
    return None

cmd.extend('comp_obj_seq', comp_obj_seq)

