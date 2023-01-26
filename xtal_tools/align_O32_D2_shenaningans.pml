### make aligned input monomer into a D2
maked2(sel="UWN-553_C2_ank1-G3-r5_full_m0_chA_C2_dummy_asu_D22F_66_1_0037", name="model_D2")

### make 4 copies of D2 to align in crystal
cmd.copy("D2_1", "model_D2")
cmd.copy("D2_2", "model_D2")
cmd.copy("D2_3", "model_D2")
#cmd.copy("D2_4", "model_D2")

### make 4 copies of designed cage output
cmd.set_name("rpxdock_top0_0__compA_1BH_69_chA__compB_UWN-553_rot45_asu", "cage")
cmd.copy("cage_1", "cage")
cmd.copy("cage_2", "cage")
cmd.copy("cage_3", "cage")
#cmd.copy("cage_4", "cage")

### align D2_1 and D2_2 to cage_1
super D2_1 and chain A, cage_1 and chain B
super D2_2 and chain A, cage_1 and chain N

###
super cage_2 and chain B, D2_1 and chain C
super cage_3 and chain B, D2_2 and chain C

###
#super D2_3 and chain A, cage_3 and chain B
super D2_3 and chain A, cage_3 and chain J

###
#super cage_4 and chain A, D2_3 and chain E
