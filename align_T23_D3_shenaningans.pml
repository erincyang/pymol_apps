### make aligned input monomer into a D3
#maked3(sel="D3_1552_1", name="model_D3")

### make 6 copies of D3 to align in crystal (hard coded)
cmd.copy("model_D3","UWN_455_C3_1na0-1_m0_chA_C2_dummy_asu_D32F_27_1_0025")

cmd.copy("D3_1", "model_D3")
cmd.copy("D3_2", "model_D3")
cmd.copy("D3_3", "model_D3")
cmd.copy("D3_4", "model_D3")
cmd.copy("D3_5", "model_D3")
cmd.copy("D3_6", "model_D3")


### make 6 copies of hard coded designed cage output
cmd.set_name("T23d_top0_0__compA_C2_3L6H-4_1__compB_455_asym_rot15_0002", "cage")

cmd.copy("cage_1", "cage")
cmd.copy("cage_2", "cage")
cmd.copy("cage_3", "cage")
cmd.copy("cage_4", "cage")
cmd.copy("cage_5", "cage")
cmd.copy("cage_6", "cage")

### align D3_1 and D3_2 to cage_1
super D3_1 and chain A, cage_1 and chain B
super D3_2 and chain A, cage_1 and chain X

###
super cage_2 and chain B, D3_1 and chain E
super cage_3 and chain B, D3_2 and chain E

###
super D3_4 and chain A, cage_2 and chain X
super D3_3 and chain A, cage_3 and chain N

###
super cage_4 and chain B, D3_3 and chain E
super cage_5 and chain B, D3_4 and chain E

###
super D3_5 and chain A, cage_4 and chain R
super D3_6 and chain A, cage_5 and chain N

###
super cage_6 and chain B, D3_5 and chain E

