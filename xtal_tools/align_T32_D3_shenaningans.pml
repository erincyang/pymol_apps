### make aligned input monomer into a D3
#maked3(sel="D3_1552_1", name="model_D3")

### make 6 copies of D3 to align in crystal (hard coded)
cmd.copy("model_D3","D3_1wa3-33_rot15")

cmd.copy("D3_1", "model_D3")
cmd.copy("D3_2", "model_D3")
cmd.copy("D3_3", "model_D3")
cmd.copy("D3_4", "model_D3")
cmd.copy("D3_5", "model_D3")
cmd.copy("D3_6", "model_D3")


### make 6 copies of hard coded designed cage output
cmd.set_name("TET", "cage")

cmd.copy("cage_1", "cage")
cmd.copy("cage_2", "cage")
cmd.copy("cage_3", "cage")
cmd.copy("cage_4", "cage")
cmd.copy("cage_5", "cage")
cmd.copy("cage_6", "cage")

### align D3_1 and D3_2 to cage_1
super D3_1 and chain A, cage_1 and chain P
super D3_2 and chain A, cage_1 and chain L

###
super cage_2 and chain B, D3_1 and chain D
super cage_3 and chain B, D3_2 and chain D

###
super D3_4 and chain A, cage_2 and chain J
super D3_3 and chain A, cage_3 and chain X

###
super cage_4 and chain B, D3_3 and chain D
super cage_5 and chain B, D3_4 and chain D

###
super D3_5 and chain A, cage_4 and chain X
super D3_6 and chain A, cage_5 and chain X

###
super cage_6 and chain B, D3_5 and chain D

