### make aligned input monomer into a D3
maked3(sel="3sba_chain_A_align_model", name="model_D3")

### make 4 copies of D3 to align in crystal
cmd.copy("D3_1", "model_D3")
cmd.copy("D3_2", "model_D3")
cmd.copy("D3_3", "model_D3")
cmd.copy("D3_4", "model_D3")


### make 4 copies of designed cage output
cmd.set_name("combined_3sba_chain_A_align_model_C2_3hm4_1_O32F_10_1_0001", "cage")
cmd.copy("cage_1", "cage")
cmd.copy("cage_2", "cage")
cmd.copy("cage_3", "cage")
cmd.copy("cage_4", "cage")

### align D3_1 and D3_2 to cage_1
super D3_1 and chain A, cage_1 and chain A
super D3_2 and chain A, cage_1 and chain 6

###
super cage_2 and chain A, D3_1 and chain D
super cage_3 and chain A, D3_2 and chain D

###
super D3_4 and chain A, cage_2 and chain Q
super D3_3 and chain A, cage_3 and chain I

###
super cage_4 and chain A, D3_3 and chain E
