### make aligned input monomer into a D3
maked3(sel="*model", name="model_D3")

### make 4 copies of D3 to align in crystal
cmd.copy("d3_1", "model_D3")
cmd.copy("d3_2", "model_D3")
cmd.copy("d3_3", "model_D3")
cmd.copy("d3_4", "model_D3")

### make 4 copies of designed cage output
cmd.set_name("combined_1ehw_chain_A_align_model_C2_jf21_1_O32F_7_1_0086", "cage")
cmd.copy("cage_1", "cage")
cmd.copy("cage_2", "cage")
cmd.copy("cage_3", "cage")
cmd.copy("cage_4", "cage")

### align d3_1 and d3_2 to cage_1
super d3_1 and chain A, cage_1 and chain A
super d3_2 and chain A, cage_1 and chain O

###
super cage_2 and chain A, d3_1 and chain E
super cage_3 and chain A, d3_2 and chain E

###
super d3_3 and chain A, cage_3 and chain U
super d3_4 and chain A, cage_2 and chain 0

###
super cage_4 and chain A, d3_3 and chain E