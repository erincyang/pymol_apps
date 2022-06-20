hide everything
show ribbon
show lines
hide lines, name C+O+N+*H*
show lines, (resn PRO and name N)+(name NH1+NH2)+(name OH)+(name CH2)
show cartoon
set cartoon_transparency, 0.75
util.cbc(quiet=1)
color nitrogen, name N*
color oxygen, name O*
color sulfur, name S*
color hydrogen, hydrogens
set ray_trace_mode, 1
set antialias, 2
set ray_shadows, 0
set specular, 0
set use_shaders, 1
