########## USEFUL COMMANDS ##########
#alignaxis("vis", Vec(0,0,1), com("vis))
#print com("vis")
#translate [0,0,Z], selection=vis, camera=0
#trans("vis",Vec(0,0,Z))
#aligncx(sele="rop4",nfold="2",chains="AB")
#set cartoon_cylindrical_helices, 1
#alter 1cll, resi=str(int(resi)-4)
#
#flip something 180
#rotate x, 180, camera=0, object=OBJECT, origin=[0,0,0]
#####################################

util.cbc
alias xxx, util.cbc; util.cnc
util.cnc
alias zzz, util.cnc

set all_states, off 
set seq_view, 0


########## Set up aliases ##########
alias reloadrc, @~/.pymolrc
alias ard, show sphere, sele + (sele around 6) and visible 
alias hs, hide sphere; hide surface 
alias sc, show cartoon
alias hc, hide cartoon 
alias sortobj, order *, yes
alias himain, hide line, backbone; show line, name CA
alias show main, show line, backbone + name CB; hide line, hydro
alias fasta, print cmd.get_fastastr("sele")
alias listed, list=[]; iterate ("sele" and name ca),list.append((resi,resn)); print(list)
alias findterm, show spheres, name 1H + name OXT
alias viewx, center sele; orient sele; zoom; center visible

alias selabc, sele chainA, chain A and vis; sele chain B, chain B and vis; sele chain C, chain C and vis

alias cyl, hide all; show car; set cartoon_transparency, 0;  set cartoon_cylindrical_helices, 1

alias fastaa, sele chain A; fasta; sele chain B; fasta

alias render, set specular, 0; set ray_shadow, off; set valence, off; set antialias, 2; set ray_trace_mode, 1; set ray_trace_disco_factor, 1; set ray_trace_gain, 0.1; set power, 0.2; set ambient, 0.4; set cartoon_transparency, 0.2

alias figure, hide lines; hide spheres; hide labels; set cartoon_ladder_mode, 1; set dash_color, yellow; set cartoon_ring_mode, 3; set ray_trace_mode, 1; set ray_trace_color, black; set cartoon_fancy_helices=1; set cartoon_cylindrical_helices, 2; set cartoon_loop_radius=0.1; set cartoon_smooth_loops = 0; set light, (0,-5,-20); set direct = .5; set ray_shadows=0; set ray_opaque_background, off; set ray_trace_mode,  1; set ambient, 1; set reflect, 0; set two_sided_lighting, on; cartoon oval; set cartoon_oval_width, 0.3; set cartoon_oval_length, 1.1; ray

alias 2figure, hide lines; hide labels; set transparency, 0; set sphere_transparency, 1; set cartoon_transparency, 0; cartoon loop; set cartoon_fancy_helices, 0; set cartoon_loop_radius, 2; set cartoon_flat_sheets, 0; set cartoon_smooth_loops, 0; ray

alias fig, hide lines; hide spheres; hide labels; set dash_color, yellow; set ray_trace_mode, 1; set ray_trace_color, black; set light, (0,-5,-20); set direct = .5; set ray_shadows=0; set ray_opaque_background, off; set ray_trace_mode,  1; set ambient, 1; set reflect, 0; set two_sided_lighting, on; ray

############ Color by feature ############

#color by polar/non-polar
alias clrd, color grey60; color green, resn ASP+GLU+ASN+GLN+HIS+LYS+ARG; color palegreen, resn TYR+TRP; color marine, resn SER+THR; color pink, resn GLY, color teal, resn PRO; color white, resn ALA; util.cnc

#color by ss
alias clrss, color firebrick, ss H; color marine, ss S; color grey60, ss L
#hide not_designs
alias drk, select not_design_resis, ! design_resis; color grey20, not_design_resis; util.cnc
alias drkk, select not_design_resis, ! PDBinfo-LABEL; color grey20, not_design_resis; util.cnc

#find polar contacts
alias fspc, sele vis and sc.; delete polar_contacts; distance polar_contacts, sele, visible, mode=2; show line, don. extend 1 and h.; hide label

#hide sidechains
alias dockmode, set all_states, on; hide line; show line, motif*; show ribbon; showmain; util.cbc

alias cm, run ~/pymol_apps/cage_interface_viewer.pml

alias hydrogens, hide (h.)

#set transparency_mode to 1 for ray trace
set transparency_mode, 2

#set surface
	set transparency, 0.6

show lines
	hide (hydro)

show ribbon
	set ribbon_width, 9
	#set ribbon_smooth, 1

#show cartoon
	set cartoon_transparency, 0.0
	set cartoon_flat_sheets, 0
	set cartoon_smooth_loops, 0
	set cartoon_fancy_helices, 0

#show sphere
	set sphere_transparency, 0.5

#raytrace mode
	set ray_trace_mode, 1
	set ray_shadows, 0
	set ray_shadow, 0
	set ray_trace_gain, 0.01
	set antialias, 2
	set ray_opaque_background, off

#glow on atoms 
	set specular, 0

set orthoscopic,  on
	set line_smooth, 1

#set surface cavity mode
	set surface_cavity_mode, 2

#shows NT
	show spheres, name 1H
#shows CT
	show spheres, name OXT

######## Run pymol scripts #######
run /home/erinyang/pymol_apps/pymol/obj_arrows.py
run /home/erinyang/pymol_apps/pymol/intsele.py
run /home/erinyang/pymol_apps/pymol/center_of_mass.py
run /home/erinyang/pymol_apps/pymol/drawline.py
run /home/erinyang/pymol_apps/pymol/cgocircle.py
run /home/erinyang/pymol_apps/pymol/marshmallow.py
run /home/erinyang/pymol_apps/pymol/save_image.py
run /home/erinyang/pymol_apps/pymol/mesh_view.py
run /home/erinyang/pymol_apps/pymol/cage_format.py
run /home/erinyang/pymol_apps/pymol/xyzMath.py
run /home/erinyang/pymol_apps/pymol/xyzGeom.py
run /home/erinyang/pymol_apps/pymol/pymol_util.py
run /home/erinyang/pymol_apps/pymol/sym_util.py
run /home/erinyang/pymol_apps/pymol/symgen.py
run /home/erinyang/pymol_apps/pymol/axes.py
run /home/erinyang/pymol_apps/pymol/util.py
run /home/erinyang/pymol_apps/pymol/PyMOL-RosettaServer.py
run /home/erinyang/pymol_apps/pymol/highlight_dif.py
run /home/erinyang/pymol_apps/pymol/symgen_goldberg_polyhedra.py
run /home/erinyang/pymol_apps/pymol/align_util.py 

run /home/erinyang/pymol_apps/annotate/annotate.py
cmd.set_key('F10', annotate_bad)
cmd.set_key('F12', annotate_good)

run /home/erinyang/pymol_apps/pymol/Cycler.py
run /home/erinyang/pymol_apps/pymol/colors.py
run /home/erinyang/pymol_apps/pymol/mass_alignment.py
run /home/erinyang/pymol_apps/pymol/GenUtils.py
run /home/erinyang/pymol_apps/pymol/mass_good_rainbow.py

hide everything

zoom
sortobj
#clrd
show cartoon
show dash
show labels

#set mouse to 3-button
mouse three_button_viewing
