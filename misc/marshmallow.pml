remove (hydro)
bg_color white
hide everything
set surface_quality, 1
alter all, b=50
alter all, q=1
set gaussian_resolution, 12
map_new mapA, gaussian, 1, chain A and name ca, 12
map_new mapC, gaussian, 1, chain C and name ca, 12
map_new mapE, gaussian, 1, chain E and name ca, 12
map_new mapM, gaussian, 1, chain M and name ca, 12
map_new mapQ, gaussian, 1, chain Q and name ca, 12
map_new mapO, gaussian, 1, chain O and name ca, 12
map_new mapI, gaussian, 1, chain I and name ca, 12
map_new mapG, gaussian, 1, chain G and name ca, 12
map_new mapK, gaussian, 1, chain K and name ca, 12
map_new mapU, gaussian, 1, chain U and name ca, 12
map_new mapW, gaussian, 1, chain W and name ca, 12
map_new mapS, gaussian, 1, chain S and name ca, 12
map_new mapL, gaussian, 1, chain L and name ca, 12
map_new mapV, gaussian, 1, chain V and name ca, 12
map_new mapT, gaussian, 1, chain T and name ca, 12
map_new mapF, gaussian, 1, chain F and name ca, 12
map_new mapD, gaussian, 1, chain D and name ca, 12
map_new mapN, gaussian, 1, chain N and name ca, 12
map_new mapJ, gaussian, 1, chain J and name ca, 12
map_new mapP, gaussian, 1, chain P and name ca, 12
map_new mapX, gaussian, 1, chain X and name ca, 12
map_new mapR, gaussian, 1, chain R and name ca, 12
map_new mapH, gaussian, 1, chain H and name ca, 12
map_new mapB, gaussian, 1, chain B and name ca, 12
isosurface surfA, mapA
isosurface surfC, mapC
isosurface surfE, mapE
isosurface surfM, mapM
isosurface surfQ, mapQ
isosurface surfO, mapO
isosurface surfI, mapI
isosurface surfG, mapG
isosurface surfK, mapK
isosurface surfU, mapU
isosurface surfW, mapW
isosurface surfS, mapS
isosurface surfL, mapL
isosurface surfV, mapV
isosurface surfT, mapT
isosurface surfF, mapF
isosurface surfD, mapD
isosurface surfN, mapN
isosurface surfJ, mapJ
isosurface surfP, mapP
isosurface surfX, mapX
isosurface surfR, mapR
isosurface surfH, mapH
isosurface surfB, mapB
color cyan, surfA
color cyan, surfC
color cyan, surfE
color cyan, surfM
color cyan, surfQ
color cyan, surfO
color cyan, surfI
color cyan, surfG
color cyan, surfK
color cyan, surfU
color cyan, surfW
color cyan, surfS
color purple, surfL
color purple, surfV
color purple, surfT
color purple, surfF
color purple, surfD
color purple, surfN
color purple, surfJ
color purple, surfP
color purple, surfX
color purple, surfR
color purple, surfH
color purple, surfB
set antialias, 2
set ray_trace_gain, 0.4
set ray_shadows, 0
set specular, 0
show surface, surf*
