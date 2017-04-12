remove (hydro)
bg_color white
hide everything
set surface_quality, 1
alter all, b=50
alter all, q=1
set gaussian_resolution, 12

map_new map0, gaussian, 1, chain 0 and name ca, 12
isosurface surf0, map0
color blue, surf0
map_new map1, gaussian, 1, chain 1 and name ca, 12
isosurface surf1, map1
color green, surf1
map_new map2, gaussian, 1, chain 2 and name ca, 12
isosurface surf2, map2
color blue, surf2
map_new map3, gaussian, 1, chain 3 and name ca, 12
isosurface surf3, map3
color green, surf3
map_new map4, gaussian, 1, chain 4 and name ca, 12
isosurface surf4, map4
color blue, surf4
map_new map5, gaussian, 1, chain 5 and name ca, 12
isosurface surf5, map5
color green, surf5
map_new map6, gaussian, 1, chain 6 and name ca, 12
isosurface surf6, map6
color blue, surf6
map_new map7, gaussian, 1, chain 7 and name ca, 12
isosurface surf7, map7
color green, surf7
map_new map8, gaussian, 1, chain 8 and name ca, 12
isosurface surf8, map8
color blue, surf8
map_new map9, gaussian, 1, chain 9 and name ca, 12
isosurface surf9, map9
color green, surf9
map_new mapA, gaussian, 1, chain A and name ca, 12
isosurface surfA, mapA
color blue, surfA
map_new mapB, gaussian, 1, chain B and name ca, 12
isosurface surfB, mapB
color green, surfB
map_new mapC, gaussian, 1, chain C and name ca, 12
isosurface surfC, mapC
color blue, surfC
map_new mapD, gaussian, 1, chain D and name ca, 12
isosurface surfD, mapD
color green, surfD
map_new mapE, gaussian, 1, chain E and name ca, 12
isosurface surfE, mapE
color blue, surfE
map_new mapF, gaussian, 1, chain F and name ca, 12
isosurface surfF, mapF
color green, surfF
map_new mapG, gaussian, 1, chain G and name ca, 12
isosurface surfG, mapG
color green, surfG
map_new mapH, gaussian, 1, chain H and name ca, 12
isosurface surfH, mapH
color blue, surfH
map_new mapI, gaussian, 1, chain I and name ca, 12
isosurface surfI, mapI
color green, surfI
map_new mapJ, gaussian, 1, chain J and name ca, 12
isosurface surfJ, mapJ
color blue, surfJ
map_new mapK, gaussian, 1, chain K and name ca, 12
isosurface surfK, mapK
color green, surfK
map_new mapL, gaussian, 1, chain L and name ca, 12
isosurface surfL, mapL
color blue, surfL
map_new mapM, gaussian, 1, chain M and name ca, 12
isosurface surfM, mapM
color green, surfM
map_new mapN, gaussian, 1, chain N and name ca, 12
isosurface surfN, mapN
color blue, surfN
map_new mapO, gaussian, 1, chain O and name ca, 12
isosurface surfO, mapO
color green, surfO
map_new mapP, gaussian, 1, chain P and name ca, 12
isosurface surfP, mapP
color blue, surfP
map_new mapQ, gaussian, 1, chain Q and name ca, 12
isosurface surfQ, mapQ
color green, surfQ
map_new mapR, gaussian, 1, chain R and name ca, 12
isosurface surfR, mapR
color blue, surfR
map_new mapS, gaussian, 1, chain S and name ca, 12
isosurface surfS, mapS
color green, surfS
map_new mapT, gaussian, 1, chain T and name ca, 12
isosurface surfT, mapT
color blue, surfT
map_new mapU, gaussian, 1, chain U and name ca, 12
isosurface surfU, mapU
color green, surfU
map_new mapV, gaussian, 1, chain V and name ca, 12
isosurface surfV, mapV
color blue, surfV
map_new mapW, gaussian, 1, chain W and name ca, 12
isosurface surfW, mapW
color green, surfW
map_new mapX, gaussian, 1, chain X and name ca, 12
isosurface surfX, mapX
color blue, surfX
map_new mapY, gaussian, 1, chain Y and name ca, 12
isosurface surfY, mapY
color green, surfY
map_new mapZ, gaussian, 1, chain Z and name ca, 12
isosurface surfZ, mapZ
color blue, surfZ
map_new mapa, gaussian, 1, chain a and name ca, 12
isosurface surfa, mapa
color green, surfa
map_new mapb, gaussian, 1, chain b and name ca, 12
isosurface surfb, mapb
color blue, surfb
map_new mapc, gaussian, 1, chain c and name ca, 12
isosurface surfc, mapc
color green, surfc
map_new mapd, gaussian, 1, chain d and name ca, 12
isosurface surfd, mapd
color blue, surfd
map_new mape, gaussian, 1, chain e and name ca, 12
isosurface surfe, mape
color green, surfe
map_new mapf, gaussian, 1, chain f and name ca, 12
isosurface surff, mapf
color blue, surff
map_new mapg, gaussian, 1, chain g and name ca, 12
isosurface surfg, mapg
color green, surfg
map_new maph, gaussian, 1, chain h and name ca, 12
isosurface surfh, maph
color blue, surfh
map_new mapi, gaussian, 1, chain i and name ca, 12
isosurface surfi, mapi
color green, surfi
map_new mapj, gaussian, 1, chain j and name ca, 12
isosurface surfj, mapj
color blue, surfj
map_new mapk, gaussian, 1, chain k and name ca, 12
isosurface surfk, mapk
color green, surfk
map_new mapl, gaussian, 1, chain l and name ca, 12
isosurface surfl, mapl
color blue, surfl


set antialias, 2
set ray_trace_gain, 0.4
set ray_shadows, 0
set specular, 0
show surface, surf*
