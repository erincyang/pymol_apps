remove (hydro)
#bg_color white
hide everything
set surface_quality, 1
alter all, b=50
alter all, q=1
set gaussian_resolution, 12

map_new mapB, gaussian, 1, chain B and name ca, 12
map_new mapD, gaussian, 1, chain D and name ca, 12
map_new mapF, gaussian, 1, chain F and name ca, 12
map_new mapH, gaussian, 1, chain H and name ca, 12
map_new mapJ, gaussian, 1, chain J and name ca, 12
map_new mapL, gaussian, 1, chain L and name ca, 12
map_new mapN, gaussian, 1, chain N and name ca, 12
map_new mapP, gaussian, 1, chain P and name ca, 12

map_new map1, gaussian, 1, chain 1 and name ca, 12
map_new map3, gaussian, 1, chain 3 and name ca, 12
map_new map5, gaussian, 1, chain 5 and name ca, 12
map_new map7, gaussian, 1, chain 7 and name ca, 12
map_new map9, gaussian, 1, chain 9 and name ca, 12
map_new mapC, gaussian, 1, chain C and name ca, 12
map_new mapG, gaussian, 1, chain G and name ca, 12
map_new mapK, gaussian, 1, chain K and name ca, 12
map_new mapO, gaussian, 1, chain O and name ca, 12
map_new mapQ, gaussian, 1, chain Q and name ca, 12
map_new mapS, gaussian, 1, chain S and name ca, 12
map_new mapU, gaussian, 1, chain U and name ca, 12
map_new mapW, gaussian, 1, chain W and name ca, 12
map_new mapY, gaussian, 1, chain Y and name ca, 12
map_new mapa, gaussian, 1, chain a and name ca, 12
map_new mapc, gaussian, 1, chain c and name ca, 12
map_new mape, gaussian, 1, chain e and name ca, 12
map_new mapg, gaussian, 1, chain g and name ca, 12
map_new mapi, gaussian, 1, chain i and name ca, 12
map_new mapk, gaussian, 1, chain k and name ca, 12
map_new mapm, gaussian, 1, chain m and name ca, 12
map_new mapo, gaussian, 1, chain o and name ca, 12
map_new mapq, gaussian, 1, chain q and name ca, 12
map_new maps, gaussian, 1, chain s and name ca, 12

isosurface surfB, mapB
isosurface surfD, mapD
isosurface surfF, mapF
isosurface surfH, mapH
isosurface surfJ, mapJ
isosurface surfL, mapL
isosurface surfN, mapN
isosurface surfP, mapP

isosurface surf1, map1
isosurface surf3, map3
isosurface surf5, map5
isosurface surf7, map7
isosurface surf9, map9
isosurface surfC, mapC
isosurface surfG, mapG
isosurface surfK, mapK
isosurface surfO, mapO
isosurface surfQ, mapQ
isosurface surfS, mapS
isosurface surfU, mapU
isosurface surfW, mapW
isosurface surfY, mapY
isosurface surfa, mapa
isosurface surfc, mapc
isosurface surfe, mape
isosurface surfg, mapg
isosurface surfi, mapi
isosurface surfk, mapk
isosurface surfm, mapm
isosurface surfo, mapo
isosurface surfq, mapq
isosurface surfs, maps

color cyan, surfB
color cyan, surfD
color cyan, surfF
color cyan, surfH
color cyan, surfJ
color cyan, surfL
color cyan, surfN
color cyan, surfP


color magenta, surf1
color magenta, surf3
color magenta, surf5
color magenta, surf7
color magenta, surf9
color magenta, surfC
color magenta, surfG
color magenta, surfK
color magenta, surfO
color magenta, surfQ
color magenta, surfS
color magenta, surfU
color magenta, surfW
color magenta, surfY
color magenta, surfa
color magenta, surfc
color magenta, surfe
color magenta, surfg
color magenta, surfi
color magenta, surfk
color magenta, surfm
color magenta, surfo
color magenta, surfq
color magenta, surfs

set antialias, 2
set ray_trace_gain, 0.4
set ray_shadows, 0
set specular, 0
show surface, surf*