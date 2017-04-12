remove (hydro)
bg_color white
hide everything
set surface_quality, 1
alter all, b=50
alter all, q=1
set gaussian_resolution, 12

map_new mapA, gaussian, 1, chain A and name ca, 12
map_new mapC, gaussian, 1, chain B and name ca, 12

isosurface surfA, mapA
isosurface surfC, mapC

color cyan, surfA
color magenta, surfC

set antialias, 2
set ray_trace_gain, 0.4
set ray_shadows, 0
set specular, 0
show surface, surf*