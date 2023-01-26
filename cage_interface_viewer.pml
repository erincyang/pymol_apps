remove hydro
hide all
sele Aint1, br. chain A around 8 and not chain C+E+F+G+I+K+M+N+O+Q+S+U+W+Y
sele Bint1, br. Aint1 around 8 and chain A
sele interface, Aint1 or Bint1
show lines, interface
show ribbon
center interface
sele comp1, chain A+C+E+G+I+K+M+O+Q+S+U+W+Y+a+c+e+g+i+k+m+o+q+s+u+w+y+0+2+4+6+8
sele comp2, chain B+D+F+H+J+L+N+P+R+T+V+X+Z+b+d+f+h+j+l+n+p+r+t+v+x+z+1+3+5+7+9
#problem chains for comp1
sele chainp1, chain !
sele chainp2, chain &
sele chainp3, chain #
sele chainp4, chain <
sele comp1, comp1 or chainp1 or chainp2 or chainp3 or chainp4
sele chainp1, chain $
sele chainp2, chain .
sele chainp3, chain @
sele chainp4, chain >
sele comp2, comp2 or chainp1 or chainp2 or chainp3 or chainp4

delete chainp1
delete chainp2
delete chainp3
delete chainp4

set_color ipd_magenta=[0.710 , 0.118 , 0.518]
set_color ipd_blue=[0.263 , 0.592 , 0.710]
set_color ipd_purple=[0.216 , 0.141 , 0.420]
color ipd_magenta, comp1 and name C*
color ipd_blue, comp2 and name C*

