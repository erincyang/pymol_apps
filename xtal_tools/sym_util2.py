	# -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
import sys,os,inspect,functools
newpath = os.path.dirname(inspect.getfile(inspect.currentframe())) # script directory
if not newpath in sys.path: sys.path.append(newpath)
import string,re,gzip,itertools
from pymol_util import *
import operator as op

def get_xforms_by_chain(sele="all",verbose=False,userms=False):
	v = cmd.get_view()
	cen = com("("+sele+") and (name CA and not HET)")
	chains = cmd.get_chains(sele)
	xforms = dict()
	maxrms = 0.0
	for c1,c2 in filter( lambda t: True, product(chains,chains) ):
		refsele = "((%s) and chain %s and name CA and not HET)"%(sele,c1)
		movsele = "((%s) and chain %s and name CA and not HET)"%(sele,c2)
		if userms: x,rms = getrelframe_rmsalign( movsele, refsele, Xform(-cen) )
		else:      x,rms = getrelframe( movsele, refsele, Xform(-cen)), 0.0
		maxrms = max(maxrms,rms)
		xforms[c1,c2] = x
		#if c1 in "AB" and c2 in "AB":
		#	print movsele
		#	print refsele
		#	print x.pretty()
		#	print
	# raise Exception
	cmd.set_view(v)
	return xforms, maxrms

def find_symelems(sele_or_xforms="all",verbose=False):
	xforms = sele_or_xforms
	if isinstance(sele_or_xforms,basestring): xforms, maxrms = get_xforms_by_chain(sele_or_xforms,verbose=True)
	elif not isinstance(sele_or_xforms,dict): raise ValueError
	symelems = list()
	maxangerr = 0.0
	for c,x in xforms.items():
		assert len(c)==2
		assert isinstance(x,Xform)
		if c[0]==c[1]: continue
		dis = x.t.length()
		if dis > 5.0: continue
		axis,ang = x.rotation_axis()
		nfold = round(math.pi*2.0/ang)
		angerr = abs(ang-math.pi*2.0/nfold)*180.0/math.pi
		if verbose: print("candidate symelem:",nfold, c, angerr, axis)
		if angerr > 360.0/nfold/8.0: continue # require unambiguous symelems
		maxangerr = max(maxangerr,angerr*nfold)
		symelems.append( (nfold,axis,c,angerr) )
	symelemdis = lambda x,y: line_line_angle_degrees(x[1],y[1]) if x[0]==y[0] else 9e9
	if verbose:
		for se1,se2 in filter( lambda t: t[0]<t[1], product(symelems,symelems) ):
			if se1[0]==se2[0]:
				print(se1)
				print(se)
				print(symelemdis(se1,se2), "degrees")
				#print
	hier = HierarchicalClustering(symelems, symelemdis )
	thresh = 6.0
	clusters = hier.getlevel(thresh);
	print("number of symmetry element clusters at threshold",thresh,"degrees is" , len(clusters))
	centers0 = list()
	maxaxiserr = 0.0
	for clust in clusters:
		print("symelem cluster:",clust)
		center = list(clust[0])
		center[2] = list((center[2],))
		for i in range(1,len(clust)):
			ax = clust[i][1]
			center[1] = center[1] + ( ax if ax.dot(center[1]) > 0 else -ax )
			center[2].append(clust[i][2])
			center[3] = max(center[3],clust[i][3])
		center[1].normalize()
		centers0.append(center)
		axiserr = 0.0
		for c in clust:	axiserr = max( axiserr, 1.0-abs(center[1].dot(c[1])) )
		maxaxiserr = max(maxaxiserr,axiserr)
	# sort on nfold, then on number of chain pairs in cluster
	centers0 = sorted( centers0, cmp = lambda x,y: cmp(y[0],x[0]) if x[0]!=y[0] else cmp(len(y[2]),len(x[2])) )
	centers = list()
	for center in centers0:
		if verbose: print("DEBUG prune center:",center)
		seenit = False
		for censeen in centers:
			remainder = abs( ( censeen[0] / center[0] ) % 1.0)
			if verbose: print("   ",remainder,censeen)
			if remainder > 0.01: continue # not a symmetry multiple
			if 1.0-abs(center[1].dot(censeen[1])) < 0.01:
				seenit = True # axis are same
		if not seenit:
			centers.append(center)
	print("centers:")
	cen_of_geom = com("("+sele_or_xforms+") and (name CA and not HET)")
	for center in centers:
		print(center)
		# if center[0]>2.1: continue
		#showvecfrompoint(50*center[1],cen_of_geom)
	return centers, maxrms, maxangerr, maxaxiserr

def guessdxaxes(sele="all",verbose=False):
	nfold = len(cmd.get_chains(sele))
	assert nfold % 2 is 0
	nfold /= 2
	symelems, maxrms, angerr, axiserr = find_symelems(sele,verbose=verbose)
	assert len(symelems) > 1
	assert symelems[0][0] == float(nfold)
	assert symelems[1][0] == float(2)
	axis_high = symelems[0][1]
	axis_low  = symelems[1][1]
	return axis_high, axis_low, maxrms, angerr, axiserr

def aligndx(sele='all',verbose=False):
	trans(sele,-com(sele+" and name CA and not HET"))
	haxis, laxis, maxrms, angerr, axiserr = guessdxaxes(sele,verbose=verbose)
	xalign = alignvectors( haxis, laxis, Uz, Ux )
	xform(sele,xalign)
	return maxrms, angerr, axiserr

def guesscxaxis(sele,nfold,chains=None):	
	sele = "(("+sele+") and (name CA))"
	if not chains: chains = cmd.get_chains(sele)
	if len(chains) != nfold: 
		print("num chains != n-fold")
		return None
	atoms = cmd.get_model(sele).atom
	idx = {}
	for i,c in enumerate(chains): idx[c] = i
	coords = [list() for c in chains]
	for a in atoms:
		if a.chain in idx:
			coords[idx[a.chain]].append(Vec(a.coord))
	if max(len(x) for x in coords) != min(len(x) for x in coords):
		print("chains not same length",)
		return None
	cen = reduce(op.add,(x for c in coords for x in c),V0) / len(coords[0]*len(coords))
	guesses = []
	axis = Vec(0,0,0)
	diserr = 0
	for xyzs in izip(*coords):
		a = reduce(op.add,(xyz-cen for xyz in xyzs))/len(xyzs)
		da = a if a.dot(Uz) > 0 else -a
		axis += da
		guesses.append(da)
		dis = [a.distance(cen) for a in xyzs]
		avgdis = sum(dis)/len(dis)
		diserr += sum((x-avgdis)**2 for x in dis)
	axis.normalize()
	diserr = math.sqrt(diserr/len(coords)/len(coords[0]))
	angerr = math.sqrt(sum(projperp(axis,a).length_squared() for a in guesses)/len(guesses))
	return axis,cen,diserr,angerr

def aligncx(sele,nfold,alignsele=None,tgtaxis=Uz,chains=None):
	if not alignsele: alignsele = sele
	tmp = guesscxaxis(sele,nfold,chains)
	if not tmp: return None
	axis,cen,diserr,angerr = guesscxaxis(sele,nfold,chains)
	trans(sele,-cen)
	alignaxis(sele,tgtaxis,axis,xyz.Vec(0,0,0))
	return tmp

guessc2axis = functools.partial(guesscxaxis,nfold=2)
guessc3axis = functools.partial(guesscxaxis,nfold=3)
guessc4axis = functools.partial(guesscxaxis,nfold=4)
guessc5axis = functools.partial(guesscxaxis,nfold=5)
guessc6axis = functools.partial(guesscxaxis,nfold=6)
alignc2 = functools.partial(aligncx,nfold=2)
alignc3 = functools.partial(aligncx,nfold=3)
alignc4 = functools.partial(aligncx,nfold=4)
alignc5 = functools.partial(aligncx,nfold=5)
alignc6 = functools.partial(aligncx,nfold=6)

# def c2axis(sele,alignsele=None,chains=["A","B"]):
# 	if alignsele is None: alignsele = sele
#         #cmd.create('tmp98367598',sele)
#         #sele = 'tmp98367598'
# 	cmd.remove(sele+" and resn HOH")
# 	cen = com(alignsele)
# 	a = cmd.get_model(alignsele+" and chain "+chains[0]+" and name CA").atom
# 	b = cmd.get_model(alignsele+" and chain "+chains[1]+" and name CA").atom
# 	# print len(a),len(b)
# 	if len(a) != len(b) or len(a) == 0:
# 		print "ERROR on %s: subunits are not the same size!"%alignsele
# 		return False
# 	axis = xyz.Vec(0,0,0)
# 	for i in range(len(a)):
# 		axis1 = ( xyz.Vec(a[i].coord)+xyz.Vec(b[i].coord)-2*cen ).normalized()
# 		if axis.length() > 0.0001 and axis.dot(axis1) < 0:
# 			axis1 *= -1
# 		axis += axis1
# 		# print axis1
# 	axis.normalize()
#         #cmd.delete('tmp98367598')
# 	return axis

# def alignc2(sele,alignsele=None,tgtaxis=xyz.Vec(0,0,1),chains=["A","B"]):
# 	if alignsele is None: alignsele = sele
# 	axis = c2axis(sele,alignsele,chains)
# 	if not axis: return -1
# 	# print "axis of rotation:",axis
# 	alignaxis(sele,tgtaxis,axis,xyz.Vec(0,0,0))	
# 	# seleA = "("+alignsele+") and chain %s"%chains[0]
# 	# seleB = "("+alignsele+") and chain %s"%chains[1]
# 	# print seleA
# 	# print seleB
# 	# rot(seleA,tgtaxis,180.0)
# 	# r = cmd.rms_cur(seleA,seleB)
# 	# rot(seleA,tgtaxis,180.0)	
# 	return 0

# def c3axis(sele,alignsele=None,chains=["A","B","C"]):
# 	if alignsele is None: alignsele = sele
#         #cmd.create('tmp98367598',sele)
#         #sele = 'tmp98367598'
# 	cmd.remove(sele+" and resn HOH")
# 	cen = com(alignsele)
# 	a = cmd.get_model(alignsele+" and chain "+chains[0]+" and name CA").atom
# 	b = cmd.get_model(alignsele+" and chain "+chains[1]+" and name CA").atom
# 	c = cmd.get_model(alignsele+" and chain "+chains[2]+" and name CA").atom
# 	# print "subunit lengths:",len(a),len(b),len(c)
# 	if len(a) != len(b) or len(a) != len(c) or len(a) == 0:
# 		print "ERROR on %s: subunits are not the same size!"%alignsele
# 		return False
# 	axis = xyz.Vec(0,0,0)
# 	for i in range(len(a)):
# 		axis1 = ( xyz.Vec(a[i].coord)+xyz.Vec(b[i].coord)+xyz.Vec(c[i].coord) - 3*cen ).normalized()
# 		if axis.length() > 0.0001 and axis.dot(axis1) < 0:
# 			axis1 *= -1
# 		axis += axis1
# 		# print axis1
# 	axis.normalize()
#         #cmd.delete('tmp98367598')
# 	return axis

# def alignc3(sele,alignsele=None,tgtaxis=xyz.Vec(0,0,1),chains=["A","B","C"]):
# 	if alignsele is None: alignsele = sele
# 	cmd.remove(sele+" and resn HOH")
# 	axis = c3axis(sele,alignsele,chains)
# 	# print "axis of rotation:",axis
# 	alignaxis(sele,tgtaxis,axis,xyz.Vec(0,0,0))
# 	return True

# def c4axis(sele,alignsele=None,chains=["A","B","C","D"]):
# 	if alignsele is None: alignsele = sele
# 	cmd.remove(sele+" and resn HOH")
# 	trans(sele,-com(alignsele))
# 	a = cmd.get_model(alignsele+" and chain "+chains[0]+" and name CA").atom
# 	b = cmd.get_model(alignsele+" and chain "+chains[1]+" and name CA").atom
# 	c = cmd.get_model(alignsele+" and chain "+chains[2]+" and name CA").atom
# 	d = cmd.get_model(alignsele+" and chain "+chains[3]+" and name CA").atom
# 	# print "subunit lengths:",len(a),len(b),len(c)
# 	if len(a) != len(b) or len(a) != len(c) or len(a) == 0 or len(d) != len(a):
# 		print "ERROR on %s: subunits are not the same size!"%alignsele
# 		return False
# 	axis = xyz.Vec(0,0,0)
# 	for i in range(len(a)):
# 		axis1 = ( xyz.Vec(a[i].coord)+xyz.Vec(b[i].coord)+xyz.Vec(c[i].coord)+xyz.Vec(d[i].coord) ).normalized()
# 		if axis.length() > 0.0001 and axis.dot(axis1) < 0:
# 			axis1 *= -1
# 		axis += axis1
# 		# print axis1
# 	axis.normalize()
# 	return axis

# def alignc4(sele,alignsele=None,tgtaxis=xyz.Vec(0,0,1),chains=["A","B","C","D"]):
# 	if alignsele is None: alignsele = sele
# 	cmd.remove(sele+" and resn HOH")
# 	axis = c3axis(sele,alignsele,chains)
# 	# print "axis of rotation:",axis
# 	alignaxis(sele,tgtaxis,axis,xyz.Vec(0,0,0))
# 	return True

# def c5axis(sele,alignsele=None,chains=["A","B","C","D","E"]):
# 	if alignsele is None: alignsele = sele
# 	cmd.remove(sele+" and resn HOH")
# 	trans(sele,-com(alignsele))
# 	a = cmd.get_model(alignsele+" and chain "+chains[0]+" and name CA").atom
# 	b = cmd.get_model(alignsele+" and chain "+chains[1]+" and name CA").atom
# 	c = cmd.get_model(alignsele+" and chain "+chains[2]+" and name CA").atom
# 	d = cmd.get_model(alignsele+" and chain "+chains[3]+" and name CA").atom
# 	e = cmd.get_model(alignsele+" and chain "+chains[4]+" and name CA").atom
# 	# print "subunit lengths:",len(a),len(b),len(c)
# 	if len(a) != len(b) or len(a) != len(c) or len(a) == 0:
# 		print "ERROR on %s: subunits are not the same size!"%alignsele
# 		return False
# 	axis = xyz.Vec(0,0,0)
# 	for i in range(len(a)):
# 		axis1 = ( xyz.Vec(a[i].coord)+xyz.Vec(b[i].coord)+xyz.Vec(c[i].coord)+xyz.Vec(d[i].coord)+xyz.Vec(e[i].coord) ).normalized()
# 		if axis.length() > 0.0001 and axis.dot(axis1) < 0:
# 			axis1 *= -1
# 		axis += axis1
# 		# print axis1
# 	axis.normalize()
# 	return axis

# def alignc5(sele,alignsele=None,tgtaxis=xyz.Vec(0,0,1),chains=["A","B","C","D","E"]):
# 	if alignsele is None: alignsele = sele
# 	cmd.remove(sele+" and resn HOH")
# 	axis = c5axis(sele=sele,alignsele=alignsele,chains=chains)
# 	print "axis of rotation:",axis,"to",tgtaxis
# 	alignaxis(sele,tgtaxis,axis,xyz.Vec(0,0,0))
# 	return True

def myint(s):
   i = len(s)
   while i > 0 and not s[:i].isdigit(): i -= 1
   if not i: return None
   return int(s[:i])

def mki213(N, sel = 'all'):
	v = cmd.get_view()
	cmd.delete("i213_*")
	cmd.delete('base80345769083457')
	cmd.delete('tmp80345769083457')
	c2 = com(sel)
	c3 = xyz.Vec(0, 0, 0)
	cmd.create( 'tmp80345769083457', sel)
	a2 = c2axis('tmp80345769083457')
	cmd.delete( 'tmp80345769083457')
	a3 = xyz.Vec(0, 0, 1)
	cmd.create('base80345769083457', sel+" and chain A and visible")
	seenit = []
	R2 = [xyz.rotation_matrix_degrees(a2, 0), xyz.rotation_matrix_degrees(a2, 180), ]
	R3 = [xyz.rotation_matrix_degrees(a3, 0), xyz.rotation_matrix_degrees(a3, 120), xyz.rotation_matrix_degrees(a3, 240), ]
	C = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	print(a2, c2, a3, c3)
	for i21 in range(2):
		for i32 in range(3 if N > 1 else 1):
			for i22 in range(2 if N > 2 else 1):
				for i33 in range(3 if N > 3 else 1):
					for i23 in range(2 if N > 4 else 1):
						for i34 in range(3 if N > 5 else 1):
							for i24 in range(2 if N > 6 else 1):
								for i35 in range(3 if N > 7 else 1):
									for i25 in range(2 if N > 8 else 1):
										test = xyz.Vec(0, 0, 0)
										test = R2[i21]*(test-c2)+c2
										test = R3[i32]*(test-c3)+c3
										test = R2[i22]*(test-c2)+c2
										test = R3[i33]*(test-c3)+c3
										test = R2[i23]*(test-c2)+c2
										test = R3[i34]*(test-c3)+c3
										test = R2[i24]*(test-c2)+c2
										test = R3[i35]*(test-c3)+c3
										test = R2[i25]*(test-c2)+c2
										#print test
										seen = False
										for xs in seenit:
											if (xs-test).length() < 0.1:
												seen = True
												break
										if seen: continue
										else: seenit.append(test)
										n = "i213_%i%i%i%i%i%i%i%i%i"%(i25, i35, i24, i34, i23, i33, i22, i32, i21)
										cmd.create(n, 'base80345769083457')
										rot(n, a2, i21*180.0, c2)
										rot(n, a3, i32*120.0, c3)
										rot(n, a2, i22*180.0, c2)
										rot(n, a3, i33*120.0, c3)
										rot(n, a2, i23*180.0, c2)
										rot(n, a3, i34*120.0, c3)
										rot(n, a2, i24*180.0, c2)
										rot(n, a3, i35*120.0, c3)
										rot(n, a2, i25*180.0, c2)
	print(len(seenit))
	cmd.delete('base80345769083457')
	cmd.set_view(v)

def viewi213(sel = "all"):
	cmd.hide('ev')
	cmd.show('rib')
	mki213(sel)
	cmd.show('car', 'not i213*')
	cmd.hide('rib', 'not i213*')
	cmd.show('lines', '(byres (%s and not i213* and chain A) within 7.0 of (%s and not i213* and chain B))'%(sel, sel))
	cmd.show('lines', '(byres (%s and not i213* and chain B) within 7.0 of (%s and not i213* and chain A))'%(sel, sel))









def mkp23(N, R=43.5, i=0, sel = 'all'):
	v = cmd.get_view()
	cmd.delete("p23_*")
	cmd.delete('base80345769083457')
	cmd.delete('tmp80345769083457')
	c2 = xyz.Vec(0, 0, 0)
	c3 = xyz.Vec(R,R,-R)
	cmd.create( 'tmp80345769083457', sel)
	cmd.delete( 'tmp80345769083457')
	a3 = [xyz.Vec(0,0,0),xyz.Vec(1,1,1),xyz.Vec(-1,-1,-1)]
	a2 = [xyz.Vec(0,0,0),xyz.Vec(1,0,0),xyz.Vec(0,1,0),xyz.Vec(0,0,1)]
	cmd.create('base80345769083457', sel+" and visible")
	seenit = []
	R2 = [xyz.rotation_matrix_degrees(a2[1],  0), # hack
	      xyz.rotation_matrix_degrees(a2[1],180),
	      xyz.rotation_matrix_degrees(a2[2],180),
	      xyz.rotation_matrix_degrees(a2[3],180) ]
	R3 = [xyz.rotation_matrix_degrees(a3[1],  0), # hack!
	      xyz.rotation_matrix_degrees(a3[1],120),
	      xyz.rotation_matrix_degrees(a3[2],120), ]
	C = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	print(a2, c2, a3, c3)
	for i21 in range(4):
		for i32 in range(3 if N > 1 else 1):
			for i22 in range(4 if N > 2 else 1):
				for i33 in range(3 if N > 3 else 1):
					for i23 in range(4 if N > 4 else 1):
						for i34 in range(3 if N > 5 else 1):
							for i24 in range(4 if N > 6 else 1):
								for i35 in range(3 if N > 7 else 1):
									for i25 in range(4 if N > 8 else 1):
										test = xyz.Vec(1, 1, 1)
										test = R2[i21]*(test-c2)+c2
										test = R3[i32]*(test-c3)+c3
										test = R2[i22]*(test-c2)+c2
										test = R3[i33]*(test-c3)+c3
										test = R2[i23]*(test-c2)+c2
										test = R3[i34]*(test-c3)+c3
										test = R2[i24]*(test-c2)+c2
										test = R3[i35]*(test-c3)+c3
										test = R2[i25]*(test-c2)+c2
										#print test
										seen = False
										for xs in seenit:
											if (xs-test).length() < 0.1:
												seen = True
												break
										if seen: continue
										else: seenit.append(test)
										n = "p23_%i%i%i%i%i%i%i%i%i"%(i25, i35, i24, i34, i23, i33, i22, i32, i21)
										cmd.create(n, 'base80345769083457')
										if i21 > 0: rot(n, a2[i21], 180.0, c2)
										if i32 > 0: rot(n, a3[i32], 120.0, c3)
										if i22 > 0: rot(n, a2[i22], 180.0, c2)
										if i33 > 0: rot(n, a3[i33], 120.0, c3)
										if i23 > 0: rot(n, a2[i23], 180.0, c2)
										if i34 > 0: rot(n, a3[i34], 120.0, c3)
										if i24 > 0: rot(n, a2[i24], 180.0, c2)
										if i35 > 0: rot(n, a3[i35], 120.0, c3)
										if i25 > 0: rot(n, a2[i25], 180.0, c2)
	print("seen:",len(seenit))
	cmd.delete('base80345769083457')
	cmd.set_view(v)

def selbycomp(trn=0):
	cmd.select("TRI1","TRI and chain A+B+C")
	cmd.select("TRI2","TRI and chain D+E+F")
	cmd.select("TRI3","TRI and chain G+H+I")
	cmd.select("TRI4","TRI and chain J+K+L")
	cmd.select("TRI5","TRI and chain xyz.Mat+N+O")
	cmd.select("TRI6","TRI and chain P+Q+R")
	cmd.select("TRI7","TRI and chain S+T+U")
	cmd.select("TRI8","TRI and chain xyz.Vec+W+Ux")
	cmd.select("DIM1","DIM and chain A+D")
	cmd.select("DIM2","DIM and chain B+G")
	cmd.select("DIM3","DIM and chain C+J")
	cmd.select("DIM4","DIM and chain E+U")
	cmd.select("DIM5","DIM and chain F+R")
	cmd.select("DIM6","DIM and chain H+T")
	cmd.select("DIM7","DIM and chain I+O")
	cmd.select("DIM8","DIM and chain K+Q")
	cmd.select("DIM9","DIM and chain L+N")
	cmd.select("DIM10","DIM and chain xyz.Mat+xyz.Vec")
	cmd.select("DIM11","DIM and chain P+W")
	cmd.select("DIM12","DIM and chain Ux+S")
	cmd.delete("LINE*")
	cmd.delete("serf*")

	cmd.do("""alter all, b=50
	alter all, q=1
	set gaussian_resolution,8""")
	ISO="""map_new map%s, gaussian, 2, %s, 10
	isosurface surf%s, map%s"""


	for i in range(1, 9):
		cmd.do(ISO%(("TRI%i"%i,)*4))
		cmd.color(COLORS[i-1],"surfTRI%i"%i)
		c = com("TRI%i"%i)
		# trans("TRI%i"%i,trn*c.normalized())
		obj = [
			cgo.CYLINDER,
		   	0.0, 0.0, 0.0,
		   	1.6*c.x, 1.6*c.y, 1.6*c.z,
			1.5,
			0.1,0.1,0.1,0.1,0.1,0.1,
		]
		cmd.load_cgo(obj,'LINETRI%i'%i)
	for i in range(1,13):
		cmd.do(ISO%(("DIM%i"%i,)*4))
		cmd.color(COLORS[i+7],"surfDIM%i"%i)
		c = com("DIM%i"%i)
		# trans("DIM%i"%i,trn*com("DIM%i"%i).normalized())
		obj = [
			cgo.CYLINDER,
		   	0.0, 0.0, 0.0,
		   	1.3*c.x, 1.3*c.y, 1.3*c.z,
			1.0,
			0,0,1,0,0,1
		]
		cmd.load_cgo(obj,'LINEDIM%i'%i)

def getframe(obj):
	m = cmd.get_model(obj)
	x = xyz.Vec(m.atom[       0     ].coord)
	y = xyz.Vec(m.atom[len(m.atom)/2].coord)
	z = xyz.Vec(m.atom[      -1     ].coord)
	return xyz.stub(x,y,z)

def getrelframe(newobj,refobj,Forigin=None):
	if Forigin is None: Forigin = xyz.Xform(xyz.Imat,xyz.Vec(0,0,0))
	Fref = getframe(refobj)
	Fnew = getframe(newobj)
	Fdelta = Fnew * ~Fref
	#print (Fdelta * Forigin).pretty()
	return Fdelta * Forigin



def rechain(sel,nres):
	chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890abcdefghijklmnopqrstuvwxyz"
	ntot = len(getres(sel))
	assert ntot % nres == 0
	for i in range(ntot/nres):
		cmd.alter("resi %i-%i"%( nres*i+1,nres*(i+1)),"chain='%s'"%chains[i])

def makecx(sel = 'all', n = 5):
	v = cmd.get_view()
	cmd.delete("C%i_*"%n)
	chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
	for i in range(n): cmd.create("C%i_%i"%(n, i), sel+" and (not C%i_*)"%n)
	for i in range(n): rot("C%i_%i"%(n, i), Uz, 360.0*float(i)/float(n))
	for i in range(n): cmd.alter("C%i_%i"%(n, i), "chain = '%s'"%chains[i])
	util.cbc("C*")
	cmd.set_view(v)

def makedx(sel = 'all', n = 5):
	v = cmd.get_view()
	cmd.delete("D%i_*"%n)
	chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
	for i in range(n):
		dsel  = "D%i_%i"%(n,  i)
		dsel2 = "D%i_%i"%(n,n+i)
		cmd.create(dsel , sel+" and (not D%i_*)"%n)
		rot       (dsel , Uz, 360.0*float(i)/float(n))
		cmd.create(dsel2, dsel )
		rot       (dsel2, Uy, 180.0 )
		cmd.alter (dsel , "chain = '%s'"%chains[i])
		cmd.alter (dsel2, "chain = '%s'"%chains[i+n])
	util.cbc("D*")
	cmd.set_view(v)

for i in range(2,21):
		globals()['makec%i'%i] = partial(makecx,n=i)

for i in range(2,21):
		globals()['maked%i'%i] = partial(makedx,n=i)

def makecxauto():
	for o in cmd.get_object_list():
		n = int(re.search("_C\d+_", o).group(0)[2:-1])
		makecx(o, n)

def makekinwire(sel,movres,fixres):
	v = cmd.get_view()
	cmd.delete("ha"); cmd.create("ha",sel); cmd.alter("ha","chain='A'")
	cmd.delete("hb"); cmd.create("hb",sel); cmd.alter("hb","chain='B'")
	cmd.delete("hc"); cmd.create("hc",sel); cmd.alter("hc","chain='C'")
	cmd.delete("hd"); cmd.create("hd",sel); cmd.alter("hd","chain='D'")
	cmd.delete("he"); cmd.create("he",sel); cmd.alter("he","chain='E'")
	cmd.align("hb and resi %i"%movres,"ha and resi %i"%fixres);
	cmd.align("hc and resi %i"%movres,"hb and resi %i"%fixres);
	cmd.align("hd and resi %i"%movres,"hc and resi %i"%fixres);
	cmd.align("he and resi %i"%movres,"hd and resi %i"%fixres);
	util.cbc('elem C')
	v = cmd.set_view(v)

def get_contigs(x,n=7):
	"""
	>>> test = list(range(1,8)) + list(range(20,33)) + list(range(40,44)) + list(range(49,50))+ list(range(0,8))
	>>> print test
	[1, 2, 3, 4, 5, 6, 7, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 40, 41, 42, 43, 49, 0, 1, 2, 3, 4, 5, 6, 7]
	
	>>> print get_contigs( test )
	[[1, 2, 3, 4, 5, 6, 7], [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32], [0, 1, 2, 3, 4, 5, 6, 7]]
	"""
	if type(x) is type(""):
		x = [x[1] for x in getres(x)]
	x.append(-123456789)
	contigs = [[],]
	for i in range(len(x)-1):
		contigs[-1].append(x[i])
		if x[i]+1 is not x[i+1]:
			contigs.append(list())
	return [c for c in contigs if len(c) >= n]

def get_contigs_termini(x,n=7):
	"""
	>>> test = list(range(1,8)) + list(range(20,33)) + list(range(40,44)) + list(range(49,50))+ list(range(0,8))
	>>> print test
	[1, 2, 3, 4, 5, 6, 7, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 40, 41, 42, 43, 49, 0, 1, 2, 3, 4, 5, 6, 7]
	
	>>> print get_contigs_termini( test )
	[[1, 2, 3, 4, 5, 6, 7], [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32], [0, 1, 2, 3, 4, 5, 6, 7]]
	"""
	pass #cend = []

def get_fixed_size_contigs(x,n=7):
	"""
	>>> test = list(range(1,8)) + list(range(20,33)) + list(range(40,44)) + list(range(49,50))+ list(range(0,8))
	>>> print(test)
	[1, 2, 3, 4, 5, 6, 7, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 40, 41, 42, 43, 49, 0, 1, 2, 3, 4, 5, 6, 7]
	
	>>> for f in get_fixed_size_contigs(test,7): print(f)
	[1, 2, 3, 4, 5, 6, 7]
	[20, 21, 22, 23, 24, 25, 26]
	[21, 22, 23, 24, 25, 26, 27]
	[22, 23, 24, 25, 26, 27, 28]
	[23, 24, 25, 26, 27, 28, 29]
	[24, 25, 26, 27, 28, 29, 30]
	[25, 26, 27, 28, 29, 30, 31]
	[26, 27, 28, 29, 30, 31, 32]
	[0, 1, 2, 3, 4, 5, 6]
	[1, 2, 3, 4, 5, 6, 7]

	>>> for f in get_fixed_size_contigs(test,9): print(f)
	[20, 21, 22, 23, 24, 25, 26, 27, 28]
	[21, 22, 23, 24, 25, 26, 27, 28, 29]
	[22, 23, 24, 25, 26, 27, 28, 29, 30]
	[23, 24, 25, 26, 27, 28, 29, 30, 31]
	[24, 25, 26, 27, 28, 29, 30, 31, 32]

	>>> print(len(get_fixed_size_contigs(test,1)))
	28

	>>> for f in get_fixed_size_contigs(test,4): print(f)
	[1, 2, 3, 4]
	[2, 3, 4, 5]
	[3, 4, 5, 6]
	[4, 5, 6, 7]
	[20, 21, 22, 23]
	[21, 22, 23, 24]
	[22, 23, 24, 25]
	[23, 24, 25, 26]
	[24, 25, 26, 27]
	[25, 26, 27, 28]
	[26, 27, 28, 29]
	[27, 28, 29, 30]
	[28, 29, 30, 31]
	[29, 30, 31, 32]
	[0, 1, 2, 3]
	[1, 2, 3, 4]
	[2, 3, 4, 5]
	[3, 4, 5, 6]
	[4, 5, 6, 7]
	"""
	f = []
	for c in get_contigs(x):
		for i in range(0,len(c)-n+1):
			f.append(range(c[i],c[i]+n))
	return f

def tmpname():
	return "TEMPORARY_"+str(random.random())

def gen_helical_alignments(sele1,sele2,pref="HALN"):
	cmd.delete(pref+"_*")
	cmd.alter('(%s) or (%s)'%(sele1,sele2),'resn="ALA"')
	cmd.alter(sele2+' and chain A','chain="Z"')
	cmd.alter(sele2+' and chain B','chain="Y"')
	cmd.alter(sele2+' and chain C','chain="X"')
	cmd.alter(sele2+' and chain D','chain="W"')
	cmd.alter(sele2+' and chain E','chain="V"')
	cmd.alter(sele2+' and chain F','chain="U"')
	cmd.alter(sele2+' and chain G','chain="T"')
	cmd.alter(sele2+' and chain H','chain="S"')
	cmd.alter(sele2+' and chain I','chain="R"')
	cmd.alter(sele2+' and chain J','chain="Q"')
	chunks1 = ["chain A and resi "+str(h)[1:-1].replace(', ','+') for h in get_fixed_size_contigs("chain A and %s and ss H"%sele1)]
	chunks2 = ["chain Z and resi "+str(h)[1:-1].replace(', ','+') for h in get_fixed_size_contigs("chain Z and %s and ss H"%sele2)]
	for i,hsel1 in enumerate(chunks1):
		name1 = pref+"_"+tmpname()
		algn1 = name1+" and "+hsel1
		cmd.create(name1,sele1)
		for j,hsel2 in enumerate(chunks2):
			name2 = pref+"_"+tmpname()
			algn2 = name2+" and "+hsel2
			cmd.create(name2,sele2)
			print(algn2,algn1)
			print(name1+" and chain A and "+hsel1)
			print(name2+" and chain Z and "+hsel2)
			# now align them
			cmd.align( name2+" and chain Z and "+hsel2, name1+" and chain A and "+hsel1 )
			name3 = pref+"_%03i_%03i"%(i,j)
			cmd.create(name3,name1+" or "+name2)
			util.cbc(name3+" and elem C")
			cmd.delete(name2)
		cmd.delete(name1)



def nulltest():
	"""
	>>> print("foo")
	foo
	"""
	return None

def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite())
    return tests

if __name__ == '__main__':
   import doctest
   r = doctest.testmod()
   print(r)


