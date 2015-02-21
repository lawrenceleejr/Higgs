#!/usr/bin/env python

from math import *

import numpy as np
import math

MZ       = 91.154
MW       = 80.370
# MC       = 1.42
# MB       = 4.49
# MT       = 172.5
MTAU     = 1.77
MMUON    = 0.105658367


MC       = 1.64
MB       = 4.87
MT       = 177.1

VHAT = 246.

g = 2*MW/VHAT
gp = sqrt(  (2*MZ/VHAT)**2 - (2*MW/VHAT)**2  )

YC = MC*sqrt(2)/VHAT
YB = MB*sqrt(2)/VHAT
YT = MT*sqrt(2)/VHAT
YTAU = MTAU*sqrt(2)/VHAT
YMUON = MMUON*sqrt(2)/VHAT

VARRAY = [100.,150.,200.,250.,300.,350.,400.,450.,500.,1000.,2000.,5000.,10000.,100000.,1000000.]
VARRAY = [50.,100.,150.,200.,210.,220.,230.,240.,246.,250.,260.,270.,280.,290.,300.]

# VARRAY = [246.]
#VARRAY = [246.]+[10.*x for x in xrange(3,100) ]+[10.**x for x in xrange(3,30)]

VARRAY = [0.05*x for x in xrange(1,200) ]

# print "V,MZ,MW,MC,MB,MT,MTAU,MMUON,GF"
# for V in VARRAY:
# 	print "%.7g,%.7g,%.7g,%.7g,%.7g,%.7g,%.7g,%.7g,%.7g"%(
# 		V,
# 		sqrt(g**2+gp**2)*V/2.,
# 		g*V/2.,
# 		YC*V/sqrt(2),
# 		YB*V/sqrt(2),
# 		YT*V/sqrt(2),
# 		YTAU*V/sqrt(2),
# 		YMUON*V/sqrt(2),
# 		(1/V**2)*(246.**2)*0.000011663
# 		)


print "V,MZ,MW,MC,MB,MT,MTAU,MMUON,GF,GAMW,GAMZ,ALPHAS"
for V in VARRAY:
	# print V
	print "%.7g,%.7g,%.7g,%.7g,%.7g,%.7g,%.7g,%.7g,%.7g,%.7g,%.7g,%.7g"%(
		V,
		MZ*V,
		MW*V,
		MC*V,
		MB*V,
		MT*V,
		MTAU*V,
		MMUON*V,
		(1/V**2)*0.000011663,
		2.08856*V,
		2.49581*V,
		 1./( (1./0.118)+(1./(2.*np.pi))*(-23./3.)*math.log(1./V) )
		)


const =  1./( (1./0.118)+(1./(2.*np.pi))*(-23./3.)*math.log(1./V) )
