#!/usr/bin/env python

from math import *

import numpy as np
import math

MZ       = 91.154
MW       = 80.370
MTAU     = 1.776
MMUON    = 0.105658367

#NNLO Pole Masses from HDECAY
MC       = 1.64
MB       = 4.87
MT       = 177.1



GF = 0.000011663

GAMW  = 2.08856
GAMZ  = 2.49581

ALPHAS = 0.118

# VHAT = 246.

# g = 2*MW/VHAT
# gp = sqrt(  (2*MZ/VHAT)**2 - (2*MW/VHAT)**2  )

# YC = MC*sqrt(2)/VHAT
# YB = MB*sqrt(2)/VHAT
# YT = MT*sqrt(2)/VHAT
# YTAU = MTAU*sqrt(2)/VHAT
# YMUON = MMUON*sqrt(2)/VHAT


# VARRAY = [0.05*x for x in xrange(1,200) ]
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
		(1/V**2)*GF, ## Should there also be a (1/math.sqrt(2) )?
		GAMW*V,
		GAMZ*V,
		 1./( (1./ALPHAS)+(1./(2.*np.pi))*(-23./3.)*math.log(1./V) )
		)


# const =  1./( (1./0.118)+(1./(2.*np.pi))*(-23./3.)*math.log(1./V) )
