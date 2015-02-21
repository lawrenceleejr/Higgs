#!/usr/bin/env python

import sys,os
import numpy as np
import shutil
import re
import random

from ROOT import *
gROOT.LoadMacro("atlasstyle/AtlasStyle.C")
# SetAtlasStyle()

# from rootpy.interactive import wait
from rootpy.plotting import Canvas, Hist, Hist2D, Hist3D, Graph, Graph2D


import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from decimal import Decimal

def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))




PDFArray = {}
PDFArray['MH'] = []
PDFArray['V'] = []
PDFArray['P'] = []

MHhatArray = {} #for each vev, I want an array of MhatValues
MHhatArray['MHhat'] = []
MHhatArray['V'] = []
MHhatArray['P'] = []


SMParams = [
'ALS(MZ)',
'MSBAR(2)',
'MC',
'MB',
'MT',
'MTAU',
'MMUON',
'1/ALPHA',
'GF',
'GAMW',
'GAMZ',
'MZ',
'MW',
'VTB',
'VTS',
'VTD',
'VCB',
'VCS',
'VCD',
'VUB',
'VUS',
'VUD'
]

	
VEVSamplePoints = [

#### Nominal 
# [246.000000,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],


# #### ALS(MZ)
# [0.06  ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],
# [0.07  ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],
# [0.08  ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],
# [0.09  ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],
# [0.1   ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],
# [0.11  ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],
# [0.115 ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],
# [0.116 ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],
# [0.117 ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],
# [0.118 ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],
# [0.119 ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],
# [0.120 ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],
# [0.121 ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],
# [0.122 ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],
# [0.125 ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],
# [0.13  ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],
# [0.14  ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],
# [0.15  ,91.154000,80.370000,1.420000,7.490000,172.500000,1.770000,0.105658],

# [0.000000000001  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.00000000001  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.0000000001  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.000000001  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.00000001  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.0000001  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.000001  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.00001  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.0001  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.001  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.01  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.02  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.03  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.04  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.05  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.06  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.07  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.08  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.09  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.1   ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.11  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.115 ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.116 ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.117 ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.118 ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.119 ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.120 ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.121 ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.122 ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.125 ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.13  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.14  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.15  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],

# [0.16  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.17  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.18  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],
# [0.19  ,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],



#V,MZ,MW,MC,MB,MT,MTAU,MMUON,GF, GAMW,GAMZ
# [0.05,4.5577,4.0185,0.082,0.2435,8.855,0.0885,0.005282918,0.0046652,0.104428,0.1247905,0.2075025],
# [0.1,9.1154,8.037,0.164,0.487,17.71,0.177,0.01056584,0.0011663,0.208856,0.249581,0.1765228],
# [0.15,13.6731,12.0555,0.246,0.7305,26.565,0.2655,0.01584876,0.0005183556,0.313284,0.3743715,0.1623447],
# [0.2,18.2308,16.074,0.328,0.974,35.42,0.354,0.02113167,0.000291575,0.417712,0.499162,0.1535919],
# [0.25,22.7885,20.0925,0.41,1.2175,44.275,0.4425,0.02641459,0.000186608,0.52214,0.6239525,0.1474266],
# [0.3,27.3462,24.111,0.492,1.461,53.13,0.531,0.03169751,0.0001295889,0.626568,0.748743,0.1427449],
[0.35,31.9039,28.1295,0.574,1.7045,61.985,0.6195,0.03698043,9.520816e-05,0.730996,0.8735335,0.1390125],
[0.4,36.4616,32.148,0.656,1.948,70.84,0.708,0.04226335,7.289375e-05,0.835424,0.998324,0.1359337],
[0.45,41.0193,36.1665,0.738,2.1915,79.695,0.7965,0.04754627,5.759506e-05,0.939852,1.123115,0.1333289],
[0.5,45.577,40.185,0.82,2.435,88.55,0.885,0.05282918,4.6652e-05,1.04428,1.247905,0.1310821],
[0.55,50.1347,44.2035,0.902,2.6785,97.405,0.9735,0.0581121,3.855537e-05,1.148708,1.372696,0.1291138],
[0.6,54.6924,48.222,0.984,2.922,106.26,1.062,0.06339502,3.239722e-05,1.253136,1.497486,0.1273679],
[0.65,59.2501,52.2405,1.066,3.1655,115.115,1.1505,0.06867794,2.760473e-05,1.357564,1.622277,0.1258029],
# [0.7,63.8078,56.259,1.148,3.409,123.97,1.239,0.07396086,2.380204e-05,1.461992,1.747067,0.1243879],
# [0.75,68.3655,60.2775,1.23,3.6525,132.825,1.3275,0.07924378,2.073422e-05,1.56642,1.871857,0.1230989],
# [0.8,72.9232,64.296,1.312,3.896,141.68,1.416,0.08452669,1.822344e-05,1.670848,1.996648,0.121917],
# [0.85,77.4809,68.3145,1.394,4.1395,150.535,1.5045,0.08980961,1.614256e-05,1.775276,2.121439,0.1208273],
# [0.9,82.0386,72.333,1.476,4.383,159.39,1.593,0.09509253,1.439877e-05,1.879704,2.246229,0.1198176],
# [0.95,86.5963,76.3515,1.558,4.6265,168.245,1.6815,0.1003754,1.292299e-05,1.984132,2.37102,0.118878],
[1,91.154,80.37,1.64,4.87,177.1,1.77,0.1056584,1.1663e-05,2.08856,2.49581,0.118],

# [1.05,95.7117,84.3885,1.722,5.1135,185.955,1.8585,0.1109413,1.057868e-05,2.192988,2.620601,0.1171768],
# [1.1,100.2694,88.407,1.804,5.357,194.81,1.947,0.1162242,9.638843e-06,2.297416,2.745391,0.1164026],
# [1.15,104.8271,92.4255,1.886,5.6005,203.665,2.0355,0.1215071,8.818904e-06,2.401844,2.870182,0.1156723],
# [1.2,109.3848,96.444,1.968,5.844,212.52,2.124,0.12679,8.099306e-06,2.506272,2.994972,0.1149816],
# [1.25,113.9425,100.4625,2.05,6.0875,221.375,2.2125,0.132073,7.46432e-06,2.6107,3.119763,0.1143268],
# [1.3,118.5002,104.481,2.132,6.331,230.23,2.301,0.1373559,6.901183e-06,2.715128,3.244553,0.1137047],
# [1.35,123.0579,108.4995,2.214,6.5745,239.085,2.3895,0.1426388,6.399451e-06,2.819556,3.369344,0.1131124],
# [1.4,127.6156,112.518,2.296,6.818,247.94,2.478,0.1479217,5.95051e-06,2.923984,3.494134,0.1125475],
# [1.45,132.1733,116.5365,2.378,7.0615,256.795,2.5665,0.1532046,5.547206e-06,3.028412,3.618925,0.1120078],
# [1.5,136.731,120.555,2.46,7.305,265.65,2.655,0.1584876,5.183556e-06,3.13284,3.743715,0.1114912],
# [1.55,141.2887,124.5735,2.542,7.5485,274.505,2.7435,0.1637705,4.854527e-06,3.237268,3.868506,0.110996],
# [1.6,145.8464,128.592,2.624,7.792,283.36,2.832,0.1690534,4.555859e-06,3.341696,3.993296,0.1105208],
# [1.65,150.4041,132.6105,2.706,8.0355,292.215,2.9205,0.1743363,4.28393e-06,3.446124,4.118087,0.1100641],
# [1.7,154.9618,136.629,2.788,8.279,301.07,3.009,0.1796192,4.03564e-06,3.550552,4.242877,0.1096246],
# [1.75,159.5195,140.6475,2.87,8.5225,309.925,3.0975,0.1849021,3.808327e-06,3.65498,4.367668,0.1092012],
# [1.8,164.0772,144.666,2.952,8.766,318.78,3.186,0.1901851,3.599691e-06,3.759408,4.492458,0.1087928],
# [1.85,168.6349,148.6845,3.034,9.0095,327.635,3.2745,0.195468,3.407743e-06,3.863836,4.617249,0.1083985],
# [1.9,173.1926,152.703,3.116,9.253,336.49,3.363,0.2007509,3.230748e-06,3.968264,4.742039,0.1080175],
# [1.95,177.7503,156.7215,3.198,9.4965,345.345,3.4515,0.2060338,3.067193e-06,4.072692,4.86683,0.107649],
# [2,182.308,160.74,3.28,9.74,354.2,3.54,0.2113167,2.91575e-06,4.17712,4.99162,0.1072922],
# [2,182.308,160.74,3.28,9.74,354.2,3.54,0.2113167,2.91575e-06,4.17712,4.99162,0.1072922],
# [2.5,227.885,200.925,4.1,12.175,442.75,4.425,0.2641459,1.86608e-06,5.2214,6.239525,0.1042468],
# [3,273.462,241.11,4.92,14.61,531.3,5.31,0.3169751,1.295889e-06,6.26568,7.48743,0.1018839],
# [3.5,319.039,281.295,5.74,17.045,619.85,6.195,0.3698043,9.520816e-07,7.30996,8.735335,0.09996819],
# [4,364.616,321.48,6.56,19.48,708.4,7.08,0.4226335,7.289375e-07,8.35424,9.98324,0.09836598],
# [4.5,410.193,361.665,7.38,21.915,796.95,7.965,0.4754627,5.759506e-07,9.39852,11.23114,0.09699478],
# [5,455.77,401.85,8.2,24.35,885.5,8.85,0.5282918,4.6652e-07,10.4428,12.47905,0.09580019],
# [5.5,501.347,442.035,9.02,26.785,974.05,9.735,0.581121,3.855537e-07,11.48708,13.72696,0.09474462],


]

# VARNAME = "MB"
VARNAME = "VEV"

BRGraph = {}

BRGraph['scanParam'] = []
BRGraph['BB'] = []
BRGraph['TAU_TAU'] = []
BRGraph['MU_MU'] = []
BRGraph['E_E'] = []
BRGraph['U_U'] = []
BRGraph['D_D'] = []
BRGraph['SS'] = []
BRGraph['CC'] = []
BRGraph['GG'] = []
BRGraph['GAM_GAM'] = []
BRGraph['Z_GAM'] = []
BRGraph['WW'] = []
BRGraph['ZZ'] = []



for (V,MZ,MW,MC,MB,MT,MTAU,MMUON,GF,GAMZ,GAMW,ALPHAS) in VEVSamplePoints:
# for (V,MZ,MW,MC,MB,MT,MTAU,MMUON) in VEVSamplePoints:

	print (V,MZ,MW,MC,MB,MT,MTAU,MMUON)

	# higgsMassStart = "85.D0"
	# higgsMassEnd = "185.D0"
	# nHiggsMassSteps = "6"

	higgsMassStart = "  1.D0"
	higgsMassEnd = "125.D0"
	nHiggsMassSteps = "50"
	
	# if VARNAME=="VEV":
	# 	if float(V) > 700.:
	# 		higgsMassStart = str(int(V/10.))+".D0"
	# 		higgsMassEnd = str(int(9.*V/10.))+".D0"
	# 		nHiggsMassSteps = "101"
	# 	if float(V) < 110.:
	# 		higgsMassStart = str(int(V/3.))+".D0"
	# 		higgsMassEnd = str(int(2.*V/3.))+".D0"
	# 		nHiggsMassSteps = "51"
	# 	if float(V) < 50.:
	# 		higgsMassStart = str(int(3.*V/5.))+".D0"
	# 		higgsMassEnd = str(int(4.*V/5.))+".D0"
	# 		nHiggsMassSteps = "51"
	# if VARNAME=="MFermion":
	# 	if float(V) < 10.:
	# 		higgsMassStart = "80.D0"
	# 		higgsMassEnd = "180.D0"
	# 		nHiggsMassSteps = "51"


	# if float(V) < 10.:
	# 	higgsMassStart = "80.D0"
	# 	higgsMassEnd = "180.D0"
	# 	nHiggsMassSteps = "11"



	# higgsMassStart = "100.D0"
	# higgsMassEnd = "200.D0"
	# nHiggsMassSteps = "4"


	scanParamList = []
	# ALPHAS = V

	for iSample in xrange(1):
		scanParam = {}

		# scanParam['MZ']      = str( Decimal( float(MZ)  *random.gauss(1,0.042/2.085)  )   )+"D0"
		# scanParam['MW']      = str( Decimal( float(MW)  *random.gauss(1,0.042/2.085)  )   )+"D0"
		# scanParam['MC']      = str( Decimal( float(MC)  *random.gauss(1,0.03/1.42)    )   )+"D0"
		# scanParam['MB']      = str( Decimal( float(MB)  *random.gauss(1,0.06/4.49)    )   )+"D0"
		# scanParam['MT']      = str( Decimal( float(MT)  *random.gauss(1,2.5/172.5)    )    )+"D0"
		# scanParam['MTAU']    = str( Decimal( float(MTAU)*random.gauss(1,0.16/1776.82) )    )+"D0"
		# scanParam['MMUON']   = str( Decimal( float(MMUON)   ) )+"D0"
		# # scanParam['GF'] = str( Decimal(GF)  )+"D0"
		# # scanParam['GAMW'] = str( Decimal(GAMW)  )+"D0"
		# # scanParam['GAMZ'] = str( Decimal(GAMZ)  )+"D0"
		# scanParam['ALS(MZ)'] = str( Decimal(ALPHAS*random.gauss(1,0.002/0.119 ) )  )+"D0"

		scanParam['V'] = float(V)
		scanParam['VARNAME'] = VARNAME

		scanParam['MZ']      = str( Decimal( float(MZ)   )   )+"D0"
		scanParam['MW']      = str( Decimal( float(MW)   )   )+"D0"
		scanParam['MC']      = str( Decimal( float(MC)   )   )+"D0"
		scanParam['MB']      = str( Decimal( float(MB)   )   )+"D0"
		scanParam['MT']      = str( Decimal( float(MT)   )    )+"D0"
		scanParam['MTAU']    = str( Decimal( float(MTAU) )    )+"D0"
		scanParam['MMUON']   = str( Decimal( float(MMUON)   ) )+"D0"
		# scanParam['GF'] = str( Decimal(GF)  )+"D0"
		# scanParam['GAMW'] = str( Decimal(GAMW)  )+"D0"
		# scanParam['GAMZ'] = str( Decimal(GAMZ)  )+"D0"
		scanParam['ALS(MZ)'] = '{0:.100f}'.format(ALPHAS)+"D0"


		# print "check format of GF"
		# print scanParam['GF']
		# if VARNAME=="ALS(MZ)":
		# 	scanParam['ALS(MZ)'] = str( float(V)  *  1.  )+"D0"
		# if VARNAME=="1/ALPHA":
		# 	scanParam['1\/ALPHA'] = str( float(V)  *  1.  )+"D0"

		scanParamList.append(scanParam)

	hdecayExe = "HDECAY/run"
	originalTemplateFile = "HDECAY/hdecay.in"
	templateFile = "hdecay.in.sm"

	tmpInputFile = "hdecay.in"

	# Ok let's fuck with the HDECAY input
	# Start by giving me a combo of SM params

	for i,iScan in enumerate(scanParamList):
		print iScan
		shutil.copyfile(templateFile,tmpInputFile)

		os.system("sed -i.bkup 's/^MABEG.*/MABEG%s= %s/' %s "%(" "*(9-len(higgsMassStart)) ,higgsMassStart,tmpInputFile)  )
		os.system("sed -i.bkup 's/^MAEND.*/MAEND%s= %s/' %s "%(" "*(9-len(higgsMassStart)) , higgsMassEnd,tmpInputFile)  )
		os.system("sed -i.bkup 's/^NMA.*/NMA%s= %s/' %s "%(" "*(8-len(nHiggsMassSteps )) ,nHiggsMassSteps,tmpInputFile)  )

		for iVar in iScan:
			if iVar=="V" or iVar=="VARNAME":
				continue
			os.system("sed -i.bkup 's/^%s .*/%s%s= %s/' %s "%(iVar,iVar," "*(9-len(iVar) ),iScan[iVar],tmpInputFile)  )

		os.system("./HDECAY/run")
		try:
			os.mkdir("Point_%s"%'{0:.40f}'.format(V))
		except:
			print "."
		os.system("cp br.* Point_%s/."%'{0:.40f}'.format(V))


		table1File = open("Point_%s/br.sm1"%'{0:.40f}'.format(V),'r')
		for j,line in enumerate(table1File):
		    if j is 1: break
		    varNames = re.split(r'\s\s+',line)[1:-1]

		table1 = np.recfromtxt("Point_%s/br.sm1"%'{0:.40f}'.format(V), skiprows=3, names=varNames )
		table1File.close()

		table2File = open("Point_%s/br.sm2"%'{0:.40f}'.format(V),'r')
		for j,line in enumerate(table2File):
		    if j is 1: break
		    varNames = re.split(r'\s\s+',line)[1:]
		    print varNames

		table2 = np.recfromtxt("Point_%s/br.sm2"%'{0:.40f}'.format(V), skiprows=3, names=varNames )
		table2File.close()



		data = []
		for iHiggsMass in xrange(len(table1) ):
			# print table1
			scaledBRToUU = (3./100.)**2*table1[iHiggsMass]['SS']
			scaledBRToDD = (5./100.)**2*table1[iHiggsMass]['SS']
			scaledBRToEE = (0.511/105.7)**2*table1[iHiggsMass]['MU_MU']

			productOfBRs = \
				table1[iHiggsMass]['BB']* \
				table1[iHiggsMass]['TAU_TAU']* \
				table1[iHiggsMass]['MU_MU']* \
				scaledBRToEE* \
				scaledBRToUU* \
				scaledBRToDD* \
				table1[iHiggsMass]['SS']* \
				table1[iHiggsMass]['CC']* \
				table1[iHiggsMass]['TT']* \
				table2[iHiggsMass]['GG']* \
				table2[iHiggsMass]['GAM_GAM']* \
				table2[iHiggsMass]['Z_GAM']* \
				table2[iHiggsMass]['WW']* \
				table2[iHiggsMass]['ZZ']


			productOfBRsExclTop = \
				table1[iHiggsMass]['BB']* \
				table1[iHiggsMass]['TAU_TAU']* \
				table1[iHiggsMass]['MU_MU']* \
				scaledBRToEE* \
				scaledBRToUU* \
				scaledBRToDD* \
				table1[iHiggsMass]['SS']* \
				table1[iHiggsMass]['CC']* \
				table2[iHiggsMass]['GG']* \
				table2[iHiggsMass]['GAM_GAM']* \
				table2[iHiggsMass]['Z_GAM']* \
				table2[iHiggsMass]['WW']* \
				table2[iHiggsMass]['ZZ']

			print "TEST"
			print table1[iHiggsMass]['BB']
			print table1[iHiggsMass]['TAU_TAU']
			print table1[iHiggsMass]['MU_MU']
			print scaledBRToEE
			print scaledBRToUU
			print scaledBRToDD
			print table1[iHiggsMass]['SS']
			print table1[iHiggsMass]['CC']
			print table2[iHiggsMass]['GG']
			print table2[iHiggsMass]['GAM_GAM']
			print table2[iHiggsMass]['Z_GAM']
			print table2[iHiggsMass]['WW']
			print table2[iHiggsMass]['ZZ']
			print productOfBRsExclTop

			# productOfBRs = abs(productOfBRs)
			# productOfBRsExclTop = abs(productOfBRsExclTop)

			if productOfBRsExclTop<0:
				continue

			# print "PRODUCT IS............."
			# print productOfBRsExclTop
			# print table2[iHiggsMass]['GG']

			if table1[iHiggsMass]['MHSM'] < 2.*MT :
				if productOfBRsExclTop:
					data.append( (table1[iHiggsMass]['MHSM'], productOfBRsExclTop)  )
			else:
				if productOfBRs:
					data.append( (table1[iHiggsMass]['MHSM'], productOfBRs)  )

		dataArray = np.array(data)
		arrayHiggsMass = dataArray[:,0]
		arrayProductOfBRs = dataArray[:,1]
		# arrayProductOfBRsExclTop = dataArray[:,2]

		# print arrayHiggsMass
		# print arrayProductOfBRsExclTop

		graph = Graph(dataArray.shape[0] )
		for iPoint, (xx, yy) in enumerate(zip(arrayHiggsMass,arrayProductOfBRs) ):

			graph.SetPoint(iPoint, xx, yy)
			PDFArray['MH'].append(xx)
			PDFArray['V'].append(V)
			PDFArray['P'].append(yy)


		tmpintegral = graph.Integral()
		tmpintegral = 0
		for ipoint in xrange(graph.GetN()):
			xx,yy = Double(),Double()
			graph.GetPoint(ipoint,xx,yy)
			tmpintegral+=yy

		print "Integral is: %.128f"%tmpintegral
		if tmpintegral!=tmpintegral:
			continue



		# canvas = Canvas()
		# graph.Draw("APL")
		# print arrayHiggsMass
		# fitFcn = TF1("fitFcn","[0]*exp(-(x-[1])**2/(2*[2]**2))",0,500.)
		# fitFcn = TF1("fitFcn","gaus",arrayHiggsMass[0],arrayHiggsMass[-1])
		# fitFcn.SetParameters(1.0e-24,V/2.,V)
		# fitFcn.SetParLimits(0, 0., 1.e-5);
		# fitFcn.SetParLimits(1, 
		# 	float(arrayHiggsMass[0]), 
		# 	float(arrayHiggsMass[-1])  );
		# fitFcn.SetParLimits(2, V*0.01, V*10.);


		graph.Fit("gaus","")

		fit = graph.GetFunction("gaus")
		fitConstant = fit.GetParameter(0)
		fitConstantErr = fit.GetParError(0)
		fitMean = fit.GetParameter(1)
		fitMeanErr = fit.GetParError(1)
		fitSigma = fit.GetParameter(2)
		fitSigmaErr = fit.GetParError(2)

		print fit.GetMaximum()

		MHhatArray['MHhat'].append(fitMean)
		MHhatArray['P'].append(fit.GetMaximum() )
		MHhatArray['V'].append(V)


		print "Sanity Check me.........."
		print V
		print table1[0]['BB']

		# BRGraph['scanParam'].append(V)
		# BRGraph['BB'].append(table1[4]['BB'] )
		# BRGraph['TAU_TAU'].append(table1[4]['TAU_TAU'] )
		# BRGraph['MU_MU'].append(table1[4]['MU_MU'] )
		# BRGraph['E_E'].append((0.511/105.7)**2*table1[4]['MU_MU']  )
		# BRGraph['U_U'].append((3./100.)**2*table1[4]['SS'] )
		# BRGraph['D_D'].append((5./100.)**2*table1[4]['SS'] )
		# BRGraph['SS'].append(table1[4]['SS'] )
		# BRGraph['CC'].append(table1[4]['CC'] )
		# BRGraph['GG'].append(table2[4]['GG'] )
		# BRGraph['GAM_GAM'].append(table2[4]['GAM_GAM'] )
		# BRGraph['Z_GAM'].append(table2[4]['Z_GAM'] )
		# BRGraph['WW'].append(table2[4]['WW'] )
		# BRGraph['ZZ'].append(table2[4]['ZZ'] )

		BRGraph['scanParam'].append(V)
		BRGraph['BB'].append(table1[0]['BB'] )
		BRGraph['TAU_TAU'].append(table1[0]['TAU_TAU'] )
		BRGraph['MU_MU'].append(table1[0]['MU_MU'] )
		BRGraph['E_E'].append((0.511/105.7)**2*table1[0]['MU_MU']  )
		BRGraph['U_U'].append((3./100.)**2*table1[0]['SS'] )
		BRGraph['D_D'].append((5./100.)**2*table1[0]['SS'] )
		BRGraph['SS'].append(table1[0]['SS'] )
		BRGraph['CC'].append(table1[0]['CC'] )
		BRGraph['GG'].append(table2[0]['GG'] )
		BRGraph['GAM_GAM'].append(table2[0]['GAM_GAM'] )
		BRGraph['Z_GAM'].append(table2[0]['Z_GAM'] )
		BRGraph['WW'].append(table2[0]['WW'] )
		BRGraph['ZZ'].append(table2[0]['ZZ'] )

		fig,axes = plt.subplots()
		axes.plot(arrayHiggsMass,arrayProductOfBRs,'o', markeredgewidth=0)
		tmpx = np.linspace(float(higgsMassStart[:-3]),
			float(higgsMassEnd[:-3]),
			float(nHiggsMassSteps[:])*10.)
		# tmpx = np.linspace(120,140,100)
		axes.plot(tmpx,gauss_function(tmpx, fitConstant, fitMean, fitSigma) )

		axes.set_yscale('log')
		axes.set_xlabel(r"$M_{H}$ [GeV]")
		axes.set_ylabel(r"$\prod_i BR_i$")
		axes.annotate(r"Mean: $%.2f \pm %.2f$"%(fitMean,fitMeanErr) , xy=(0.05, 0.95), xycoords='axes fraction', fontweight='bold', fontsize=8)
		axes.annotate(r"Max  : $%.50f$"%(fit.GetMaximum() ) , xy=(0.05, 0.9), xycoords='axes fraction', fontweight='bold', fontsize=8)
		axes.annotate(r"@125: $%.50f$"%(fit.Eval(125.) ) , xy=(0.05, 0.87), xycoords='axes fraction', fontweight='bold', fontsize=8)
		plt.grid()
		plt.savefig("Point_%s.pdf"%'{0:.40f}'.format(V))

		# MHhatValues.append(fitMean)

f = open("PDFArray.txt",'w')
f2 = open("MHhatArray.txt",'w')
outputfile = TFile("output.root","RECREATE")

canvas = Canvas()
graph2D = Graph2D(len(PDFArray['V'])) 
for i in xrange(len(PDFArray['V']) ):
	graph2D.SetPoint(i,PDFArray['V'][i],PDFArray['MH'][i],PDFArray['P'][i])
	f.write('%f %f %.128f\n'%(PDFArray['V'][i],PDFArray['MH'][i],PDFArray['P'][i]) )

f.close()

graph1D = Graph(len(MHhatArray['V'])) 
for i in xrange(len(MHhatArray['V']) ):
	graph1D.SetPoint(i,MHhatArray["V"][i], MHhatArray["P"][i] )
	f2.write('%f %f \n'%(MHhatArray['V'][i],MHhatArray['MHhat'][i] ) )
f2.close()


##############################################

# modes = ['BB',
# 		'TAU_TAU',
# 		'MU_MU',
# 		'E_E',
# 		'U_U',
# 		'D_D',
# 		'SS',
# 		'CC',
# 		'GG',
# 		'GAM_GAM',
# 		'Z_GAM',
# 		'WW',
# 		'ZZ']

# BRGraphPlot = {}
# for mode in modes:
# 	BRGraphPlot[mode] = Graph(len(BRGraph['scanParam'] )) 

# print "@@@@@@@@@@@@@@@@@@@@@@@"
# print BRGraph["BB"]
# print "@@@@@@@@@@@@@@@@@@@"

# for i in xrange(len(BRGraph['scanParam'] ) ):
# 	for mode in modes:
# 		print mode
# 		print BRGraph['scanParam'][i]
# 		print BRGraph[mode][i]
# 		BRGraphPlot[mode].SetPoint(i,float(BRGraph["scanParam"][i]), float(BRGraph[mode][i]) )
# #	f2.write('%f %f \n'%(MHhatArray['V'][i],MHhatArray['MHhat'][i] ) )
# #f2.close()

# icolor = 1
# imarker = 22

# drawopts = "ALP"

# print BRGraph['scanParam']
# realAlphaSIndex = BRGraph['scanParam'].index(0.118)

# for i,mode in enumerate(modes):

# 	if icolor>5:
# 		icolor=1
# 		imarker=imarker+1

# 	x, y = Double(0),Double(0)
# 	BRGraphPlot[mode].GetPoint(realAlphaSIndex,x,y)
# 	centralValue = y 

# 	for ipoint in xrange(BRGraphPlot[mode].GetN() ):
# 		x, y = Double(0),Double(0)
# 		BRGraphPlot[mode].GetPoint(ipoint,x,y)
# 		BRGraphPlot[mode].SetPoint(ipoint,x,y/centralValue)


# 	BRGraphPlot[mode].SetLineColor(icolor)
# 	BRGraphPlot[mode].SetMarkerColor(icolor)
# 	BRGraphPlot[mode].SetMarkerStyle(imarker)
# 	icolor = icolor+1

# 	BRGraphPlot[mode].SetTitle(mode)
# 	BRGraphPlot[mode].Draw(drawopts)
# 	drawopts = "LP"

# # canvas.SetLogy()
# canvas.SetLogx()
# canvas.SetGrid()

# # BRGraphPlot["BB"].SetMinimum(0.000000000000001)
# # BRGraphPlot["BB"].SetMaximum(1)


# BRGraphPlot["BB"].SetMinimum(0)
# BRGraphPlot["BB"].SetMaximum(2)


# # BRGraphPlot["BB"].GetXaxis().SetRangeUser(0.0000000,.5)
# BRGraphPlot["BB"].GetXaxis().SetTitle("#alpha_{s}")
# BRGraphPlot["BB"].GetYaxis().SetTitle("BR/BR(#hat{#alpha_{s}})")


# canvas.BuildLegend(0.55,0.2,0.68,0.45)

# label = TLatex(0.06,10e-14,"M_{H} = 125 GeV")
# label = TLatex(0.06,0.05,"M_{H} = 125 GeV")
# label.Draw()

# canvas.SaveAs("BRvsScanParam_125.eps")

########################################################

#graph2D.Draw("colz")
graph2D_th2 = graph2D.GetHistogram()
graph2D_th2.Draw("colz")
graph2D_th2_pfx = graph2D_th2.ProfileX()
graph2D_th2_pfx.Draw("sames")

graph2D_th2_pfy = graph2D_th2.ProfileY()
graph2D_th2_pfy.Draw("sames")

gPad.SetLogz()
canvas.SaveAs("holyshit.eps")

graph1D.Draw("ALP")
graph1D.Write()
canvas.SaveAs("P_vs_V.eps")

graph2D.Write()
graph2D_th2.Write()
outputfile.Write()
outputfile.Close()

# c1 = Canvas()
# a = Hist(200, 0, 2000)
# a.fill_array(MHhatValues)
# a.Draw('hist')
# rplt.bar(a)
# plt.savefig("PMhhat.pdf")

