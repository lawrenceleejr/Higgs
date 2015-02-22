#!/usr/bin/env python

import sys,os
import numpy as np
import shutil
import re
import random
import math

from ROOT import *
gROOT.LoadMacro("atlasstyle/AtlasStyle.C")
SetAtlasStyle()

# from rootpy.interactive import wait
from rootpy.plotting import Canvas, Hist, Hist2D, Hist3D, Graph, Graph2D


import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from decimal import Decimal

from operator import mul


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


	
SamplePoints = [

#### Nominal 
# [246.000000,91.154000,80.370000,1.420000,4.490000,172.500000,1.770000,0.105658],


# #### ALS(MZ)
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



#V,MZ,MW,MC,MB,MT,MTAU,MMUON,GF, GAMW,GAMZ, ALPHAS


# [0.15,13.6731,12.0555,0.246,0.7305,26.565,0.2664,0.01584876,0.0005183556,0.313284,0.3743715,0.1623447],

[0.2,18.2308,16.074,0.328,0.974,35.42,0.3552,0.02113167,0.000291575,0.417712,0.499162,0.1535919],
# [0.25,22.7885,20.0925,0.41,1.2175,44.275,0.444,0.02641459,0.000186608,0.52214,0.6239525,0.1474266],
[0.3,27.3462,24.111,0.492,1.461,53.13,0.5328,0.03169751,0.0001295889,0.626568,0.748743,0.1427449],
# [0.35,31.9039,28.1295,0.574,1.7045,61.985,0.6216,0.03698043,9.520816e-05,0.730996,0.8735335,0.1390125],
[0.4,36.4616,32.148,0.656,1.948,70.84,0.7104,0.04226335,7.289375e-05,0.835424,0.998324,0.1359337],
# [0.45,41.0193,36.1665,0.738,2.1915,79.695,0.7992,0.04754627,5.759506e-05,0.939852,1.123115,0.1333289],
[0.5,45.577,40.185,0.82,2.435,88.55,0.888,0.05282918,4.6652e-05,1.04428,1.247905,0.1310821],
# [0.55,50.1347,44.2035,0.902,2.6785,97.405,0.9768,0.0581121,3.855537e-05,1.148708,1.372696,0.1291138],
[0.6,54.6924,48.222,0.984,2.922,106.26,1.0656,0.06339502,3.239722e-05,1.253136,1.497486,0.1273679],
# [0.65,59.2501,52.2405,1.066,3.1655,115.115,1.1544,0.06867794,2.760473e-05,1.357564,1.622277,0.1258029],
[0.7,63.8078,56.259,1.148,3.409,123.97,1.2432,0.07396086,2.380204e-05,1.461992,1.747067,0.1243879],
# [0.75,68.3655,60.2775,1.23,3.6525,132.825,1.332,0.07924378,2.073422e-05,1.56642,1.871857,0.1230989],
[0.8,72.9232,64.296,1.312,3.896,141.68,1.4208,0.08452669,1.822344e-05,1.670848,1.996648,0.121917],
# [0.85,77.4809,68.3145,1.394,4.1395,150.535,1.5096,0.08980961,1.614256e-05,1.775276,2.121439,0.1208273],
[0.9,82.0386,72.333,1.476,4.383,159.39,1.5984,0.09509253,1.439877e-05,1.879704,2.246229,0.1198176],
# [0.95,86.5963,76.3515,1.558,4.6265,168.245,1.6872,0.1003754,1.292299e-05,1.984132,2.37102,0.118878],
[1,91.154,80.37,1.64,4.87,177.1,1.776,0.1056584,1.1663e-05,2.08856,2.49581,0.118],
# [1.05,95.7117,84.3885,1.722,5.1135,185.955,1.8648,0.1109413,1.057868e-05,2.192988,2.620601,0.1171768],
[1.1,100.2694,88.407,1.804,5.357,194.81,1.9536,0.1162242,9.638843e-06,2.297416,2.745391,0.1164026],
# [1.15,104.8271,92.4255,1.886,5.6005,203.665,2.0424,0.1215071,8.818904e-06,2.401844,2.870182,0.1156723],
[1.2,109.3848,96.444,1.968,5.844,212.52,2.1312,0.12679,8.099306e-06,2.506272,2.994972,0.1149816],
# [1.25,113.9425,100.4625,2.05,6.0875,221.375,2.22,0.132073,7.46432e-06,2.6107,3.119763,0.1143268],
[1.3,118.5002,104.481,2.132,6.331,230.23,2.3088,0.1373559,6.901183e-06,2.715128,3.244553,0.1137047],
# [1.35,123.0579,108.4995,2.214,6.5745,239.085,2.3976,0.1426388,6.399451e-06,2.819556,3.369344,0.1131124],
[1.4,127.6156,112.518,2.296,6.818,247.94,2.4864,0.1479217,5.95051e-06,2.923984,3.494134,0.1125475],
# [1.45,132.1733,116.5365,2.378,7.0615,256.795,2.5752,0.1532046,5.547206e-06,3.028412,3.618925,0.1120078],
[1.5,136.731,120.555,2.46,7.305,265.65,2.664,0.1584876,5.183556e-06,3.13284,3.743715,0.1114912],
# [1.55,141.2887,124.5735,2.542,7.5485,274.505,2.7528,0.1637705,4.854527e-06,3.237268,3.868506,0.110996],
[1.6,145.8464,128.592,2.624,7.792,283.36,2.8416,0.1690534,4.555859e-06,3.341696,3.993296,0.1105208],
[1.65,150.4041,132.6105,2.706,8.0355,292.215,2.9304,0.1743363,4.28393e-06,3.446124,4.118087,0.1100641],
[1.7,154.9618,136.629,2.788,8.279,301.07,3.0192,0.1796192,4.03564e-06,3.550552,4.242877,0.1096246],
[1.75,159.5195,140.6475,2.87,8.5225,309.925,3.108,0.1849021,3.808327e-06,3.65498,4.367668,0.1092012],
[1.8,164.0772,144.666,2.952,8.766,318.78,3.1968,0.1901851,3.599691e-06,3.759408,4.492458,0.1087928],
[1.85,168.6349,148.6845,3.034,9.0095,327.635,3.2856,0.195468,3.407743e-06,3.863836,4.617249,0.1083985],
[1.9,173.1926,152.703,3.116,9.253,336.49,3.3744,0.2007509,3.230748e-06,3.968264,4.742039,0.1080175],
[1.95,177.7503,156.7215,3.198,9.4965,345.345,3.4632,0.2060338,3.067193e-06,4.072692,4.86683,0.107649],
[2,182.308,160.74,3.28,9.74,354.2,3.552,0.2113167,2.91575e-06,4.17712,4.99162,0.1072922],
# [2.05,186.8657,164.7585,3.362,9.9835,363.055,3.6408,0.2165997,2.775253e-06,4.281548,5.116411,0.1069464],


# [2.1,191.4234,168.777,3.444,10.227,371.91,3.7296,0.2218826,2.644671e-06,4.385976,5.241201,0.1066112],
# [2.15,195.9811,172.7955,3.526,10.4705,380.765,3.8184,0.2271655,2.523094e-06,4.490404,5.365991,0.1062858],
# [2.2,200.5388,176.814,3.608,10.714,389.62,3.9072,0.2324484,2.409711e-06,4.594832,5.490782,0.1059699],
# [2.25,205.0965,180.8325,3.69,10.9575,398.475,3.996,0.2377313,2.303802e-06,4.69926,5.615572,0.1056629],
# [2.3,209.6542,184.851,3.772,11.201,407.33,4.0848,0.2430142,2.204726e-06,4.803688,5.740363,0.1053643],
# [2.35,214.2119,188.8695,3.854,11.4445,416.185,4.1736,0.2482972,2.111906e-06,4.908116,5.865154,0.1050738],
# [2.4,218.7696,192.888,3.936,11.688,425.04,4.2624,0.2535801,2.024826e-06,5.012544,5.989944,0.1047909],
# [2.45,223.3273,196.9065,4.018,11.9315,433.895,4.3512,0.258863,1.943024e-06,5.116972,6.114735,0.1045154],
[2.5,227.885,200.925,4.1,12.175,442.75,4.44,0.2641459,1.86608e-06,5.2214,6.239525,0.1042468],
# [2.55,232.4427,204.9435,4.182,12.4185,451.605,4.5288,0.2694288,1.793618e-06,5.325828,6.364316,0.1039848],
# [2.6,237.0004,208.962,4.264,12.662,460.46,4.6176,0.2747118,1.725296e-06,5.430256,6.489106,0.1037293],
# [2.65,241.5581,212.9805,4.346,12.9055,469.315,4.7064,0.2799947,1.660805e-06,5.534684,6.613897,0.1034798],
# [2.7,246.1158,216.999,4.428,13.149,478.17,4.7952,0.2852776,1.599863e-06,5.639112,6.738687,0.1032361],
# [2.75,250.6735,221.0175,4.51,13.3925,487.025,4.884,0.2905605,1.542215e-06,5.74354,6.863478,0.1029981],
# [2.8,255.2312,225.036,4.592,13.636,495.88,4.9728,0.2958434,1.487628e-06,5.847968,6.988268,0.1027654],
# [2.85,259.7889,229.0545,4.674,13.8795,504.735,5.0616,0.3011263,1.435888e-06,5.952396,7.113059,0.1025378],
# [2.9,264.3466,233.073,4.756,14.123,513.59,5.1504,0.3064093,1.386801e-06,6.056824,7.237849,0.1023152],
# [2.95,268.9043,237.0915,4.838,14.3665,522.445,5.2392,0.3116922,1.34019e-06,6.161252,7.36264,0.1020973],
[3,273.462,241.11,4.92,14.61,531.3,5.328,0.3169751,1.295889e-06,6.26568,7.48743,0.1018839],
# [3.05,278.0197,245.1285,5.002,14.8535,540.155,5.4168,0.322258,1.253749e-06,6.370108,7.612221,0.101675],
# [3.1,282.5774,249.147,5.084,15.097,549.01,5.5056,0.3275409,1.213632e-06,6.474536,7.737011,0.1014703],
# [3.15,287.1351,253.1655,5.166,15.3405,557.865,5.5944,0.3328239,1.175409e-06,6.578964,7.861802,0.1012697],
# [3.2,291.6928,257.184,5.248,15.584,566.72,5.6832,0.3381068,1.138965e-06,6.683392,7.986592,0.101073],
# [3.25,296.2505,261.2025,5.33,15.8275,575.575,5.772,0.3433897,1.104189e-06,6.78782,8.111382,0.1008801],
# [3.3,300.8082,265.221,5.412,16.071,584.43,5.8608,0.3486726,1.070983e-06,6.892248,8.236173,0.1006909],
# [3.35,305.3659,269.2395,5.494,16.3145,593.285,5.9496,0.3539555,1.039252e-06,6.996676,8.360964,0.1005052],
# [3.4,309.9236,273.258,5.576,16.558,602.14,6.0384,0.3592384,1.00891e-06,7.101104,8.485754,0.1003229],
# [3.45,314.4813,277.2765,5.658,16.8015,610.995,6.1272,0.3645214,9.798782e-07,7.205532,8.610545,0.100144],
[3.5,319.039,281.295,5.74,17.045,619.85,6.216,0.3698043,9.520816e-07,7.30996,8.735335,0.09996819],
# [3.55,323.5967,285.3135,5.822,17.2885,628.705,6.3048,0.3750872,9.254513e-07,7.414388,8.860126,0.09979551],
# [3.6,328.1544,289.332,5.904,17.532,637.56,6.3936,0.3803701,8.999228e-07,7.518816,8.984916,0.09962584],
# [3.65,332.7121,293.3505,5.986,17.7755,646.415,6.4824,0.385653,8.754363e-07,7.623244,9.109707,0.09945908],
# [3.7,337.2698,297.369,6.068,18.019,655.27,6.5712,0.390936,8.519357e-07,7.727672,9.234497,0.09929512],
# [3.75,341.8275,301.3875,6.15,18.2625,664.125,6.66,0.3962189,8.293689e-07,7.8321,9.359288,0.0991339],
# [3.8,346.3852,305.406,6.232,18.506,672.98,6.7488,0.4015018,8.07687e-07,7.936528,9.484078,0.09897532],
# [3.85,350.9429,309.4245,6.314,18.7495,681.835,6.8376,0.4067847,7.868443e-07,8.040956,9.608868,0.09881932],
# [3.9,355.5006,313.443,6.396,18.993,690.69,6.9264,0.4120676,7.667982e-07,8.145384,9.733659,0.09866581],
# [3.95,360.0583,317.4615,6.478,19.2365,699.545,7.0152,0.4173505,7.475084e-07,8.249812,9.85845,0.09851472],
# [4,364.616,321.48,6.56,19.48,708.4,7.104,0.4226335,7.289375e-07,8.35424,9.98324,0.09836598],




]

# VARNAME = "MB"

ParameterNameArray = ['V','MZ','MW','MC','MB','MT','MTAU','MMUON','GF', 'GAMW','GAMZ', 'ALS(MZ)']
VARNAME = "V"

doErrors=0
FractionalErrors = {}
FractionalErrors["MZ"] = 0.042/2.085
FractionalErrors["MW"] = 0.042/2.085 
FractionalErrors["MC"] = 0.03/1.42   
FractionalErrors["MB"] = 0.06/4.49
FractionalErrors["MT"] = 2.5/172.5   
FractionalErrors["MTAU"] = 0.16/1776.82
FractionalErrors["ALPHAS"] = 0.002/0.119



modes = [
		'TAU_TAU',
		'MU_MU',
		'EE',
		'BB',
		'CC',
		'SS',
		'DD',
		'UU',
		'GG',
		'GAM_GAM',
		'Z_GAM',
		'WW',
		'ZZ']

BRGraph = {}
for mode in modes:
	BRGraph[mode] = []
BRGraph['scanParam'] = []


for SamplePoint in SamplePoints:


	ParameterDictionary = dict(zip(ParameterNameArray, SamplePoint))

	print ParameterDictionary

	endvalue = 2*125*ParameterDictionary['V']
	if endvalue<60:
		endvalue=60

	higgsMassStart = "  1.D0"
	higgsMassEnd = '{0:.10f}'.format(endvalue) + "D0" #"145.D0"
	nHiggsMassSteps = "20"
	

	scanParamList = []
	# ALPHAS = V

	for iSample in xrange(1):
		scanParam = {}

		for iParameter in ParameterNameArray:
			if doErrors and (iParameter in FractionalErrors):
				scanParam[iParameter] = '{0:.100f}'.format(  ParameterDictionary[iParameter] * random.gauss(1,FractionalErrors[iParameter])  )+"D0"
			else:
				scanParam[iParameter] = '{0:.100f}'.format(  ParameterDictionary[iParameter]  )+"D0"

		scanParamList.append(scanParam)

	hdecayExe = "HDECAY/run"
	templateFile = "hdecay.in.sm"

	tmpInputFile = "hdecay.in"

	# Ok let's fuck with the HDECAY input
	# Start by giving me a combo of SM params

	for i,iScan in enumerate(scanParamList):
		# print iScan
		shutil.copyfile(templateFile,tmpInputFile)

		os.system("sed -i.bkup 's/^MABEG.*/MABEG%s= %s/' %s "%(" "*(9-len(higgsMassStart)) ,higgsMassStart,tmpInputFile)  )
		os.system("sed -i.bkup 's/^MAEND.*/MAEND%s= %s/' %s "%(" "*(9-len(higgsMassStart)) , higgsMassEnd,tmpInputFile)  )
		os.system("sed -i.bkup 's/^NMA.*/NMA%s= %s/' %s "%(" "*(8-len(nHiggsMassSteps )) ,nHiggsMassSteps,tmpInputFile)  )

		for iVar in iScan:
			os.system("sed -i.bkup 's/^%s .*/%s%s= %s/' %s "%(iVar,iVar," "*(9-len(iVar) ),iScan[iVar],tmpInputFile)  )

		os.system("./HDECAY/run")
		directoryName = "Point_%s_%s"%( VARNAME, '{0:.40f}'.format(ParameterDictionary[VARNAME]) )
		try:
			os.mkdir(directoryName)
		except:
			print "."
		os.system("cp br.* %s/."%directoryName)


		table1File = open("%s/br.sm1"%(directoryName), 'r')
		for j,line in enumerate(table1File):
		    if j is 1: break
		    varNames = re.split(r'\s\s+',line)[1:-1]

		table1 = np.recfromtxt("%s/br.sm1"%directoryName  ,  skiprows=3, names=varNames )
		table1File.close()

		# print table1[0]

		table2File = open("%s/br.sm2"%directoryName    ,'r')
		for j,line in enumerate(table2File):
		    if j is 1: break
		    varNames = re.split(r'\s\s+',line)[1:]
		    print varNames

		table2 = np.recfromtxt("%s/br.sm2"%directoryName, skiprows=3, names=varNames )
		table2File.close()

		data = []

		tmpBRs = {}

		BRLibrary = {}

		for iHiggsMass in xrange(len(table1) ):
			# print table1

			BRLibrary[iHiggsMass] = {}

			BRLibrary[iHiggsMass]['MHSM'] = table1[iHiggsMass]['MHSM']

			# Scale to s
			# BRLibrary[iHiggsMass]['UU'] = (3./100.)**2*table1[iHiggsMass]['SS']
			# BRLibrary[iHiggsMass]['DD'] = (5./100.)**2*table1[iHiggsMass]['SS']

			# Scale to c
			# BRLibrary[iHiggsMass]['UU'] = (3./1290.)**2*table1[iHiggsMass]['CC']
			# BRLibrary[iHiggsMass]['DD'] = (5./1290.)**2*table1[iHiggsMass]['CC']
			# BRLibrary[iHiggsMass]['SS'] = (100./1290.)**2*table1[iHiggsMass]['CC']

			# Scale to b
			BRLibrary[iHiggsMass]['UU'] = (3./4190.)**2*table1[iHiggsMass]['BB']
			BRLibrary[iHiggsMass]['DD'] = (5./4190.)**2*table1[iHiggsMass]['BB']
			BRLibrary[iHiggsMass]['SS'] = (100./4190.)**2*table1[iHiggsMass]['BB']
			# BRLibrary[iHiggsMass]['CC'] = (1290./4190.)**2*table1[iHiggsMass]['BB']

			BRLibrary[iHiggsMass]['EE'] = (0.511/105.7)**2*table1[iHiggsMass]['MU_MU']

			BRLibrary[iHiggsMass]['TT']      = table1[iHiggsMass]['TT'] 
			BRLibrary[iHiggsMass]['BB']      = table1[iHiggsMass]['BB'] 
			BRLibrary[iHiggsMass]['TAU_TAU'] = table1[iHiggsMass]['TAU_TAU']      
			BRLibrary[iHiggsMass]['MU_MU']   = table1[iHiggsMass]['MU_MU']   
			# BRLibrary[iHiggsMass]['SS']      = table1[iHiggsMass]['SS']    
			BRLibrary[iHiggsMass]['CC']      = table1[iHiggsMass]['CC']    
			BRLibrary[iHiggsMass]['GG']      = table2[iHiggsMass]['GG']    
			BRLibrary[iHiggsMass]['GAM_GAM'] = table2[iHiggsMass]['GAM_GAM']         
			BRLibrary[iHiggsMass]['Z_GAM']   = table2[iHiggsMass]['Z_GAM']       
			BRLibrary[iHiggsMass]['WW']      = table2[iHiggsMass]['WW']    
			BRLibrary[iHiggsMass]['ZZ']      = table2[iHiggsMass]['ZZ']    


			productOfBRs = 1
			for mode in modes:
				# print '%s      %f'%(imode, BRLibrary[iHiggsMass][imode] )
				if BRLibrary[iHiggsMass][mode]<0:
					# print "SKIPPING---------"
					# continue
					BRLibrary[iHiggsMass][mode]=0
				# assert (0<=BRLibrary[iHiggsMass][imode]<=1)
				productOfBRs *= BRLibrary[iHiggsMass][mode]

			productOfBRsExclTop = productOfBRs
			productOfBRs = productOfBRs*table1[iHiggsMass]["TT"]

			if math.isnan(productOfBRsExclTop) or math.isnan(productOfBRs):
				productOfBRs = 0
				productOfBRsExclTop = 0


			print "PRODUCT IS............. @ %f"%BRLibrary[iHiggsMass]['MHSM'] 
			print productOfBRsExclTop

			if BRLibrary[iHiggsMass]['MHSM'] < 2.*ParameterDictionary["MT"] :
				data.append( (BRLibrary[iHiggsMass]['MHSM'], productOfBRsExclTop)  )
			else:
				data.append( (BRLibrary[iHiggsMass]['MHSM'], productOfBRs)  )

		dataArray = np.array(data)
		arrayHiggsMass = dataArray[:,0]
		arrayProductOfBRs = dataArray[:,1]

		graph = Graph(dataArray.shape[0] )
		for iPoint, (xx, yy) in enumerate(zip(arrayHiggsMass,arrayProductOfBRs) ):

			graph.SetPoint(iPoint, xx, yy)
			PDFArray['MH'].append(xx)
			PDFArray['V'].append(ParameterDictionary[VARNAME])
			PDFArray['P'].append(yy)

		# print arrayProductOfBRs
		tmpintegral = sum(arrayProductOfBRs)

		print "Integral is: %.128f"%tmpintegral
		if tmpintegral!=tmpintegral:
			continue


		graph.Fit("gaus","")

		fit = graph.GetFunction("gaus")
		fitConstant = fit.GetParameter(0)
		fitConstantErr = fit.GetParError(0)
		fitMean = fit.GetParameter(1)
		fitMeanErr = fit.GetParError(1)
		fitSigma = fit.GetParameter(2)
		fitSigmaErr = fit.GetParError(2)

		MHhatArray['MHhat'].append(fitMean)
		MHhatArray['P'].append(fit.GetMaximum() )
		MHhatArray['V'].append(ParameterDictionary[VARNAME])

		tmpMean = graph.GetMean()
		print tmpMean
		maxMHBin = -1
		tmplistofMH = [x[0] for x in table1]
		for i,tmpMH in enumerate(tmplistofMH):
			if tmpMH < tmpMean < tmplistofMH[i+1]:
				maxMHBin = i
				break
		if maxMHBin == -1:
			print "YOU'RE FUCKED. CHECK THE MAX MH"

		BRGraph['scanParam'].append(ParameterDictionary[VARNAME] )

		for mode in modes:
			BRGraph[mode].append(BRLibrary[maxMHBin][mode]  )

		# print BRGraph

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
		plt.savefig("%s.pdf"%directoryName)


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


##############################################

BRGraphPlot = {}
for mode in modes:
	BRGraphPlot[mode] = Graph(len(BRGraph['scanParam'] )) 

# print "@@@@@@@@@@@@@@@@@@@@@@@"
# print BRGraph["BB"]
# print "@@@@@@@@@@@@@@@@@@@"

for i in xrange(len(BRGraph['scanParam'] ) ):
	for mode in modes:
		# print mode
		# print BRGraph['scanParam'][i]
		# print BRGraph[mode][i]
		BRGraphPlot[mode].SetPoint(i,float(BRGraph["scanParam"][i]), float(BRGraph[mode][i]) )
#	f2.write('%f %f \n'%(MHhatArray['V'][i],MHhatArray['MHhat'][i] ) )
#f2.close()

icolor = 1
imarker = 20
markersize = 0.8

drawopts = "ALP"

print BRGraph['scanParam']
normValueIndex = BRGraph['scanParam'].index(1)

for i,mode in enumerate(modes):

	if icolor>4:
		icolor=1
		imarker=imarker+1
		markersize = 0.8

	x, y = Double(0),Double(0)
	BRGraphPlot[mode].GetPoint(normValueIndex,x,y)
	centralValue = y 

	for ipoint in xrange(BRGraphPlot[mode].GetN() ):
		x, y = Double(0),Double(0)
		BRGraphPlot[mode].GetPoint(ipoint,x,y)
		BRGraphPlot[mode].SetPoint(ipoint,x,y/centralValue)


	BRGraphPlot[mode].SetLineColor(icolor)
	BRGraphPlot[mode].SetMarkerColor(icolor)
	BRGraphPlot[mode].SetMarkerStyle(imarker)
	BRGraphPlot[mode].SetMarkerSize(markersize)
	icolor = icolor+1
	markersize -= .1

	BRGraphPlot[mode].SetTitle(mode)
	BRGraphPlot[mode].Draw(drawopts)
	drawopts = "LP"

canvas.SetLogy()
# canvas.SetLogx()
canvas.SetGrid()

# BRGraphPlot["BB"].SetMinimum(0.000000000000001)
# BRGraphPlot["BB"].SetMaximum(1)


BRGraphPlot["TAU_TAU"].SetMinimum(0.03)
BRGraphPlot["TAU_TAU"].SetMaximum(10)


# BRGraphPlot["TAU_TAU"].GetXaxis().SetRangeUser(0.0000000,.5)
# BRGraphPlot["TAU_TAU"].GetXaxis().SetTitle("#alpha_{s}")
BRGraphPlot["TAU_TAU"].GetYaxis().SetTitle("BR/#hat{BR}")


canvas.BuildLegend(0.55,0.2,0.68,0.45)

# label = TLatex(0.06,10e-14,"M_{H} = 125 GeV")
label = TLatex(0.06,0.05,"M_{H} = #hat{M}_{H}(V)")
label.Draw()

canvas.SaveAs("BRvsScanParam.eps")

########################################################
