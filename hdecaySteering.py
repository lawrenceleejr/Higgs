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
# [0.05,4.5577,4.0185,0.082,0.2435,8.855,0.0885,0.005282918,0.0046652,0.104428,0.1247905,0.2075025],
# [0.1,9.1154,8.037,0.164,0.487,17.71,0.177,0.01056584,0.0011663,0.208856,0.249581,0.1765228],
[0.15,13.6731,12.0555,0.246,0.7305,26.565,0.2655,0.01584876,0.0005183556,0.313284,0.3743715,0.1623447],
[0.2,18.2308,16.074,0.328,0.974,35.42,0.354,0.02113167,0.000291575,0.417712,0.499162,0.1535919],
[0.25,22.7885,20.0925,0.41,1.2175,44.275,0.4425,0.02641459,0.000186608,0.52214,0.6239525,0.1474266],
[0.3,27.3462,24.111,0.492,1.461,53.13,0.531,0.03169751,0.0001295889,0.626568,0.748743,0.1427449],
[0.35,31.9039,28.1295,0.574,1.7045,61.985,0.6195,0.03698043,9.520816e-05,0.730996,0.8735335,0.1390125],
[0.4,36.4616,32.148,0.656,1.948,70.84,0.708,0.04226335,7.289375e-05,0.835424,0.998324,0.1359337],
[0.45,41.0193,36.1665,0.738,2.1915,79.695,0.7965,0.04754627,5.759506e-05,0.939852,1.123115,0.1333289],


[0.4,36.4616,32.148,0.656,1.948,70.84,0.708,0.04226335,7.289375e-05,0.835424,0.998324,0.1359337],
[0.41,37.37314,32.9517,0.6724,1.9967,72.611,0.7257,0.04331993,6.938132e-05,0.8563096,1.023282,0.1353792],
[0.42,38.28468,33.7554,0.6888,2.0454,74.382,0.7434,0.04437651,6.611678e-05,0.8771952,1.04824,0.1348424],
[0.43,39.19622,34.5591,0.7052,2.0941,76.153,0.7611,0.0454331,6.307734e-05,0.8980808,1.073198,0.1343224],
[0.44,40.10776,35.3628,0.7216,2.1428,77.924,0.7788,0.04648968,6.024277e-05,0.9189664,1.098156,0.1338182],
[0.45,41.0193,36.1665,0.738,2.1915,79.695,0.7965,0.04754627,5.759506e-05,0.939852,1.123115,0.1333289],
[0.46,41.93084,36.9702,0.7544,2.2402,81.466,0.8142,0.04860285,5.511815e-05,0.9607376,1.148073,0.1328539],
[0.47,42.84238,37.7739,0.7708,2.2889,83.237,0.8319,0.04965943,5.279765e-05,0.9816232,1.173031,0.1323923],
[0.48,43.75392,38.5776,0.7872,2.3376,85.008,0.8496,0.05071602,5.062066e-05,1.002509,1.197989,0.1319436],
[0.49,44.66546,39.3813,0.8036,2.3863,86.779,0.8673,0.0517726,4.857559e-05,1.023394,1.222947,0.131507],
[0.5,45.577,40.185,0.82,2.435,88.55,0.885,0.05282918,4.6652e-05,1.04428,1.247905,0.1310821],
[0.51,46.48854,40.9887,0.8364,2.4837,90.321,0.9027,0.05388577,4.484045e-05,1.065166,1.272863,0.1306682],
[0.52,47.40008,41.7924,0.8528,2.5324,92.092,0.9204,0.05494235,4.31324e-05,1.086051,1.297821,0.1302649],
[0.53,48.31162,42.5961,0.8692,2.5811,93.863,0.9381,0.05599893,4.152011e-05,1.106937,1.322779,0.1298717],
[0.54,49.22316,43.3998,0.8856,2.6298,95.634,0.9558,0.05705552,3.999657e-05,1.127822,1.347737,0.1294882],
[0.55,50.1347,44.2035,0.902,2.6785,97.405,0.9735,0.0581121,3.855537e-05,1.148708,1.372696,0.1291138],
[0.56,51.04624,45.0072,0.9184,2.7272,99.176,0.9912,0.05916869,3.719069e-05,1.169594,1.397654,0.1287484],
[0.57,51.95778,45.8109,0.9348,2.7759,100.947,1.0089,0.06022527,3.58972e-05,1.190479,1.422612,0.1283914],
[0.58,52.86932,46.6146,0.9512,2.8246,102.718,1.0266,0.06128185,3.467004e-05,1.211365,1.44757,0.1280425],
[0.59,53.78086,47.4183,0.9676,2.8733,104.489,1.0443,0.06233844,3.350474e-05,1.23225,1.472528,0.1277014],


[0.6,54.6924,48.222,0.984,2.922,106.26,1.062,0.06339502,3.239722e-05,1.253136,1.497486,0.1273679],
[0.65,59.2501,52.2405,1.066,3.1655,115.115,1.1505,0.06867794,2.760473e-05,1.357564,1.622277,0.1258029],
[0.7,63.8078,56.259,1.148,3.409,123.97,1.239,0.07396086,2.380204e-05,1.461992,1.747067,0.1243879],
[0.75,68.3655,60.2775,1.23,3.6525,132.825,1.3275,0.07924378,2.073422e-05,1.56642,1.871857,0.1230989],
[0.8,72.9232,64.296,1.312,3.896,141.68,1.416,0.08452669,1.822344e-05,1.670848,1.996648,0.121917],
[0.85,77.4809,68.3145,1.394,4.1395,150.535,1.5045,0.08980961,1.614256e-05,1.775276,2.121439,0.1208273],
[0.9,82.0386,72.333,1.476,4.383,159.39,1.593,0.09509253,1.439877e-05,1.879704,2.246229,0.1198176],
[0.95,86.5963,76.3515,1.558,4.6265,168.245,1.6815,0.1003754,1.292299e-05,1.984132,2.37102,0.118878],
[1,91.154,80.37,1.64,4.87,177.1,1.77,0.1056584,1.1663e-05,2.08856,2.49581,0.118],

[2,182.308,160.74,3.28,9.74,354.2,3.54,0.2113167,2.91575e-06,4.17712,4.99162,0.1072922],
# [3,273.462,241.11,4.92,14.61,531.3,5.31,0.3169751,1.295889e-06,6.26568,7.48743,0.1018839],
# [4,364.616,321.48,6.56,19.48,708.4,7.08,0.4226335,7.289375e-07,8.35424,9.98324,0.09836598],
# [5,455.77,401.85,8.2,24.35,885.5,8.85,0.5282918,4.6652e-07,10.4428,12.47905,0.09580019],
# [9.95,906.9823,799.6815,16.318,48.4565,1762.145,17.6115,1.051301,1.178051e-07,20.78117,24.83331,0.08866784],


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



for SamplePoint in SamplePoints:


	ParameterDictionary = dict(zip(ParameterNameArray, SamplePoint))

	print ParameterDictionary

	endvalue = 2*125*ParameterDictionary['V']
	if endvalue<60:
		endvalue=60

	higgsMassStart = "  1.D0"
	higgsMassEnd = '{0:.10f}'.format(endvalue) + "D0" #"145.D0"
	nHiggsMassSteps = "200"
	

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

			BRLibrary[iHiggsMass]['UU'] = (3./100.)**2*table1[iHiggsMass]['SS']
			BRLibrary[iHiggsMass]['DD'] = (5./100.)**2*table1[iHiggsMass]['SS']
			BRLibrary[iHiggsMass]['EE'] = (0.511/105.7)**2*table1[iHiggsMass]['MU_MU']

			BRLibrary[iHiggsMass]['TT']      = table1[iHiggsMass]['TT'] 
			BRLibrary[iHiggsMass]['BB']      = table1[iHiggsMass]['BB'] 
			BRLibrary[iHiggsMass]['TAU_TAU'] = table1[iHiggsMass]['TAU_TAU']      
			BRLibrary[iHiggsMass]['MU_MU']   = table1[iHiggsMass]['MU_MU']   
			BRLibrary[iHiggsMass]['SS']      = table1[iHiggsMass]['SS']    
			BRLibrary[iHiggsMass]['CC']      = table1[iHiggsMass]['CC']    
			BRLibrary[iHiggsMass]['GG']      = table2[iHiggsMass]['GG']    
			BRLibrary[iHiggsMass]['GAM_GAM'] = table2[iHiggsMass]['GAM_GAM']         
			BRLibrary[iHiggsMass]['Z_GAM']   = table2[iHiggsMass]['Z_GAM']       
			BRLibrary[iHiggsMass]['WW']      = table2[iHiggsMass]['WW']    
			BRLibrary[iHiggsMass]['ZZ']      = table2[iHiggsMass]['ZZ']    


			productOfBRs = 1
			for imode in ['BB','CC','SS','TAU_TAU','MU_MU','GG','GAM_GAM','Z_GAM','WW','ZZ','UU','DD','EE']:
				# print '%s      %f'%(imode, BRLibrary[iHiggsMass][imode] )
				if BRLibrary[iHiggsMass][imode]<0:
					# print "SKIPPING---------"
					# continue
					BRLibrary[iHiggsMass][imode]=0
				# assert (0<=BRLibrary[iHiggsMass][imode]<=1)
				productOfBRs *= BRLibrary[iHiggsMass][imode]

			productOfBRsExclTop = productOfBRs
			productOfBRs = productOfBRs*table1[iHiggsMass]["TT"]

			if math.isnan(productOfBRsExclTop) or math.isnan(productOfBRs):
				productOfBRs = 0
				productOfBRsExclTop = 0


			print "PRODUCT IS............."
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

		print fit.GetMaximum()

		MHhatArray['MHhat'].append(fitMean)
		MHhatArray['P'].append(fit.GetMaximum() )
		MHhatArray['V'].append(ParameterDictionary[VARNAME])


		BRGraph['scanParam'].append(ParameterDictionary[VARNAME] )
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

