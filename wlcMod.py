# Simulation of filament as a WLC chain with heterogeneous bending stiffnessess along its length.
# Author: Ilya M. Beskin
# Date: Jun 2, 2024

from matplotlib import pyplot as plt
import numpy as np
import random as rand
from scipy.optimize import curve_fit

#This method will be further developed in the future.
def genRfromLinSens(strt_nm, decPow, numSegs):
	rArr = np.empty(numSegs)
	for i in range(numSegs):
		#if(i<=160):
		#	rArr[i] = 30e-9
		#else:
		#	rArr[i] = 30e-9-0.05*(i-160)
		rArr[i] = 1e-9*strt_nm*(((numSegs-i)/numSegs)**decPow)
		#rArr[i] = 1e-9*strt_nm*(numSegs-i)/numSegs
	return rArr

# genRfromLinDensRate generates a 1D array containing the radius [nm] at each point along the fibril assuming
#		that we are studying an alpha tip.
#	linDensRate 	->
#	numSegs 		-> Number of segments in the chain.
#	segLen			-> Length of each chain segment [m]
def genRfromLinDensRate(linDensRate, numSegs, segLen):
	
	rho = 0.3333	#Density of Collagen [kDa/nm^3]
					#linDensRate [kDa/nm^2]

	rArr = np.empty(numSegs)
	for i in np.arange(numSegs-1, -1, -1):
		rArr[i] = 1e-9*np.sqrt(linDensRate/rho*1000*(numSegs-i)*segLen*1e6/np.pi)	#[m]

	return rArr

# K_from_R generates a spring stiffness based on the Euler-Bernoulli Beam Theory prediction for a cylinder
#		 given a specific radius r.
#	gamma 	-> The Young's Modulus [MPa]
#	r 		-> Radius of cylinder [m]
#	segLen	-> Length of each chain segment [m]
def K_from_R(gamma,r, segLen):
	gamma *= 1e6	#https://doi.org/10.1089/ten.teb.2019.0243
	return gamma*np.pi*(r**4)/(4*segLen)

# genKArr generates a 1D array of spring coefficients.
#	gamma_MPa	-> Young's Modulus [MPa]
#	rArr 		-> 1D array of radii [m] along chain.
#	segLen		-> Length of each chain segment [m]
def genKArr(gamma_MPa, rArr, segLen):
	kArr = np.empty(len(rArr))
	for i in range(len(rArr)):
		kArr[i] = K_from_R(gamma_MPa, rArr[i], segLen)

	return kArr

# genDeflAng generates a deflection angle given a spring coefficient by picking
#		random value from a Gaussian distribution.
#	spK	-> spring coefficient given a specific chain length.
def genDeflAng(spK):
	kB = 1.38e-23	#Boltzmann Constant [J/K]
	T = 300			#Temperature [K]
	
	sig = np.sqrt(kB*T/spK)	#Sets width of Gaussian distribution.

	return rand.gauss(0, sig)

# genChain generates a 1D array of lateral deflections.
#	numSeg 	-> Number of segments in the chain.
#	segLen	-> Length of each chain segment [m]
#	spKArr	-> Array of spring coefficients
def genChain(numSeg, segLen, spKArr):
	chainPos = np.empty(numSeg)
	maxAng = 0.15 				#[Rad]
	endPos = numSeg*segLen 		#Stores pos at which angle of chain exceeds approx exceeds "maxAng" rad.
	endPosReached = 0

	chainPos[0] = 0
	chainPos[1] = 0
	
	for i in range(2, numSeg):
		chainPos[i] = 2*chainPos[i-1]-chainPos[i-2]+genDeflAng(spKArr[i])*segLen
		if( np.absolute((chainPos[i]-chainPos[i-1])/segLen) > maxAng and endPosReached == 0 ):
			endPos = i*segLen
			endPosReached = 1
			#print(endPos)

	return chainPos, endPos

# For fit function.
def pow(x, A, P):
	return A*(x**P) 

# genManyChains generates a 1D array containing the standard deviation of all the latteral
#		delfections at each position along the chain.
#	gamma_MPa	-> Young's Modulus [MPa]
#	numSegs 	-> Number of segments in each chain.
#	segLen		-> Length of each chain segment [m]
#	graphIt		-> Boolean variable. If true, this function will graph all the chains and the std dev array.
def genManyChains(gamma_MPa, numSegs, segLen, graphIt):
	#numSegs = 500
	#segLen = 1e-8	#10nm
	#rArr = genRfromLinSens(30, 1/2, numSegs)	#Radius at each joint.
	linDensRate = 0.2 #Rate at which linear mass density accumulates at an alpha end [kDa/nm^2]
	rArr = genRfromLinDensRate(linDensRate, numSegs, segLen)
	if(graphIt==1):
		plt.plot(rArr)
		plt.xlabel("Distance From Attachment Point [m]")
		plt.ylabel("Radius [m]")
		plt.show()
	kArr = genKArr(gamma_MPa, rArr, segLen)

	if(graphIt==1):
		plt.plot(kArr)
		plt.xlabel("Distance From Attachment Point [m]")
		plt.ylabel("Torsional Spring Coef [N*m]")
		plt.show()


	numChains = 500

	xLocs = np.arange(numSegs)*segLen

	allTheChains = np.empty([numChains, numSegs])
	endPosArr = np.empty(numChains)
	
	if(graphIt==1):
		fig, (ax1, ax2) = plt.subplots(2)
	
	for i in range(numChains):
		allTheChains[i], endPosArr[i] = genChain(numSegs, segLen, kArr)
		
		if(graphIt==1):
			ax1.plot(xLocs, allTheChains[i])

	#print(min(endPosArr))
	stdDevDisp = np.std(allTheChains, axis=0)

	#Power law fit
	fitCutoff = 4e-6 #min(endPosArr)
	popt, pcov = curve_fit(pow, xLocs[xLocs<fitCutoff], stdDevDisp[xLocs<fitCutoff])
	#print(popt)

	if(graphIt==1):
		ax2.plot(xLocs, stdDevDisp, label="S.D. from "+str(numChains)+" WLC", linewidth = 5)
		ax2.plot(xLocs[xLocs<fitCutoff], pow(xLocs[xLocs<fitCutoff], popt[0], popt[1],), 'r-', label="Power law fit: Power = "+str("{:.2f}".format(popt[1])))
		ax2.legend(loc="upper left")

		ax1.set_ylabel("Deflection [m]")
		ax2.set_xlabel("Distance From Attachment Point [m]")
		ax2.set_ylabel("Std. Dev. of Deflection [m]")

		plt.show()
		fig.savefig('STD_Dev.svg', format='svg', dpi=1200)

	return popt

# Running main finds exponential fits to the fluctuation amplitudes as a function of position
#		for a variety of elastic moduli.
def main():
	gammaArr = np.array([10, 20, 50, 100, 200, 500, 1000, 2000, 5000])
	#segNumArr = np.array([10, 20, 30, 40, 50, 60, 75, 100, 200, 300, 400, 500, 750, 1000])
	AconstArr = np.empty(gammaArr.size)
	AuncertArr = np.empty(gammaArr.size)
	BconstArr = np.empty(gammaArr.size)
	BuncertArr = np.empty(gammaArr.size)
	numLoops = 49
	fitArr = np.empty((numLoops,2))
	for j in range(gammaArr.size):
		print(gammaArr[j])
		for i in range(numLoops):
			fitArr[i] = genManyChains(gammaArr[j], 500, 1e-8, graphIt=0)
			#fitArr[i] = genManyChains(100, segNumArr[j], 5e-6/segNumArr[j], graphIt=0)
			print(fitArr[i])

		AconstArr[j] = np.mean(fitArr[:,0])
		AuncertArr[j] = np.std(fitArr[:,0])/np.sqrt(numLoops)
		BconstArr[j] = np.mean(fitArr[:,1])
		BuncertArr[j] = np.std(fitArr[:,1])/np.sqrt(numLoops)

	fig, (ax1, ax2) = plt.subplots(2)
	ax1.errorbar(gammaArr, AconstArr, AuncertArr, marker='o', color='black', capsize=5, capthick=1)
	ax2.errorbar(gammaArr, BconstArr, BuncertArr, marker='o', color='black', capsize=5, capthick=1)
	plt.show()

	data = np.transpose(np.array([gammaArr, AconstArr, AuncertArr, BconstArr, BuncertArr]))

	np.savetxt('gaussChain.csv',data,delimiter=',', header="El. Mod (MPa), A, delA, B, delB")


#main()
print(genManyChains(100, 500, 1e-8, graphIt=1))