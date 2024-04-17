#!/usr/bin/env python3
import os
import subprocess
import shutil

def createInFile(iQ2, saveName):
	initF = open("int_tmp_valerii.dat", "w")
	headLines = ['10.604\n',
	str(iQ2)+'\n',
	'1.125 4.525\n',
	'100\n', saveName + '\n'
	]
	initF.writelines(headLines)
	
	initF.close()

def main():
	Q2bins = [2.55,2.99,3.49,4.08,4.78,5.59,6.53,7.64,8.94,10.4]
	
	iQ2steps = 40
	
	for iQ2 in range(9):
		delta_Q2 = Q2bins[iQ2+1] - Q2bins[iQ2]
		for iQ2_step in range(iQ2steps):
			tmp_Q2 = Q2bins[iQ2] + delta_Q2*iQ2_step/iQ2steps
			sName = 'out_valerii/q2_' + str(iQ2+1)+'/iq2step_' + str(iQ2_step)
			createInFile(tmp_Q2,sName)
			
			subprocess.call(["./distrib_x86_64 < int_tmp_valerii.dat"], shell=True)
		
main()
