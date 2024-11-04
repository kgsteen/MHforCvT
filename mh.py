import numpy as np
import sys
import os
import re

############## all variables for initialisation

stTmp = 280		# starting temperature
endTmp = 370		# ending temperature
dT = 30			# delta temperature
tempAr = np.arange(stTmp, endTmp + 1, dT)		# create the temperature array

binNum = 1000
tail = 'biL'

mhDir = '/path/to/OUTCAR/files'		# this is where you want to output the files - mh_ttt.dat, TotMh_biL.dat, etc.

# name the total-array file that will contain the sum of all the individual histograms, formatted for for the CvT analysis program
fwTotArr = '%s/TotArr_%s.dat' % (mhDir, tail)
# name the total-array file that will contain the sum of all the individual histograms, formatted for plotting
fwTotMh = '%s/TotMh_%s.dat' % (mhDir, tail)
# name the file that will have 2 columns:  avg total energy | avg temperature 
fwAvg = '%s/Gasum_%s.dat' % (mhDir, tail)


##################################
def mkEnArr(runTmp):                          # create the energy array for each temperature

	eArr = []		# energy-without entropy array
	etArr = []		# total energy (config + kin) array
	tArr = []		# temperature array

	# path to your OUTCAR file
	frstr = f"OUTCAR_{runTmp}"

	# Define regex patterns for each of the three matches
	energy_pattern = re.compile(r'without\s+entropy=\s+(-?\d+\.\d+)')
	etotal_pattern = re.compile(r'ETOTAL\s+=\s+(-?\d+\.\d+)\s+eV')
	ekin_lat_pattern = re.compile(r'temperature\s+(\d+\.\d+)\s+K')

	# Open and read the file line-by-line
	with open(frstr, 'r') as file:
		for line in file:
			# Match each line against the patterns and extract the relevant value
			energy_match = energy_pattern.search(line)
			etotal_match = etotal_pattern.search(line)
			ekin_lat_match = ekin_lat_pattern.search(line)
	
			if energy_match:	# and re.match(r'^-?\d+(\.\d+)?$', energy_match.group(1)):	
				eArr.append(float(energy_match.group(1)))
		
			if etotal_match:
				etArr.append(float(etotal_match.group(1)))
		
			if ekin_lat_match:
				tArr.append(float(ekin_lat_match.group(1)))

		# return 'energy  without entropy', total energy and temperature arrays 
		return eArr, etArr, tArr


################################################
################  Main  ########################
################################################

eTarr = []			# create a "holding" list to store the config. energies ...
				# ... list allows for different numbers of config. energies per run
				# ... in case not all runs have the same number of time steps 
with open(fwAvg,'w') as fw:
	for i, temp in enumerate(tempAr):

		eNarr, etotArr, tmpArr = mkEnArr(temp)  # function argument is the temperature of the run
			# the return arrays are (as np.arrays):
				# eNarr = configurational energy (without entropy)
				# eTotArr = total energy (config + kinetic)
				# tmpArr = kinetic energy at each time step

		eNarr.sort()		# sort the configurational energy array 
					# ... this makes setting histogram bin-range easier later
		eTarr.append(eNarr)	# append this temperature's config energies to eTarr[i] 
					#... final array will be: eTarr[i] as a list where i is indexed by number of temperatures

		# take the average total energy and temperature and write this to Gasum file
		fw.write('% 10.6f    % 10.6f\n' % (np.mean(tmpArr), np.mean(etotArr) ))

		print(temp)			# just keeping track of progress
		sys.stdout.flush()



#####################
#####################

nT = len(tempAr)				# number of temperatures

# histogram settings
# ... the energies will be binned in a range from binStart to binEnd
# ... binStart is determined by the min energy of the lowest temp configuration energies:  eTarr[0][0]
# ... binEnd is determined by the max energy of the highest temp configuration energies:  eTarr[nT-1][nPe-1]
# both are padded by 0.5 eV
binStart = eTarr[0][0] - 0.5 	
lenPE_highT = len(eTarr[nT-1])
binEnd = eTarr[nT-1][lenPE_highT-1] + 0.5
dBin = (binEnd - binStart)/((float)(binNum))	# the deltaE for each bin is set -- binNum is set at the beginning of the file
bins = np.arange(binStart, binStart + (binNum + 1) * dBin, dBin)	# set the bins for later histogramming

TCtbin = np.zeros(len(bins)-1)		# this will contain the total counts summed over ALL temperatures

with open(fwTotArr,'w') as fwTot:	# this will be the total array file that is formatted for CvT analysis

	for i, temp in enumerate(tempAr):

		Ctbin = np.zeros(len(bins)-1)
    
		hist, _ = np.histogram(eTarr[i], bins=bins)	# histogram for PE values at the current temperature
    
		Ctbin = hist			# assign those counts to Ctbin
		TCtbin += hist			# sum the counts to accumulate over all temperatures

		fwstr = '%s/mh_%d.dat' % (mhDir, temp)
		with open(fwstr,'w') as fw:
			for ct in range(binNum - 1, -1, -1):
				# Write to `fw` with bin values
				fw.write(f"{binStart + ct * dBin}\t{Ctbin[ct]}\n")

				# Append only the value to `fwTot`
				fwTot.write(f"{Ctbin[ct]}\t")

		fwTot.write("\n")

#########################
#########################
with open(fwTotMh,'w') as fwTmh:		# this is the total array file that is formatted for plotting
	for ct in range(binNum - 1, -1, -1):
		fwTmh.write(f"{binStart + ct * dBin}\t{TCtbin[ct]}\n")





