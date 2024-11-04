#! /Library/Frameworks/Python.framework/Versions/3.12/bin/python
import numpy as np
from numpy.linalg import pinv
import matplotlib.pyplot as plt

#  Variables to set in "main":
#	rtDir = directory to find mh_xxx.dat files, TotArr and Gasum files
#	Tarr = temperature range over which MH code should be valid ... 
#		... this generally matches the simulation temp range + a bit on each side
#	binTmp = one of the "xxx"'s from the mh_xxx.dat files ... doesn't matter which
#	nat = number of atoms
#       sjInd = the indices of the configurational entropy array (Sj) -- I set this range to exclude the "noisy" edge bins 
#		... if Vj is 0-999 ... then sjInd might be [50:950] ... just something that excludes the edges


#########################################
#########################################
def canon_calc_c23(nij, Nt, Nv, Ti, Vj, sjInd, kb):
    # Implementation as before
    fc = np.zeros(Nt)
    for i in range(Nt):
        beta_i = 1.0 / (kb * Ti[i])
        for j in range(sjInd[0], sjInd[1] + 1):
            if nij[i, j] > 0.0:
                c2 = nij[i, j] * (np.log(nij[i, j]) + beta_i * Vj[j])
                denom = np.sum(nij[:, j])
                c3 = 0.0
                for k in range(Nt):
                    if nij[k, j] > 0.0:
                        beta_k = 1.0 / (kb * Ti[k])
                        num = nij[i, j] * nij[k, j] * (np.log(nij[k, j]) + beta_k * Vj[j])
                        c3 += num / denom
                fc[i] += (c2 - c3)
    return fc


#########################################
#########################################
def canon_calc_d(nij, Nt, Nv, Ti, Vj, sjInd, kb):
    # Implementation as before
    D = np.zeros((Nt, Nt))
    for k in range(Nt):
        for i in range(Nt):
            for j in range(sjInd[0], sjInd[1] + 1):
                denom = np.sum(nij[:, j]) 
                if denom > 0.0:
                    if (k==i):
                        cint = 1
                    else: 
                        cint = 0
                    D[k, i] += nij[i, j] * (cint - nij[k, j] / denom)
    return D


#########################################
#########################################
def canon_calc_sj(Ak, nij, Nt, Nv, Ti, Vj, sjInd, kb):
    # Implementation as before
    fc = np.zeros(Nv)
    for j in range(sjInd[0], sjInd[1] + 1):
        num1 = 0.0
        num2 = 0.0
        denom = 0.0
        for i in range(Nt):
            if nij[i, j] > 0.0:
                beta_i = 1.0 / (Ti[i] * kb)
                num1 += nij[i, j] * (np.log(nij[i, j]) + beta_i * Vj[j])
                num2 += nij[i, j] * Ak[i]
                denom += nij[i, j]
        if denom > 0.0:
            fc[j] = (num1 - num2) / denom
    return fc


#########################################
#########################################
def canon_calc_ut(Tarr, Zt, Nv, Sj, Vj, sjInd, kb):
    # Implementation as before
    szT = Tarr.size
    fc = np.zeros(szT)
    fc2 = np.zeros(szT)
    for i in range(szT):
        num = 0.0
        num2 = 0.0
        for j in range(sjInd[0], sjInd[1] + 1):
            num += np.exp(Sj[j] - Vj[j] / (kb * Tarr[i])) * Vj[j]
            num2 += np.exp(Sj[j] - Vj[j] / (kb * Tarr[i])) * Vj[j]**2
        fc[i] = num / Zt[i]
        fc2[i] = num2 / Zt[i]
    return fc, fc2


#########################################
#########################################
def canon_calc_zt(Tarr, Nv, Sj, Vj, sjInd, kb):
    # Implementation as before
    szT = Tarr.size
    fc = np.zeros(szT)
    for i in range(szT):
        for j in range(sjInd[0], sjInd[1] + 1):
            fc[i] += np.exp(Sj[j] - Vj[j] / (kb * Tarr[i]))
    return fc



##################################################################################
##################################################################################
# Main script 
if __name__ == "__main__":

    rtDir = '/Users/gracie/Desktop/mhCvT/pyScr'				# set this to a root directory where you will keep the files
    tail = 'triL'
    if tail == 'biL':
        Tarr = np.arange(250, 401)
        binTmp = 280
        nat = 50
        sjInd = (50, 950)
    elif tail == 'triL':
        Tarr = np.arange(370, 441)
        binTmp = 370
        nat = 96
        sjInd = (100, 900)

    kb = 8.6173303e-5
    binFil = f"{rtDir}/mh_{binTmp}.dat"					# this is just a dummy file that holds the Vj's (first column)
									# ... of this file, output from mh.py code
    tavgFil = f"{rtDir}/Gasum_{tail}.dat"				# The average energy and temperature of each run
    totNijFil = f"{rtDir}/TotArr_{tail}.dat"				# the total array file that will have the "nij" data...
									# ... also output from the mh.py code

    nij = np.loadtxt(totNijFil)  # Import nij data
    teArr = np.loadtxt(tavgFil)  # Import T_i data
    meh = np.loadtxt(binFil)     # Import Vj data

    Ei = teArr[:, 1]
    Ti = teArr[:, 0]
    Vj = meh[:, 0]
    Vj -= np.min(Vj)

    Nt, Nv = nij.shape

    # Calculate c23 and D
    c23 = canon_calc_c23(nij, Nt, Nv, Ti, Vj, sjInd, kb)
    D = canon_calc_d(nij, Nt, Nv, Ti, Vj, sjInd, kb)

    # Solve for Ak
    print(D)
    Ak = pinv(D) @ c23

    # Calculate Sj
    Sj = canon_calc_sj(Ak, nij, Nt, Nv, Ti, Vj, sjInd, kb)

    # Plotting (if needed)
    plt.plot(range(sjInd[0], sjInd[1] + 1), Sj[sjInd[0]:sjInd[1] + 1])
    plt.xlabel("Index")
    plt.ylabel("Sj")
    plt.title("Sj Plot")
    plt.show()

    # Calculate Zt, Ut, and Ut2
    Zt = canon_calc_zt(Tarr, Nv, Sj, Vj, sjInd, kb)
    Ut, Ut2 = canon_calc_ut(Tarr, Zt, Nv, Sj, Vj, sjInd, kb)

    # Calculate Cv and normalized Cv
    Cv = np.zeros(Tarr.size - 1)
    c0 = 3 * nat * kb
    for i in range(Tarr.size - 1):
        Cv[i] = 3.0 * nat * kb / 2.0 + 1.0 / (kb * Tarr[i]**2) * (Ut2[i] - Ut[i]**2)
        Cv[i] = Cv[i] / c0

    nCv = np.zeros(Tarr.size - 1)
    for i in range(Tarr.size - 1):
        nCv[i] = 3.0 * nat * kb / 2.0
        nCv[i] = nCv[i] / c0 + 0.5

    # Plot Cv
    plt.plot(Tarr[:-1], Cv)
    plt.xlabel("Temperature")
    plt.ylabel("Cv")
    plt.title("Cv vs Temperature")
    plt.show()

    # Write Cv to a file
    fwstr = f"{rtDir}/newCvT_{tail}.dat"
    with open(fwstr, 'w') as fw:
        for j in range(Tarr.size - 1):
            fw.write(f"{Tarr[j]:4.4f}   {Cv[j]:10.6f}\n")

    print("Cv calculations and output complete.")

    # Write Sj to a file
    fwstr = f"{rtDir}/Ut_{tail}.dat"
    with open(fwstr, 'w') as fw:
        for j in range(Tarr.size - 1):
            fw.write(f"{j}   {Ut[j]:10.6f}\n")

    print("Cv calculations and output complete.")







