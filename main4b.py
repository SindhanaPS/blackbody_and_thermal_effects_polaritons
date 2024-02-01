import numpy as np
import math as m
import matplotlib.pyplot as plt
from func import *
from scipy.optimize import fsolve,root_scalar,brentq

######### Parameters #########
Tout = 298            # in K
Tmol = 350            # in K

R = 10**(-2)          # in m 

rsm = 229             # thermal resistance in K.W^(-1)

kconv = 2.5           # thermal conductivity due to convection in W/(m^2.K)
fact = 1
##############################

data1 = np.loadtxt("LzEtot_5.txt")
data2 = np.loadtxt("LzEtot_8.txt")

Lz = data1[:,0]          # in mu m
wcm = data1[:,1]         # in cm^-1
Ecav = data1[:,2]
EmolCaF2 = data2[:,2]

Tsetmol = np.zeros_like(Lz)
Tsetcav = np.zeros_like(Lz)

for i in range(Lz.size):
   Tsetmol[i] = FindTset(EmolCaF2[i], kconv, Tout, Tmol, R, rsm, fact)
   Tsetcav[i] = FindTset(Ecav[i], kconv, Tout, Tmol, R, rsm, fact)

np.savetxt("Tset_mol_cav.txt",np.transpose([Lz,wcm,EmolCaF2,Ecav,Tsetmol,Tsetcav]))
