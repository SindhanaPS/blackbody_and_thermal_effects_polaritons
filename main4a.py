import numpy as np
import math as m
import matplotlib.pyplot as plt
from func import *
from scipy.optimize import fsolve,root_scalar,brentq

R = 8.3145/4184     # in kcal/(K.mol)

######### Parameters #########
k0 = 0.54                        # in M^(-1)s^(-1)
C = 3                            # concentration in M=mol/L
DHrxn = -20.5                    # enthalpy change during reaction in kcal/mol
Tout = 295                       # in K
RTout = 8.3145*Tout/4184         # in kcal/mol
Ea = 6.7                         # in kcal/mol

kconv = 2.5                      # thermal conductivity due to convection in W/(m^2.K)
fact = 1
##############################

kout = k0*np.exp(-Ea/RTout)
print(kout)

const1 = -2.7*10**(-8)*Tout**4/(k0*C**2*DHrxn)
const2 = -0.48*fact*kconv*Tout/(k0*C**2*DHrxn)
c1 = Ea/RTout

data1 = np.loadtxt("LzEtot_1.txt")
data2 = np.loadtxt("LzEtot_4.txt")

Lz = data1[:,0]          # in mu m
wcm = data1[:,1]         # in cm^-1
Ecav = data1[:,2]
EmolCaF2 = data2[:,2]

EcavbyLz = np.divide(Ecav,Lz)
EmolbyLz = np.divide(EmolCaF2,Lz)

RTmmol = np.zeros_like(Lz)
RTmcav = np.zeros_like(Lz)

for i in range(Lz.size):
   c2mol = const1*EmolbyLz[i]
   c2cav = const1*EcavbyLz[i]
   c3mc = const2*np.divide(1,Lz[i])
   RTmmol[i] = brentq(f, 1.00000000001, 1000, args=(c1, c2mol, c3mc))
   RTmcav[i] = brentq(f, 1.00000000001, 1000, args=(c1, c2cav, c3mc))

RTmmol = RTmmol*RTout
RTmcav = RTmcav*RTout

Tmmol = RTmmol/R     # in K
Tmcav = RTmcav/R     # in K

kmol = k0*np.exp(-np.divide(Ea,RTmmol))
kcav = k0*np.exp(-np.divide(Ea,RTmcav))

np.savetxt("kmol.txt",np.transpose([Lz,wcm,EmolCaF2,Tmmol,kmol]))
np.savetxt("kcav.txt",np.transpose([Lz,wcm,Ecav,Tmcav,kcav]))
