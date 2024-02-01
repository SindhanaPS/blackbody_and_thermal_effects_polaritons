import numpy as np
import math as m
import matplotlib.pyplot as plt
from func import *

flag1 = 1
flagT = 0

if flagT == 0:
   kT = 25.6                # in meV   -   T=298K
elif flagT == 1:
   kT = 30.1                # in meV   -   T=350K

############### constants ###################

hbarc = 197.3*10**(-3)              # in meV.mm
c = 299                             # in mm.GHz or mu m.THz
ccm = c*500/14.98962                # in mu m.cm^-1

#############################################

n0 = 1
nL = 1
dCaF2 = 2000           # in mu m

### Molecule
wp = 100               # in cm^-1
Gamma = 200            # in cm^-1

############################### Parameters ####################################
#### flag1 = 0    -    absorbing CaF2 and Au mirror + molecules + spacer
#### flag1 = 1    -    absorbing CaF2 and Au mirror + molecules 
#### flag1 = 2    -    absorbing CaF2 and Au mirror
#### flag1 = 3    -    absorbing CaF2 + molecules + spacer
#### flag1 = 4    -    absorbing CaF2 + molecules

#### flag1 = 5    -    non-absorbing CaF2 and mirror + molecules + spacer
#### flag1 = 6    -    non-absorbing CaF2 and mirror + molecules 
#### flag1 = 7    -    non-absorbing CaF2 and mirror
#### flag1 = 8    -    non-absorbing CaF2 + molecules + spacer
#### flag1 = 9    -    non-absorbing CaF2 + molecules

#### flag1 = 10   -    non-absorbing CaF2 and Au mirror + molecules + spacer

###############################################################################

if flag1 < 5:
   ### CaF2
   flagAbs = 1

   ### Molecule
   w0 = 2000              # in cm^-1
   einf = 1.4**2          # dimensionless

   ### Mirror
   wpM = 68557            # in cm^-1
   GammaM = 2382          # in cm^-1

   if flag1 == 0:            ####   absorbing CaF2 and Au mirror + molecules + spacer
      print("flag1 = 0    -    absorbing CaF2 and Au mirror + molecules + spacer")

      ### Molecule
      f = 20                 # dimensionless
      L0 = 3.0               # in mu m      # Molecular layer

      ### Mirror
      dM = 18*10**(-3)       # in mu m

   elif flag1 == 1:          ####   absorbing CaF2 and Au mirror + molecules
      print("flag1 = 1    -    absorbing CaF2 and Au mirror + molecules")

      ### Molecule
      f = 5                  # dimensionless

      ### Mirror
      dM = 18*10**(-3)       # in mu m

   elif flag1 == 2:          ####   absorbing CaF2 and Au mirror
      print("flag1 = 2    -    absorbing CaF2 and Au mirror")

      ### Molecule
      f = 0                  # dimensionless

      ### Mirror
      dM = 18*10**(-3)       # in mu m

   elif flag1 == 3:          ####   absorbing CaF2 + molecules + spacer
      print("flag1 = 3    -    absorbing CaF2 + molecules + spacer")

      ### Molecule
      f = 20                 # dimensionless
      L0 = 3.0               # in mu m      # Molecular layer

      ### Mirror
      dM = 0                 # in mu m

   elif flag1 == 4:          ####   absorbing CaF2 + molecules
      print("flag1 = 4    -    absorbing CaF2 + molecules")

      ### Molecule
      f = 5                  # dimensionless

      ### Mirror
      dM = 0                 # in mu m

elif flag1 >= 5 and flag1<10: 
   ### CaF2
   flagAbs = 0

   ### Mirror
   wpM = 20000            # in cm^-1
   GammaM = 0             # in cm^-1

   ### Molecule
   w0 = 1000              # in cm^-1
   einf = 1.4**2          # dimensionless

   if flag1 == 5:          ####   non-absorbing CaF2 and mirror + molecules + spacer
      print("flag1 = 5    -    non-absorbing CaF2 and mirror + molecules + spacer")

      ### Molecule
      f = 20                 # dimensionless
      L0 = 3.0               # in mu m      # Molecular layer

      ### Mirror
      dM = 18*10**(-3)       # in mu m

   elif flag1 == 6:          ####   non-absorbing CaF2 and Au mirror + molecules
      print("flag1 = 6    -    non-absorbing CaF2 and mirror + molecules")

      ### Molecule
      f = 5                  # dimensionless

      ### Mirror
      dM = 18*10**(-3)       # in mu m

   elif flag1 == 7:          ####   non-absorbing CaF2 and mirror
      print("flag1 = 7    -    non-absorbing CaF2 and mirror")

      ### Molecule
      f = 0                  # dimensionless

      ### Mirror
      dM = 18*10**(-3)       # in mu m

   elif flag1 == 8:          ####   non-absorbing CaF2 + molecules + spacer
      print("flag1 = 8    -    non-absorbing CaF2 + molecules + spacer")

      ### Molecule
      f = 20                 # dimensionless
      L0 = 3.0               # in mu m      # Molecular layer

      ### Mirror
      dM = 0                 # in mu m

   elif flag1 == 9:          ####   non-absorbing CaF2 + molecules
      print("flag1 = 9    -    non-absorbing CaF2 + molecules")

      ### Molecule
      f = 5                  # dimensionless

      ### Mirror
      dM = 0                 # in mu m

elif flag1 == 10:
   print("flag1 = 5    -    non-absorbing CaF2 and Au mirror + molecules + spacer")

   ### CaF2
   flagAbs = 0

   ### Mirror
   wpM = 68557            # in cm^-1
   GammaM = 2382          # in cm^-1
   dM = 18*10**(-3)       # in mu m

   ### Molecule
   w0 = 1000              # in cm^-1
   einf = 1.4**2          # dimensionless
   f = 20                 # dimensionless
   L0 = 3.0               # in mu m      # Molecular layer

#############################################

xmax = 17.5
x = np.linspace(0.01,xmax, num=400)

thetamax = m.pi*0.99/2
theta = np.linspace(0,thetamax, num =400)

Wclist = np.linspace(25,1150,num=50)
Llist = np.divide(ccm, 2*Wclist*1.43)
Etot = np.zeros(Llist.size, dtype=float)

k = x*kT*10**(-3)/hbarc       # in um m^-1
X,Y = np.meshgrid(theta,ccm*k/(2*m.pi))

np.savetxt('fig3_X.txt', X)
np.savetxt('fig3_Y.txt', Y)

i = 0
for Lz in Llist:

   if flag1 == 1 or flag1 == 2 or flag1 == 4 or flag1 == 6 or flag1 == 7 or flag1 == 9:
      L0 = Lz

   Etot[i], E, Es, Ep, Rs, Rp, Ts, Tp = Emissivitytot(Lz, L0, wp, f, w0, Gamma, einf, kT, wpM, GammaM, dM, n0, nL, dCaF2, flagAbs, x, theta)

   np.savetxt('fig3_E_'+str(i)+'_'+str(flag1)+'.txt', E)
   np.savetxt('fig3_Es_'+str(i)+'_'+str(flag1)+'.txt', Es)
   np.savetxt('fig3_Ep_'+str(i)+'_'+str(flag1)+'.txt', Ep)
   np.savetxt('fig3_Rs_'+str(i)+'_'+str(flag1)+'.txt', Rs)
   np.savetxt('fig3_Rp_'+str(i)+'_'+str(flag1)+'.txt', Rp)
   np.savetxt('fig3_Ts_'+str(i)+'_'+str(flag1)+'.txt', Ts)
   np.savetxt('fig3_Tp_'+str(i)+'_'+str(flag1)+'.txt', Tp)
   print(Lz,L0,Wclist[i],Etot[i])
   i = i + 1

np.savetxt('LzEtot_'+str(flag1)+'.txt',np.transpose([Llist,Wclist,Etot]))
