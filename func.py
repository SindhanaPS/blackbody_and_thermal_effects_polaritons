import numpy as np
import scipy.integrate as integrate
import scipy.special as special
import math as m
import matplotlib.pyplot as plt

hbarc = 197.3*10**(-3)        # in meV.mm
c = 300                       # in mm.GHz or mu m.THz
ccm = c*500/(14.98962)              # in mu m.cm^-1
sB = 5.67*10**(-8)            # in W.m^(-2)K^(-4)

datank = np.loadtxt('CaF2.txt')
wnk = datank[:,0]
nd = datank[:,1]
kd = datank[:,2]
datak = np.loadtxt('CaF2k.txt')
wk2 = datak[:,0]
kd2 = datak[:,1]
datan = np.loadtxt('CaF2n.txt')
wn2 = datan[:,0]
nd2 = datan[:,1]

def Theta(hw,kBT):

   y = np.divide(hw,np.exp(hw/kBT)-1)

   return(y)

def ucav(hw,kBT,n,hwcav):

   u = (n**3*hwcav/(c*m.pi**2*hbarc**2))*np.multiply(np.multiply(hw,Theta(hw,kBT)),0.5+np.floor(hw/hwcav))

   return(u)

def ufree(hw,kBT,n):

   u = (n**3/(c*m.pi**2*hbarc**2))*np.multiply(np.square(hw),Theta(hw,kBT))

   return(u)

def EAfromRT(R,T):
   A = 1 - R - T
   E = A

   return(E,A)

def Omega(k):
   # c in mm.GHz and k in mu m, omega in GHz 

   omega = ccm*m.sqrt(np.dot(k,k))/(2*m.pi)      # in cm^-1

   return(omega)

def Dielectric(w, wp, f, w0, Gamma, einf):
   # Returns a Lorentzian dielectric function

   eps = einf + np.divide(f*wp**2,w0**2-np.square(w)-1j*Gamma*w)

   return(eps)

def RefractiveIndex(w, wp, f, w0, Gamma, einf):
   # Returns the complex refractive index

   eps = Dielectric(w, wp, f, w0, Gamma, einf)
   n = np.sqrt(eps)

   return(n)

def RefractiveIndexDrude(w, wp, Gamma):
   # Returns the complex refractive index

   eps = 1 - np.divide(wp**2,np.square(w)+1j*Gamma*w)

   n = np.sqrt(eps)

   return(n)

def RefractiveIndexCaF2(w):
   # Returns complex refractive index of CaF2

   dw = 0.241057
   klow = m.floor(w/dw) 

   for k in range(klow, nd.size):
      if wnk[k]>=w:
         n = nd[k] + 1j*kd[k]
         break

   if w>897:
      n = 0

      for k in range(nd2.size):
         if wn2[k]>=w:
            n = nd2[k]
            break

      for k in range(kd2.size):
         if wk2[k]>=w:
            n = n + 1j*kd2[k]
            break
 
   return(n)

def Deltaij(x, n1, n2, ctheta1, ctheta2):
   # if x = 0 TE and x = 1 TM, interface matrix
   
   if x == 0:
      dp = 0.5*(1 + n2*ctheta2/(n1*ctheta1))
      dm = 0.5*(1 - n2*ctheta2/(n1*ctheta1))
   elif x == 1:
      dp = 0.5*(ctheta2/ctheta1 + n2/n1)
      dm = 0.5*(ctheta2/ctheta1 - n2/n1)

   M = np.array([[dp,dm],[dm,dp]])

   return(M)

def Piij(k ,L, ctheta, n):
   # Propagation matrix

   phi = 1j*n*m.sqrt(np.dot(k,k))*L*ctheta

   if phi.real < -150:
      a = np.exp(-1j*phi.imag)*1e+150
      b = np.exp(1j*phi.imag)*1e-150
   else:
      a = np.exp(-phi)
      b = np.exp(phi)

   M = np.array([[a,0],[0,b]])

   return(M)

def FindCtheta(n0, stheta0, n):

   n0sintheta2 = np.square(np.multiply(n0,stheta0))
   co = np.divide(np.sqrt(np.square(n)-n0sintheta2),n)

   nc = np.multiply(n,co)
   ctheta = np.zeros_like(nc)   

   if nc.imag>=0:
      ctheta = co
   else:
      ctheta = -co

   return(ctheta)


def RTfromM(M):

   T = np.square(1/np.abs(M[0,0]))
   R = np.square(np.abs(M[1,0]/M[0,0]))

   return(R,T)

def AuMirror(flag0, stheta0, kr, wpM, GammaM, dM, n0, n1, n2):
   # flag0 = 0 s-polarized
   # falg0 = 1 p-polarized

   # Refractive index calculation
   w = Omega(kr)
   nM = RefractiveIndexDrude(w, wpM, GammaM)

   # cos theta
   ctheta1 = FindCtheta(n0, stheta0, n1)
   cthetaM = FindCtheta(n0, stheta0, nM)
   ctheta2 = FindCtheta(n0, stheta0, n2)

   # Propagation matrix
   PM = Piij(kr, dM, cthetaM, nM)

   # 1M interface matrix
   D1M = Deltaij(flag0, n1, nM, ctheta1, cthetaM)
   # M2 interface matrix
   DM2 = Deltaij(flag0, nM, n2, cthetaM, ctheta2)

   # Total transfer matrix
   M = D1M @ PM @ DM2

   return(M)

def IncohConv(M):

   MI = np.zeros_like(M, dtype=float)

   MI[0,0] = np.abs(M[0,0])**2
   MI[0,1] = -np.abs(M[0,1])**2
   MI[1,0] = np.abs(M[1,0])**2
   MI[1,1] = (np.abs(np.linalg.det(M))**2-np.abs(M[0,1]*M[1,0])**2)/np.abs(M[0,0])**2

   return(MI)

def IncohPiij(k ,L, ctheta, n):
   # Propagation matrix

   phi = 1j*n*m.sqrt(np.dot(k,k))*L*ctheta

   if phi.real < -70:
      a = 1e+140
      b = 1e-140
   else:
      a = np.abs(np.exp(-phi))**2
      b = np.abs(np.exp(phi))**2

   M = np.array([[a,0],[0,b]])

   return(M)

def RTfromMIncoh(MI):

   T = 1/MI[0,0]
   R = MI[1,0]/MI[0,0]

   return(R,T)
 
def SpaceIncohEmissivityAumirrorCaF2Thetakr(theta, kr, Lz, L0, wp, f, w0, Gamma, einf, wpM, GammaM, dM, n0, nL, dCaF2, flagAbs):
   # flagAbs = 0    non-absorbing CaF2 layer
   # flagAbs = 1    absorbing CaF2 layer

   # layers
   # layer 0 air
   n0 = n0
   # layer 1 CaF2
   L1 = dCaF2
   # layer 2 non-absorbing
   n2 = 1.43
   L2 = (Lz - L0)/2
   # layer 3 molecule
   L3 = L0
   # layer 4 non-absorbing
   n4 = 1.43
   L4 = (Lz - L0)/2
   # layer 5 CaF2
   L5 = dCaF2
   # layer 6 air
   n6 = nL

   # 0 is for s and 1 is for p
   Ri = np.zeros(2)
   Ti = np.zeros(2)
   Ai = np.zeros(2)
   Ei = np.zeros(2)

   stheta0 = m.sin(theta)
   ctheta0 = m.cos(theta)

   # Refractive index calculation
   w = Omega(kr)
   n3 = RefractiveIndex(w, wp, f, w0, Gamma, einf)

   if flagAbs ==0:
      n1 = 1.43
   elif flagAbs == 1:
      n1 = RefractiveIndexCaF2(w)

   n5 = n1

   # cos theta
   ctheta1 = FindCtheta(n0, stheta0, n1)
   ctheta2 = FindCtheta(n0, stheta0, n2)
   ctheta3 = FindCtheta(n0, stheta0, n3)
   ctheta4 = FindCtheta(n0, stheta0, n4)
   ctheta5 = FindCtheta(n0, stheta0, n5)
   ctheta6 = FindCtheta(n0, stheta0, n6)


   for i in range(2):         # i=0 s-polarization and i=1 p-polarization

      #################### Incoherent interfaces #################
      # 01 interface matrix
      D01 = Deltaij(i, n0, n1, ctheta0, ctheta1)
      I01 = IncohConv(D01)

      # 56 interface matrix
      D56 = Deltaij(i, n5, n6, ctheta5, ctheta6)
      I56 = IncohConv(D56)

      ############ Propagaation in incoherent layers #############
      # Propagation matrix
      PI1 = IncohPiij(kr ,L1, ctheta1, n1)
      PI5 = IncohPiij(kr ,L5, ctheta5, n5)

      #################### Coherent layers #######################

      # mirror 12
      Phi12 = AuMirror(i, stheta0, kr, wpM, GammaM, dM, n0, n1, n2)
      # mirror 45
      Phi45 = AuMirror(i, stheta0, kr, wpM, GammaM, dM, n0, n4, n5)

      # 23 interface
      D23 = Deltaij(i, n2, n3, ctheta2, ctheta3)
      # 34 interface
      D34 = Deltaij(i, n3, n4, ctheta3, ctheta4)

      # Propagation matrix
      P2 = Piij(kr ,L2, ctheta2, n2)
      P3 = Piij(kr ,L3, ctheta3, n3)
      P4 = Piij(kr ,L4, ctheta4, n4)

      # Coherent transfer matrix
      Mc = Phi12 @ P2 @ D23 @ P3 @ D34 @ P4 @ Phi45

      I15 = IncohConv(Mc)

      ############################################################

      # Total energy transfer matrix
      M = I01 @ PI1 @ I15 @ PI5 @ I56

      Ri[i],Ti[i] = RTfromMIncoh(M)
      Ei[i],Ai[i] = EAfromRT(Ri[i],Ti[i])

   E = (Ei[0] + Ei[1])/2

   return(E, Ei[0], Ei[1], Ri[0], Ri[1], Ti[0], Ti[1])

def EmissivityAumirrorThetakr(theta, kr, Lz, wp, f, w0, Gamma, einf, wpM, GammaM, dM, n0, nL):

   # layer 1 molecule
   L1 = Lz

   # 0 is for s and 1 is for p
   Ri = np.zeros(2)
   Ti = np.zeros(2)
   Ai = np.zeros(2)
   Ei = np.zeros(2)

   stheta0 = m.sin(theta)
   ctheta0 = m.cos(theta)

   # Refractive index calculation
   w = Omega(kr)
   n1 = RefractiveIndex(w, wp, f, w0, Gamma, einf)
   n2 = nL

   # cos theta
   ctheta1 = FindCtheta(n0, stheta0, n1)
   ctheta2 = FindCtheta(n0, stheta0, n2)

   # Propagation matrix
   P1 = Piij(kr ,L1, ctheta1, n1)

   for i in range(2):         # i=0 s-polarization and i=1 p-polarization

      # mirror 01
      Phi01 = AuMirror(i, stheta0, kr, wpM, GammaM, dM, n0, n0, n1)
      # mirror 12
      Phi12 = AuMirror(i, stheta0, kr, wpM, GammaM, dM, n0, n1, n2)

      # Total transfer matrix
      M = Phi01 @ P1 @ Phi12

      Ri[i],Ti[i] = RTfromM(M)
      Ei[i],Ai[i] = EAfromRT(Ri[i],Ti[i])

   E = (Ei[0] + Ei[1])/2

   return(E, Ei[0], Ei[1])

def Emissivitytot(Lz, L0, wp, f, w0, Gamma, einf, kT, wpM, GammaM, dM, n0, nL, dCaF2, flagAbs, x, theta):
   # flagAbs = 0    non-absorbing CaF2 layer
   # flagAbs = 1    absorbing CaF2 layer

   # integrate over x=\hbar\omega/kT
   dx = x[1] - x[0]
   Etot = 0

   #Integrate over theta
   dtheta = theta[1]-theta[0]

   E = np.zeros((x.size,theta.size), dtype=float)
   Es = np.zeros((x.size,theta.size), dtype=float)
   Ep = np.zeros((x.size,theta.size), dtype=float)
   Rs = np.zeros((x.size,theta.size), dtype=float)
   Rp = np.zeros((x.size,theta.size), dtype=float)
   Ts = np.zeros((x.size,theta.size), dtype=float)
   Tp = np.zeros((x.size,theta.size), dtype=float)

   for i in range(x.size):
      xi = x[i]

      # kT in meV
      kr = xi*kT*10**(-3)/hbarc       # in um m^-1

      Ekr = 0

      for j in range(theta.size):
         thetaj = theta[j]
         E[i,j],Es[i,j],Ep[i,j],Rs[i,j],Rp[i,j],Ts[i,j],Tp[i,j] = SpaceIncohEmissivityAumirrorCaF2Thetakr(thetaj, kr, Lz, L0, wp, f, w0, Gamma, einf, wpM, GammaM, dM, n0, nL, dCaF2, flagAbs)

         integ1 = E[i,j]*m.sin(2*thetaj)
         Ekr = Ekr + integ1*dtheta

      integ2 = xi**3*Ekr/(m.exp(xi)-1)

      Etot = Etot + integ2*dx

   Etot = Etot*15/(m.pi)**4

   return(Etot, E, Es, Ep, Rs, Rp, Ts, Tp)

def FindTset(Etot, kconv, Tout, Tmol, R, rsm, fact):

   A = 2*m.pi*R**2
   Tset = Tmol + rsm*A*(Etot*sB*(Tmol**4-Tout**4) + kconv*fact*(Tmol-Tout))
 
   return(Tset)

def f(x, c1, c2, c3):
   y = np.exp(-c1/x) - c2*(x**4-1) - c3*(x-1)
   return(y)
