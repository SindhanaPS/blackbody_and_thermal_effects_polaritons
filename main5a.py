import numpy as np
from func import *
import math as m
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import rcParams

############### constants ###################

hbarc = 197.3*10**(-3)        # in meV.mm
R = 10                        # in mm
kBT = 25.7                    # in meV
n = 1

#############################################

hw = np.linspace(0.1,350, num=100)    # in meV
hwcIR = 61.99210                 # 500 cm^-1 or 14.98962 THz or 61.99210 meV
hwcTHz = 0.04136                 # 0.33356 cm^-1 or 10 GHz or 0.04136 meV

uf = ufree(hw,kBT,n)
uIR = ucav(hw,kBT,n,hwcIR)
uTHz = ucav(hw,kBT,n,hwcTHz)

hwcm = hw*500/61.99210
hwGHz = hw*10/0.04136

np.savetxt("fig5a.txt",np.transpose([hwcm,hwGHz,uf,uIR,uTHz]))
