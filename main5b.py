import numpy as np
from func import *
import math as m


hw = np.linspace(0.0012,350, num=100)    # in meV
hwcIR = 61.99210                 # 500 cm^-1 or 14.98962 THz or 61.99210 meV
hwcTHz = 0.04136                 # 0.33356 cm^-1 or 10 GHz or 0.04136 meV

ratiocav1 = np.divide(hwcIR,hw)*(np.floor(hw/hwcIR)+0.5)
ratiocav2 = np.divide(hwcTHz,hw)*(np.floor(hw/hwcTHz)+0.5)

hwcm = hw*500/61.99210
hwGHz = hw*10/0.04136

np.savetxt("fig5b.txt",np.transpose([hwcm,hwGHz,ratiocav1,ratiocav2]))
