import numpy as np
import math as m
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import LogNorm
from scipy.signal import argrelextrema
from mpl_toolkits.axes_grid1 import make_axes_locatable

flagi = 7

#####################################################
#                Importing data                     #
#####################################################
if flagi == 1:
   n = 1.43
   Lz = 5                              # in mu m

elif flagi == 3:
    data = np.loadtxt('LzEtot_0.txt')
    Lz = np.zeros((10,data[:,0].size),dtype=float)
    wc = np.zeros((10,data[:,0].size),dtype=float)
    Etot = np.zeros((10,data[:,0].size),dtype=float)

    for flag1 in range(6):
       data = np.loadtxt('LzEtot_'+str(flag1)+'.txt')

       Lz[flag1] = data[:,0]
       wc[flag1] = data[:,1]
       Etot[flag1] = data[:,2]

    for flag1 in range(7,9):
       data = np.loadtxt('LzEtot_'+str(flag1)+'.txt')

       Lz[flag1] = data[:,0]
       wc[flag1] = data[:,1]
       Etot[flag1] = data[:,2]

    X = np.loadtxt("fig3_X.txt")
    Y = np.loadtxt("fig3_Y.txt")

    ilist5 = np.array([333])

    E5 = np.zeros((ilist5.size,np.size(X,0),np.size(X,1)))
    T5 = np.zeros((ilist5.size,np.size(X,0),np.size(X,1)))

    flag1 = 5
    j = 0
    for i in ilist5:
       E5[j] = np.loadtxt('fig3_E_'+str(i)+'_'+str(flag1)+'.txt')
       T5[j] = np.loadtxt('fig3_Ts_'+str(i)+'_'+str(flag1)+'.txt')
       j = j+1

    ilist0 = np.array([500])

    E0 = np.zeros((ilist0.size,np.size(X,0),np.size(X,1)))
    T0 = np.zeros((ilist0.size,np.size(X,0),np.size(X,1)))
   
    flag1 = 0
    j = 0
    for i in ilist0:
       E0[j] = np.loadtxt('fig3_E_'+str(i)+'_'+str(flag1)+'.txt')
       T0[j] = np.loadtxt('fig3_Ts_'+str(i)+'_'+str(flag1)+'.txt')
       j = j+1

elif flagi == 4:
   data1 = np.loadtxt("kcav.txt")
   data2 = np.loadtxt("kmol.txt")

   Lz = data1[:,0]
   wcm = data1[:,1]         # in cm^-1
   Ecav = data1[:,2]
   Tmcav = data1[:,3]
   kcav = data1[:,4]
   Emol = data2[:,2]
   Tmmol = data2[:,3]
   kmol = data2[:,4]

   Tout = 295            # in K

elif flagi == 5:
   data1 = np.loadtxt("Tset_mol_cav.txt")

   Lz = data1[:,0]
   wcm = data1[:,1]         # in cm^-1
   Emol = data1[:,2]
   Ecav = data1[:,3]
   Tsetmol = data1[:,4]
   Tsetcav = data1[:,5]

elif flagi == 6:
   data = np.loadtxt("fig5a.txt")

   wcm = data[:,0]         # in cm^-1
   wGHz = data[:,1]        # in GHz
   uf = data[:,2]          # free space
   uIR = data[:,3]         # IR cavity
   uTHz = data[:,4]        # THz cavity

   data = np.loadtxt("fig5b.txt")

   wcm = data[:,0]         # in cm^-1
   wGHz = data[:,1]        # in GHz
   ratiocav1 = data[:,2]    # ratio of photon density of states
   ratiocav2 = data[:,3]    # ratio of photon density of states
   hwcav1 = 500
   hwcav2 = 0.33356

elif flagi == 7:
    data = np.loadtxt('LzEtot_5.txt')
    Lz = np.zeros((10,data[:,0].size),dtype=float)
    wc = np.zeros((10,data[:,0].size),dtype=float)
    Etot = np.zeros((10,data[:,0].size),dtype=float)

    for flag1 in [5,8]:
       data = np.loadtxt('LzEtot_'+str(flag1)+'.txt')

       Lz[flag1] = data[:,0]
       wc[flag1] = data[:,1]
       Etot[flag1] = data[:,2]

#####################################################

######################################################
#                 Formatting                         #
######################################################

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 16}

rc('font', **font)

rcParams['text.latex.preamble'] = [
       r'\usepackage{physics}',
       r'\usepackage{amsmath}',
]

rcParams['axes.linewidth'] = 1

plt.rc('text.latex', preamble=r'\usepackage[cm]{sfmath}')
plt.rc('text', usetex=True)
######################################################

if flagi == 1:
   fig1, ax1 = plt.subplots(1)

   c = 300                       # in mm.GHz or mu m.THz
   ccm = c*500/(14.98962)              # in mu m.cm^-1 
   kpar = np.linspace(-1,1,num=200)    # in mu m^-1
   mmax = 2


   for i in range(mmax+1):
      w = (ccm/n)*np.sqrt(np.square(kpar)+(i*m.pi/Lz)**2)/(2*m.pi)
      ax1.plot(kpar,w, color = 'lightcoral', linewidth=3, zorder=-1)
      if i>0:
         ax1.text(x=1,y=w[199],s=r'TE$_%d$'%i+r'TM$_%d$'%i,fontsize=20)
      else:
         ax1.text(x=1,y=w[199],s=r'TM$_%d$'%i,fontsize=20)

   # Move left y-axis and bottom x-axis to centre, passing through (0,0)
   ax1.spines['left'].set_position('center')

   # Eliminate upper and right axes
   ax1.spines['right'].set_color('none')
   ax1.spines['top'].set_color('none')

   # Show ticks in the left and lower axes only
   ax1.xaxis.set_ticks_position('bottom')
   ax1.yaxis.set_ticks_position('left')

   ax1.set_ylim(bottom=0)
   ax1.set_xlabel(r'$k_{||}$  $(\mathrm{in}$ $\mu\mathrm{m}^{-1})$',zorder=1)
   ax1.set_ylabel(r'$\omega$  $(\mathrm{in}$ $\mathrm{cm}^{-1})$')

   ratio = 1.1
   x_left, x_right = ax1.get_xlim()
   y_low, y_high = ax1.get_ylim()
   ax1.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)

   plt.savefig('fig1b.pdf', bbox_inches='tight')
   plt.show()

elif flagi == 3:
   fig1, ax1 = plt.subplots(1)
   
   hw0 = 2000

   ax1.plot(wc[0], Etot[0], color = '#91bfdb',linewidth=3, label = 'cavity+molecules')
   ax1.plot(wc[0], Etot[2], color = '#fc8d59',linewidth=3, linestyle='dashed', label = 'cavity')

   ax2 = ax1.twiny()
   ax2.plot(wc[0]*0.02998, Etot[3], color = '#fc8d59',linewidth=3)
   ax1.plot(wc[0], Etot[3], color = '#fc8d59', linewidth=3, label = 'molecules')

   ax2.set_zorder(-1)

   ax2.set_xlabel(r'$\omega_{\mathrm{cav}}$ $(\mathrm{in}$ $\mathrm{THz})$',fontsize=20)
   ax1.set_xlabel(r'$\omega_{\mathrm{cav}}$ $(\mathrm{in}$ $\mathrm{cm}^{-1})$', fontsize=20)
   ax1.set_ylabel(r'$\mathrm{Emissivity}$  $\mathcal{E}_{\mathrm{tot}}$', fontsize=20)
   ax1.set_xlim(right=1200)

   ax1.axvline(x=1000,linestyle=(0, (5, 10)),color='black',linewidth=2)
   ax1.axvline(x=666.6,linestyle=(0, (5, 10)),color='black',linewidth=2)
   ax1.axvline(x=500,linestyle=(0, (5, 10)),color='black',linewidth=2)

   ax1.legend(fontsize=12, loc='lower right')
   plt.savefig('fig3b.pdf', bbox_inches='tight')
   plt.show()

   rc('xtick', labelsize=35) 
   rc('ytick', labelsize=35) 

   j = 0
   flag1 = 0
   for i in ilist0:
      fig1, ax1 = plt.subplots(1)
      ax1.set_ylabel(r'$\omega$ $(\mathrm{in}$ $\mathrm{cm}^{-1})$', fontsize=40)
      ax1.set_xlabel(r'$\theta$', fontsize=40)
      ax1.pcolor(X,Y,E0[j], rasterized=True)
      im = ax1.pcolor(X,Y,E0[j],rasterized=True)
      ax1.axhline(y=hw0,color='red',linewidth=2)
      ratio = 0.75
      x_leftb, x_rightb = ax1.get_xlim()
      y_lowb, y_highb = ax1.get_ylim()
      ax1.set_aspect(abs((x_rightb-x_leftb)/(y_lowb-y_highb))*ratio)
      cbar = plt.colorbar(im,shrink=0.8)
      plt.title(r'$\mathcal{E}(\mathbf{k})$',fontsize=40, pad=20)
      plt.savefig('fig3b1.png', bbox_inches='tight')
      plt.savefig('fig3b1.pdf', bbox_inches='tight')
      plt.show()

      plt.plot(Y[:,0],T0[j,:,0],color='black')
      ax1 = plt.gca()
      ax1.axvline(x=hw0,color='red',linewidth=2)
      ax1.set_ylabel(r'$\mathcal{T}(\theta=0)$', fontsize=40)
      ax1.set_xlabel(r'$\omega$ $(\mathrm{in}$ $\mathrm{cm}^{-1})$', fontsize=40)
      ratio = 0.35
      x_leftb, x_rightb = ax1.get_xlim()
      y_lowb, y_highb = ax1.get_ylim()
      ax1.set_aspect(abs((x_rightb-x_leftb)/(y_lowb-y_highb))*ratio)
      plt.savefig('fig3b2.pdf', bbox_inches='tight')
      plt.show()

      j = j+1

   rc('xtick', labelsize=16) 
   rc('ytick', labelsize=16) 

   fig1, ax1 = plt.subplots(1)

   ax1.plot(wc[0], Etot[1], color = '#91bfdb',linewidth=3, label = 'cavity+molecules')
   ax2 = ax1.twiny()
   ax2.plot(wc[0]*0.02998, Etot[4], color = '#fc8d59',linewidth=3)
   ax1.plot(wc[0], Etot[4], color = '#fc8d59', linewidth=3, label = 'molecules')

   ax2.set_zorder(-1)

   ax2.set_xlabel(r'$\omega_{\mathrm{cav}}$ $(\mathrm{in}$ $\mathrm{THz})$',fontsize=20)
   ax1.set_xlabel(r'$\omega_{\mathrm{cav}}$ $(\mathrm{in}$ $\mathrm{cm}^{-1})$', fontsize=20)
   ax1.set_ylabel(r'$\mathrm{Emissivity}$  $\mathcal{E}_{\mathrm{tot}}$', fontsize=20)
   ax1.legend(fontsize=17)
   plt.savefig('figS1a.pdf', bbox_inches='tight')
   plt.show()

   fig1, ax1 = plt.subplots(1)

   hw0 = 1000
   ax1.plot(wc[5], Etot[5], color = '#91bfdb',linewidth=3, label = 'cavity+molecules')
   ax1.plot(wc[5], Etot[7], color = '#fc8d59',linewidth=3, linestyle='dashed', label = 'cavity')
   ax2 = ax1.twiny()
   ax2.plot(wc[5]*0.02998, Etot[8], color = '#fc8d59',linewidth=3)
   ax1.plot(wc[5], Etot[8], color = '#fc8d59', linewidth=3, label = 'molecules')
   
   ax2.set_zorder(-1)

   ax1.axvline(x=1000,linestyle=(0, (5, 10)),color='black',linewidth=2)
   ax1.axvline(x=500,linestyle=(0, (5, 10)),color='black',linewidth=2)
   ax1.axvline(x=333.3,linestyle=(0, (5, 10)),color='black',linewidth=2)

   ax2.set_xlabel(r'$\omega_{\mathrm{cav}}$ $(\mathrm{in}$ $\mathrm{THz})$',fontsize=20)
   ax1.set_xlabel(r'$\omega_{\mathrm{cav}}$ $(\mathrm{in}$ $\mathrm{cm}^{-1})$', fontsize=20)
   ax1.set_ylabel(r'$\mathrm{Emissivity}$  $\mathcal{E}_{\mathrm{tot}}$', fontsize=20)
   ax1.set_xlim(right=1200)
   ax1.legend(fontsize=12,loc = 'upper right')

   plt.savefig('fig3a.pdf', bbox_inches='tight')
   plt.show()

   rc('xtick', labelsize=35) 
   rc('ytick', labelsize=35) 

   flag1 = 5
   j = 0
   for i in ilist5:
      fig1, ax1 = plt.subplots(1)
      ax1.set_ylabel(r'$\omega$ $(\mathrm{in}$ $\mathrm{cm}^{-1})$', fontsize=40)
      ax1.set_xlabel(r'$\theta$', fontsize=40)
      ax1.pcolor(X,Y,E5[j],rasterized=True)
      im = ax1.pcolor(X,Y,E5[j], rasterized=True)
      ax1.axhline(y=hw0,color='red',linewidth=2)
      ratio = 0.75
      x_leftb, x_rightb = ax1.get_xlim()
      y_lowb, y_highb = ax1.get_ylim()
      ax1.set_aspect(abs((x_rightb-x_leftb)/(y_lowb-y_highb))*ratio)
      cbar = plt.colorbar(im,shrink=0.8)
      plt.title(r'$\mathcal{E}(\mathbf{k})$',fontsize=40, pad=20)
      plt.savefig('fig3a1.png', bbox_inches='tight')
      plt.savefig('fig3a1.pdf', bbox_inches='tight')
      plt.show()

      plt.plot(Y[:,0],T5[j,:,0],color='black')
      ax1 = plt.gca()
      ax1.axvline(x=hw0,color='red',linewidth=2)
      ax1.set_ylabel(r'$\mathcal{T}(\theta=0)$', fontsize=40)
      ax1.set_xlabel(r'$\omega$ $(\mathrm{in}$ $\mathrm{cm}^{-1})$', fontsize=40)
      ratio = 0.35
      x_leftb, x_rightb = ax1.get_xlim()
      y_lowb, y_highb = ax1.get_ylim()
      ax1.set_aspect(abs((x_rightb-x_leftb)/(y_lowb-y_highb))*ratio)
      plt.savefig('fig3a2.pdf', bbox_inches='tight')
      plt.show()

      j = j+1

   rc('xtick', labelsize=16) 
   rc('ytick', labelsize=16) 

elif flagi == 4:
   fig,ax1 = plt.subplots()

   ax1.set_ylim([294.9,295.1])
   ax1.set_ylabel(r'$T_{\mathrm{mol}}$ (in K)')
   ax1.set_xlabel(r'$\omega_{\mathrm{cav}}$ $(\mathrm{in}$ $\mathrm{cm}^{-1})$')
   ax1.plot(wcm,Tmcav,color='limegreen',linewidth=2,label='cavity+molecules')
   ax1.plot(wcm,Tmmol,linestyle='dashed',color='limegreen',linewidth=2, label='molecules')
   ax2 = ax1.twiny()
   ax2.plot(wcm*0.02998, Tmmol,linestyle='dashed',color='limegreen',linewidth=2)
   ax1.axhline(y=Tout,linestyle='dashed',color='black',linewidth=2)
   ax2.set_xlabel(r'$\omega_{\mathrm{cav}}$ $(\mathrm{in}$ $\mathrm{THz})$')
   ax1.text(x=1000,y=Tout-0.02,s=r'$T_{\mathrm{out}}$',fontsize=15)
   ax1.legend(fontsize="14")

   plt.savefig("fig4a.pdf",bbox_inches='tight',dpi=100)
   plt.show()

elif flagi == 5:
   fig1,ax1 = plt.subplots()

   ax1.set_ylabel(r'$T_{\mathrm{set}}$ (in K)',fontsize=25)
   ax1.set_xlabel(r'$\omega_{\mathrm{cav}}$ $(\mathrm{in}$ $\mathrm{cm}^{-1})$', fontsize=25)
   ax1.plot(wcm,Tsetcav,color='limegreen',linewidth=3,label='cavity+molecules')
   ax1.plot(wcm,Tsetmol,linestyle='dashed',color='limegreen',linewidth=3, label='molecules')

   ax2 = ax1.twiny()
   ax2.plot(wcm*0.02998, Tsetmol,linestyle='dashed',color='limegreen',linewidth=3)

   ax2.set_xlabel(r'$\omega_{\mathrm{cav}}$ $(\mathrm{in}$ $\mathrm{THz})$',fontsize=25)
   ax2.set_zorder(-1)

   zoom = 2
   w, h = fig1.get_size_inches()
   fig1.set_size_inches(w * zoom, h * zoom)

   ax1.tick_params(axis='both',labelsize=25)
   ax2.tick_params(axis='both',labelsize=25)

   w, h = fig1.get_size_inches()
   t = fig1.subplotpars.top
   b = fig1.subplotpars.bottom
   axs = h*(t-b)
   l = (1.-axs/w)*0.2/2
   fig1.subplots_adjust(left=l, right=1-l, top=0.8*t, bottom=b)

   ax1.legend(fontsize="25", loc = 'upper left')

   plt.savefig("fig4b.pdf",bbox_inches='tight',dpi=100)
   plt.show()

elif flagi == 6:

   fig1, ax1 = plt.subplots(1)

   ax1.plot(wcm,uTHz,linewidth=2, color = 'turquoise', label = r'$\omega_{\mathrm{cav}}=0.33\mathrm{cm}^{-1}=10\mathrm{GHz}$')
   ax1.plot(wcm,uIR,linewidth=2, color = 'red', label = r'$\omega_{\mathrm{cav}}=500\mathrm{cm}^{-1}=15\mathrm{THz}$')
   ax1.scatter(wcm,uf, color = 'black', label = 'free', marker='o', s=10)

   ax2 = ax1.twiny()
   ax2.plot(wGHz/10**3,uTHz,linewidth=2, color = 'turquoise', label = r'$\omega_{\mathrm{cav}}=0.33\mathrm{cm}^{-1}=10\mathrm{GHz}$')
   ax2.scatter(wGHz/10**3,uf, color = 'black', label = 'free', s=10, zorder=10)

   ax2.set_xlabel(r'$\omega$  $(\mathrm{in}$ $\mathrm{THz})$', fontsize=20)
   ax1.set_xlabel(r'$\omega$  $(\mathrm{in}$ $\mathrm{cm}^{-1})$', fontsize=20)
   ax1.set_ylabel(r'$u(\omega)$  $(\mathrm{in}$ $\mathrm{meV.mm}^{-3}\mathrm{GHz}^{-1})$', fontsize=20)

   ax1.legend(fontsize=14, loc='upper right')
   plt.savefig('fig5a.pdf', bbox_inches='tight')
   plt.show()

   fig1, ax1 = plt.subplots(1)

   ax1.semilogy(wcm,ratiocav2, color='turquoise', linewidth=2, label = r'$\omega_{\mathrm{cav}}=0.33\mathrm{cm}^{-1}=10\mathrm{GHz}$')
   x_left, x_right = ax1.get_xlim()
   y_low2, y_high2 = ax1.get_ylim()
   ax1.semilogy(wcm,ratiocav1, color='red', linewidth=2, label = r'$\omega_{\mathrm{cav}}=500\mathrm{cm}^{-1}=15\mathrm{THz}$')
   y_low1, y_high1 = ax1.get_ylim()
   ax1.set_ylim([y_low1,y_high2])
   ax1.axvline(x=hwcav1,linestyle=(0, (5, 10)),color='black',linewidth=2)
   ax1.text(x=hwcav1+50,y=10,s=r'$\omega=500 \mathrm{cm}^{-1}$',fontsize=15)
   ax1.axvline(x=2*hwcav1,linestyle=(0, (5, 10)),color='black',linewidth=2)
   ax1.text(x=2*hwcav1+50,y=7,s=r'$\omega=2\times 500 \mathrm{cm}^{-1}$',fontsize=15)

   ratio = 5
   ax1.set_aspect(abs((x_right-x_left)/(y_low1-y_high2))*ratio)

   ax1.legend(fontsize=12,loc='center right')
   ax1.set_xlabel(r'$\omega$  $(\mathrm{in}$ $\mathrm{cm}^{-1})$', fontsize=20)
   ax1.set_ylabel(r'$\frac{\rho_{\mathrm{cav}}(\omega)}{\rho_{\mathrm{free}}(\omega)}$', fontsize=20, rotation='horizontal')
   plt.savefig('fig5b.pdf', bbox_inches='tight')
   plt.show()

elif flagi == 7:

   fig1, ax1 = plt.subplots(1)

   hw0 = 1000
   ax1.plot(wc[5], Etot[5], color = '#91bfdb',linewidth=3, label = 'cavity+molecules')
   ax2 = ax1.twiny()
   ax2.plot(wc[5]*0.02998, Etot[8], color = '#fc8d59',linewidth=3)
   ax1.plot(wc[5], Etot[8], color = '#fc8d59', linewidth=3, label = 'molecules')

   ax2.set_zorder(-1)

   ax2.set_xlabel(r'$\omega_{\mathrm{cav}}$ $(\mathrm{in}$ $\mathrm{THz})$',fontsize=20)
   ax1.set_xlabel(r'$\omega_{\mathrm{cav}}$ $(\mathrm{in}$ $\mathrm{cm}^{-1})$', fontsize=20)
   ax1.set_ylabel(r'$\mathrm{Emissivity}$  $\mathcal{E}_{\mathrm{tot}}$', fontsize=20)
   ax1.set_xlim(right=1200)
   ax1.legend(fontsize=17)
#   ax1.legend(fontsize=12,loc = 'upper right')

   plt.savefig('figS1b.pdf', bbox_inches='tight')
   plt.show()
