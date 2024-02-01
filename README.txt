########################################################################################
Code accompanying manuscript "Blackbody radiation and thermal effects on chemical 
reactions and phase transitions in cavities" (https://doi.org/...)

The files CaF2.txt contains refractive index of CaF2 taken from https://doi.org/10.1016/j.optmat.2017.06.052 
for w<=900 cm^{-1}, CaF2k.txt and CaF2n.txt contain refractive indices taken from 
https://doi.org/10.1177/0003702816657568 for w>900 cm^{-1}.

a) The code 'main3.py' with:
    1) flagi = 0,1,2,3,4,5,7,8 and flagT=0 (T=298K) generates 'LzEtot_<flagi>.txt' and 
       'fig3_X.txt', 'fig3_Y.txt', 'fig3_E_*_<flagi>.txt' and 'fig3_Ts_*_<flagi>.txt'.
       used in Fig. 3, Fig. 4a and Fig. S1a.
    1) flagi = 5,8 and flagT=1 (T=350K) generates 'LzEtot_<flagi>.txt' and 
       'fig3_X.txt', 'fig3_Y.txt', 'fig3_E_*_<flagi>.txt' and 'fig3_Ts_*_<flagi>.txt'.
       used in Fig. 4b and Fig.S1b.

b) The code 'main4a.py' uses 'LzEtot_1.txt' and 'LzEtot_4.txt' at 298K to  generate 'kmol.txt' and 
   'kcav.txt' for Fig. 4a

c) The code 'main4b.py' uses 'LzEtot_5.txt' and 'LzEtot_8.txt' at 350K to generate 'Tset_mol_cav.txt' 
   for Fig. 4b

d) The code 'main5a.py' generates 'fig5a.txt' for Fig. 5a
e) The code 'main5b.py' generates 'fig5b.txt' for Fig. 5b

f) The code 'plot.py' with:
    1) flagi = 1 generates Fig. 1b
    2) flagi = 3 uses 
           a) 'LzEtot_5.txt', 'LzEtot_7.txt' and 'LzEtot_8.txt' at 298K to generate the Etot plot 
              'fig3a.pdf' of Fig. 3a
           b) 'fig3_E_333_5.txt', 'fig3_X.txt' and 'fig3_Y.txt' at 298K to generate the E(k) inset 
              plot 'fig3a1.pdf' of Fig. 3a
           c) 'fig3_Ts_333_5.txt', 'fig3_X.txt' and 'fig3_Y.txt' at 298K to generate the T(\theta=0) 
              inset plot 'fig3a2.pdf of Fig. 3a
           d) 'LzEtot_0.txt', 'LzEtot_2.txt' and 'LzEtot_3.txt' at 298K to generate the Etot plot 
              'fig3b.pdf' of Fig. 3b 
           e) 'fig3_E_500_0.txt', 'fig3_X.txt' and 'fig3_Y.txt' at 298K to generate the E(k) inset 
              plot 'fig3b1.pdf' of Fig. 3b
           f) 'fig3_Ts_500_0.txt', 'fig3_X.txt' and 'fig3_Y.txt' at 298K to generate the T(\theta=0) 
              inset plot 'fig3b2.pdf' of Fig. 3b
           g) 'LzEtot_1.txt' and 'LzEtot_4.txt' at 298K to generate the Etot plot 'figS1a.pdf' of 
              Fig. S1a
    3) flagi = 4 uses 'kmol.txt' and 'kcav.txt' to generate 'fig4a.pdf' - Fig. 4a
    4) flagi = 5 uses 'Tset_mol_cav.txt' to  generate 'fig4b.pdf' - Fig. 4b
    5) flagi = 6 uses 'fig5a.txt' and 'fig5b.txt' to  generate 'fig5a.pdf' and 'fig5b.pdf',
       respectively. - Fig. 5a and 5b
    6) flagi = 7 uses 'LzEtot_5.txt' and 'LzEtot_8.txt' at T=350K to generate 'figS1b.pdf'

#######################################################################################
