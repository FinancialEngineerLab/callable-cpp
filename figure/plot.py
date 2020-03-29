import numpy as np
import matplotlib.pyplot as plt

expiries = [0.25, 0.5, 1, 2, 3, 4, 5, 7, 10,  15 , 20,  25,  30  ]
vols1 = [0.298953,  0.299422,  0.299633,  0.299836,  0.299782,  0.299858,  0.299912,  0.29982,   0.299933,  0.299917,  0.299818,  0.299932,  0.299837  ]
vols2 = [0.315094,  0.320979,  0.327881,  0.336034,  0.34115,   0.345087,  0.348207,  0.352728,  0.357405,  0.361271,  0.362436,  0.36245,   0.361504 ]
vols3 = [0.32609,   0.335716,  0.347209,  0.360639,  0.368946,  0.375023,  0.379608,  0.385731,  0.390817,  0.393045,  0.391759,  0.389326,  0.386182 ]
vols4 = [0.354366,  0.373722,  0.396909,  0.422808,  0.437504,  0.446862,  0.452849,  0.458632,  0.45976,   0.453763,  0.445212,  0.436815,  0.428868  ]

plt.figure(1)
plt.title("Panel A, p = 2")
plt.plot(expiries, vols1, 'k-',  label = 'c=0')
plt.plot(expiries, vols2, 'ob-', label = 'c=3%')
plt.plot(expiries, vols3, '^g-', label = 'c=5%')
plt.plot(expiries, vols4, '+r-', label = 'c=10%')

legend = plt.legend(loc='lower center', shadow=True, fontsize='small')
plt.ylim(.2, .5)
plt.xlim(0., 30.)
plt.xlabel('Option maturity, T (year)')
plt.ylabel('Implied volatility')

plt.grid()
plt.savefig("fig_2_panel_A.png")

plt.figure(2)
plt.title("Panel B, c = 5%")
vols1 = [0.329346, 0.341793, 0.358261, 0.380009, 0.395328, 0.407568, 0.417755, 0.434035, 0.45266, 0.474817, 0.491001, 0.504062, 0.515175]
vols2 = [0.328407, 0.339939, 0.354624, 0.372966, 0.385055, 0.39422, 0.40146, 0.412156, 0.423045, 0.433453, 0.438816, 0.441607, 0.442584]
vols3 = [0.32609, 0.335716, 0.347209, 0.360639, 0.368946, 0.375023, 0.379608, 0.385731, 0.390817, 0.393045, 0.391759, 0.389326, 0.386182]

plt.plot(expiries, vols1, '+r-', label = 'p=0')
plt.plot(expiries, vols2, 'ob-', label = 'p=1/2')
plt.plot(expiries, vols3, '^g-', label = 'p=2')

legend = plt.legend(loc='lower center', shadow=True, fontsize='small')
plt.ylim(.25, .55)
plt.xlim(0., 30.)
plt.xlabel('Option maturity, T (year)')
plt.ylabel('Implied volatility')

plt.grid()
plt.savefig("fig_2_panel_B.png")
