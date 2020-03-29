import numpy as np
import matplotlib.pyplot as plt

plt.figure(1)
plt.title("Panel A, p = 2")

expiries = [0.25, 0.5, 1, 2, 3, 4, 5, 7, 10,  15 , 20,  25,  30  ]
vols1 = [0.298953,  0.299422,  0.299633,  0.299836,  0.299782,  0.299858,  0.299912,  0.29982,   0.299933,  0.299917,  0.299818,  0.299932,  0.299837  ]
vols2 = [0.315094,  0.320979,  0.327881,  0.336034,  0.34115,   0.345087,  0.348207,  0.352728,  0.357405,  0.361271,  0.362436,  0.36245,   0.361504 ]
vols3 = [0.32609,   0.335716,  0.347209,  0.360639,  0.368946,  0.375023,  0.379608,  0.385731,  0.390817,  0.393045,  0.391759,  0.389326,  0.386182 ]
vols4 = [0.354366,  0.373722,  0.396909,  0.422808,  0.437504,  0.446862,  0.452849,  0.458632,  0.45976,   0.453763,  0.445212,  0.436815,  0.428868  ]

plt.plot(expiries, vols1, 'k-',  label = 'c=0')
plt.plot(expiries, vols2, 'ob-', label = 'c=3%')
plt.plot(expiries, vols3, '^g-', label = 'c=5%')
plt.plot(expiries, vols4, '+r-', label = 'c=10%')

legend = plt.legend(loc='lower center', shadow=True, fontsize='small')
plt.ylim(.2, .5)
plt.xlim(0., 30.)
plt.xlabel('Option maturity, T (years)')
plt.ylabel('Implied volatility')

plt.grid()
plt.savefig("fig_2_panel_A.png")

plt.figure(2)
plt.title("Panel B, c = 5%")

expiries = [0.25, 0.5, 1, 2, 3, 4, 5, 7, 10,  15 , 20,  25,  30  ]
vols1 = [0.329346, 0.341793, 0.358261, 0.380009, 0.395328, 0.407568, 0.417755, 0.434035, 0.45266, 0.474817, 0.491001, 0.504062, 0.515175]
vols2 = [0.328407, 0.339939, 0.354624, 0.372966, 0.385055, 0.39422, 0.40146, 0.412156, 0.423045, 0.433453, 0.438816, 0.441607, 0.442584]
vols3 = [0.32609, 0.335716, 0.347209, 0.360639, 0.368946, 0.375023, 0.379608, 0.385731, 0.390817, 0.393045, 0.391759, 0.389326, 0.386182]

plt.plot(expiries, vols1, '+r-', label = 'p=0')
plt.plot(expiries, vols2, 'ob-', label = 'p=1/2')
plt.plot(expiries, vols3, '^g-', label = 'p=2')

legend = plt.legend(loc='lower center', shadow=True, fontsize='small')
plt.ylim(.25, .55)
plt.xlim(0., 30.)
plt.xlabel('Option maturity, T (years)')
plt.ylabel('Implied volatility')

plt.grid()
plt.savefig("fig_2_panel_B.png")

plt.figure(3)
plt.title("Panel A, T = 0.5")

moneynesses = [0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5 ]
vols1 = [0.466032, 0.385233, 0.341766, 0.322486, 0.313718, 0.309462, 0.307205, 0.306164, 0.305627, 0.305543 ]
vols2 = [0.574747, 0.469502, 0.397063, 0.356571, 0.335716, 0.324748, 0.318546, 0.315024, 0.312824, 0.311571 ]
vols3 = [0.699762, 0.576714, 0.478693, 0.412741, 0.373722, 0.351497, 0.338386, 0.330459, 0.325301, 0.321973 ]

plt.plot(moneynesses, vols1, '+k-', label = 'c=2%')
plt.plot(moneynesses, vols2, 'ob-', label = 'c=5%')
plt.plot(moneynesses, vols3, '^g-', label = 'c=10%')

legend = plt.legend(loc='upper right', shadow=True, fontsize='small')
plt.ylim(.25, .75)
plt.xlim(0.55, 1.55)
plt.xlabel('Strike (% of forward)')
plt.ylabel('Implied volatility')

plt.grid()
plt.savefig("fig_3_panel_A.png")

plt.figure(4)
plt.title("Panel B, T = 5")

moneynesses = [0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5 ]
vols1 = [0.375463, 0.358324, 0.346468, 0.338187, 0.332204, 0.327587, 0.324023, 0.321261, 0.319123, 0.317383 ]
vols2 = [0.466252, 0.434929, 0.411173, 0.393263, 0.379608, 0.368874, 0.360397, 0.353651, 0.348253, 0.343809 ]
vols3 = [0.577991, 0.537242, 0.503293, 0.475485, 0.452849, 0.434264, 0.419029, 0.406516, 0.396208, 0.387578 ]

plt.plot(moneynesses, vols1, '+k-', label = 'c=2%')
plt.plot(moneynesses, vols2, 'ob-', label = 'c=5%')
plt.plot(moneynesses, vols3, '^g-', label = 'c=10%')

legend = plt.legend(loc='upper right', shadow=True, fontsize='small')
plt.ylim(.25, .6)
plt.xlim(0.55, 1.55)
plt.xlabel('Strike (% of forward)')
plt.ylabel('Implied volatility')

plt.grid()
plt.savefig("fig_3_panel_B.png")

plt.figure(5)

moneynesses = [0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5 ]
vols1 = [0.508949, 0.444747, 0.402117, 0.375254, 0.358261, 0.347307, 0.339799, 0.334552, 0.330683, 0.327848 ]
vols2 = [0.511976, 0.442016, 0.394826, 0.365368, 0.347209, 0.335873, 0.328354, 0.323268, 0.319629, 0.317051 ]

plt.plot(moneynesses, vols1, '+k-', label = 'p=0')
plt.plot(moneynesses, vols2, 'ob-', label = 'p=2')

legend = plt.legend(loc='upper right', shadow=True, fontsize='small')
plt.ylim(.3, .55)
plt.xlim(0.55, 1.55)
plt.xlabel('Strike (% of forward)')
plt.ylabel('Implied volatility')

plt.grid()
plt.savefig("fig_4.png")

plt.figure(6)
plt.title("Panel A, p = 2")

expiries = [0.25, 0.5, 1, 2, 3, 4, 5, 7, 10,  15 , 20,  25,  30  ]
vols1 = [0.0204762, 0.0209529, 0.0218993, 0.0237012, 0.0252822, 0.0265834, 0.0276058, 0.0289493, 0.0297788, 0.0296364, 0.0287963, 0.0277764, 0.0267505 ]
vols2 = [0.0507761, 0.0514984, 0.0527421, 0.0543924, 0.055074, 0.0550737, 0.0546402, 0.0531043, 0.050268, 0.0457708, 0.0420359, 0.0389844, 0.0364622 ]
vols3 = [0.100192, 0.100218, 0.0996805, 0.0969985, 0.0933325, 0.0894644, 0.0857223, 0.079002, 0.0708598, 0.0610573, 0.0541758, 0.0490496, 0.0450583 ]

plt.plot(expiries, vols1, 'ob-', label = 'c=2%')
plt.plot(expiries, vols2, '^g-', label = 'c=5%')
plt.plot(expiries, vols3, '+r-', label = 'c=10%')

legend = plt.legend(loc='lower center', shadow=True, fontsize='small')
plt.ylim(0, .12)
plt.xlim(0., 30.)
plt.xlabel('Bond maturity, T (years)')
plt.ylabel('Risky spread, s')

plt.grid()
plt.savefig("fig_5_panel_A.png")

plt.figure(7)
plt.title("Panel B, c = 5%")

expiries = [0.25, 0.5, 1, 2, 3, 4, 5, 7, 10,  15 , 20,  25,  30  ]
vols1 = [0.0499723, 0.0499728, 0.0499729, 0.0499729, 0.0499729, 0.0499729, 0.0499729, 0.0499729, 0.0499729, 0.0499728, 0.0499728, 0.0499727, 0.0499726 ]
vols2 = [0.0507761, 0.0514984, 0.0527421, 0.0543924, 0.055074, 0.0550737, 0.0546402, 0.0531043, 0.050268, 0.0457708, 0.0420359, 0.0389844, 0.0364622 ]
vols3 = [0.0520171, 0.0538322, 0.0567695, 0.0598742, 0.0604903, 0.0598745, 0.0587001, 0.0558203, 0.0515766, 0.0458271, 0.0415035, 0.0381584, 0.0354868 ]

plt.plot(expiries, vols1, '+r-', label = 'p=0')
plt.plot(expiries, vols2, 'ob-', label = 'p=2')
plt.plot(expiries, vols3, '^g-', label = 'p=3')

legend = plt.legend(loc='lower center', shadow=True, fontsize='small')
plt.ylim(.03, .07)
plt.xlim(0., 30.)
plt.xlabel('Bond maturity, T (years)')
plt.ylabel('Risky spread, s')

plt.grid()
plt.savefig("fig_5_panel_B.png")

plt.figure(8)

expiries = [0.25, 0.5, 1, 2, 3, 4, 5, 7, 10,  15 , 20,  25,  30  ]
vols1 = [0.0977648, 0.0957332, 0.0918934, 0.0850622, 0.0792179, 0.0741862, 0.0698168, 0.0626059, 0.0544538, 0.0450907, 0.0387056, 0.0340313, 0.0304399 ]
vols2 = [0.100192, 0.100218, 0.0996805, 0.0969985, 0.0933325, 0.0894644, 0.0857223, 0.079002, 0.0708598, 0.0610573, 0.0541758, 0.0490496, 0.0450583 ]
vols3 = [0.106237, 0.111713, 0.119969, 0.127103, 0.127481, 0.125234, 0.122057, 0.115254, 0.106155, 0.0946601, 0.0863993, 0.080171, 0.075281 ]

plt.plot(expiries, vols1, '+r-', label = 'sigma=15%')
plt.plot(expiries, vols2, 'ob-', label = 'sigma=30%')
plt.plot(expiries, vols3, '^g-', label = 'sigma=50%')

legend = plt.legend(loc='lower center', shadow=True, fontsize='small')
plt.ylim(.0, .14)
plt.xlim(0., 30.)
plt.xlabel('Bond maturity, T (years)')
plt.ylabel('Risky spread, s')

plt.grid()
plt.savefig("fig_6.png")
