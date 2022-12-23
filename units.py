import numpy as np

"""
NEURON:

Time		t 	[ms]
Voltage 	v 	[mV]
Current 	i 	[mA/cm2] (distributed)
				[nA] (point process)

Concentration 	ko, ki, etc. 	[mM]

Specific capacitance 	cm 	[uf/cm2]

Length 	diam, L 	[um]

Conductance 	g 	[S/cm2] (distributed)
					[uS] (point process)

Cytoplasmic resistivity 	Ra 	[ohm cm]

Resistance 	Ri( ) 			[10**6 ohm]
"""


"""
Sterratt's et al book:
Table 2.3 p. 34


d 	Diameter of neurite 	[μm]
l 	Length of compartment 	[μm]
R_m 	Specific membrane resistance [ohm cm2]
C_m 	Specific membrane capacitance [μF cm−2]
R_a 	Specific axial resistance (resistivity) [ohm cm]
r_m 	Membrane resistance per inverse unit length [ohm cm] r_m = R_m/(πd)
c_m 	Membrane capacitance per unit length	[μF cm−1] c_m = C_m πd
r_a		Axial resistance per unit length [ohm/ cm] r_a = 4R_a/(πd**2)
V 		Membrane potential [mV]
E_m 	Leakage reversal potential due to different ions [mV]
I 		Membrane current density [μA cm−2]
I_e Injected current [nA]
I_c Capacitive current [nA]
I_i Ionic current [nA]


fig 2.12 p. 30
Conductance 	g 	[mS/cm2] 
					
"""





def Sterratt_2NRN():



def NRN_2Sterratt():
	t = 