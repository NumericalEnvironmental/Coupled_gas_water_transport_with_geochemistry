# define any up-front keyword blocks for PHREEQC input

input_string = """
PHASES
CO2(g)
	CO2 = CO2
	-log_k	-1.468
	-delta_h -4.776 kcal
	-analytic   10.5624  -2.3547e-2  -3972.8  0  5.8746e5  1.9194e-5
Ntg(g)
	Ntg = Ntg
	-analytic -58.453 1.81800e-3  3199  17.909 -27460
"""

