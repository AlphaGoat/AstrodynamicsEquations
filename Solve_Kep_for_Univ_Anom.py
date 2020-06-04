import numpy as np

from Universal_Constants import muE

from Astr_Eqtns import Stumpff

#Solves the universal Kepler's Equation for the Universal Anomaly (X)
#given time passed since X = 0, initial position, radial velocity, and the 
#reciprocal of the semi-major axis at X = 0

def Solve_Kep_for_Univ_Anom(del_t,r0,vr0,alpha,eps=1e-8):
	# del_t: time passed since X = 0
	# r0: position when X = 0
	# vr0: radial velocity when X = 0
	# alpha: reciprocal of the semi-major axis at X = 0
	# eps: chosen error tolerance

	#Obtain initial estimate for X
	X = np.sqrt(muE)*np.abs(alpha)*del_t

	#Iterate with Newton's Method to solve for X within specified tolerance

	#Set initial values for F, dF for use in Newton's method
	S = Stumpff(alpha*X^2)
	F = (r0*vr0)/np.sqrt(muE)*X^2*S.StumpffC() + /
			(1-alpha*r0)*X^3*S.StumpffS() + /
			r0*X - np.sqrt(muE)*del_t
	dF = (r0*v0)/np.sqrt(muE)*X*(1 - alpha*X^2*S.StumpffS()) + /
			(1 - alpha*r0)*X^2*S.StumpffC() + r0

	#Iterate to acceptable values for F/dF
	while F/dF > eps:
		X = X - F/dF
		S = Stumpff(alpha*X^2)
		F = (r0*vr0)/np.sqrt(muE)*X^2*S.StumpffC() + /
				(1-alpha*r0)*X^3*S.StumpffS() + /
				r0*X - np.sqrt(muE)*del_t
		dF = (r0*v0)/np.sqrt(muE)*X*(1 - alpha*X^2*S.StumpffS()) + /
				(1 - alpha*r0)*X^2*S.StumpffC() + r0
	
	X = X - F/dF

	return X


