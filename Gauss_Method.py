import numpy as np

import Solve_Kep_for_Univ_Anom

from Astro_Eqtns import Stumpff

from Universal_Constants import muE

#Implements the Gaussian method to angles-only orbit determination
#Details can be found on pg. 236 - 242 of Curtis's "Orbital Mechanics for 
#Engineering Students"

#Finds the orbit for an object after three successive observations at t1,t2,
#and t3

#Need the cosine vectors and the position vectors at the three successive times

def Gaussian_Method(Q,R,t,improved_orbit=False):
	
	# Q: 3x3 matrix consisting of three cosine vectors
	# R: 3x3 matrix consisting of the three position vectors
	# t: 1x3 array consisting of times t1,t2,t3
	# improved_orbit: flag for more accurate orbit estimate	

	#Calculate the time intervals between successive measurements
	tau1 = t[1] - t[2]
	tau3 = t[3] - t[2]
	tau = tau3 - tau1

	#Calculate the cross-products of the cosine vectors
	p1 = np.cross(Q[2],Q[3])
	p2 = np.cross(Q[1],Q[3])
	p3 = np.cross(Q[1],Q[2])

	#Calculate scalar triple product Q1(Q2xQ3) (D0)
	D0 = Q[1].dot(p1)

	#Compute the six scalar quantities
	D11 = np.dot(R[1],p1)
	D12 = np.dot(R[1],p2)
	D13 = np.dot(R[1],p3)
	D21 = np.dot(R[2],p1)
	D22 = np.dot(R[2],p2)
	D23 = np.dot(R[2],p3)
	D31 = np.dot(R[3],p1)
	D32 = np.dot(R[3],p2)
	D33 = np.dot(R[3],p3)

	#Calculate factors A and B
	A = (1/D0)*(-D12*(tau3/tau) + D22 + D32*(tau1/tau))
	B = (1/(6*D0))*(D12*(tau3^2 - tau^2)*(tau3/tau) + D32*(tau^2 -
						tau1^2)*(tau1/tau))

	#Calculate factor E
	E = np.dot(R[2],Q[2])

	#Calculate coefficients to eight order polynomial a,b,c
	a = -(A62 + 2*A*E + R[2]^2)
	b = -2*muE*B*(A + E)
	c = -muE^2*B^2

	#Find the roots of eight order polynomial, and choose the most reasonable
	#answer. This will give r2

	#First, iterate values for r2 near a zero point for the polynomial, 
	#starting at a value of 200km (LEO) and ending at 40000 (GEO)
	r2 = 200
	F = r2^8 + a*r2^6 + b*r2^3 + c
	zeros = []
	i = 0
	for r2 in range(300,40100,100):
		if np.sign(r2^8 + a*r2^6 + b*r2^3 + c) != np.sign(F):
			F = r2^8 + a*r2^6 + b*r2^3 + c
			zeros.append(r2 - 100)
			i += 1
	#Then, use Newton's method to determine the actual zero point near
	#each collected initial estimates (values in the zeros array), within
	#a tolerance of 0.001 km
	epsilon = 0.001
	#list to hold possible answers determined by Newton's Method
	poss_ans = []
	for r2 in zeros:
		F = r2^8 + a*r2^6 + b*r2^3 + c
		dF = 8*r2^7 + 6*a*r2^5 + 3*b*r2^2
		i = 0
		while F/dF < epsilon or i < 100
			r2 = r2 - F/dF
			i += 1
			F = r2^8 + a*r2^6 + b*r2^3 + c
			dF = 8*r2^7 + 6*a*r2^5 + 3*b*r2^2
		poss_ans.append(r2 - F/DF)
	np.asarray(poss_ans)	
	
	#Calculate the slant range at t2
	Q2 = A + (muE*B)/r2^3
	
	#Determine the value of r2 closest to Q2 (most likely value)
	idx = np.abs(poss_ans-Q2).argmin()
	r2  = poss_ans[idx]

	#Calculate the slant ranges at t1, t3
	factor11 = (6*(D31*(tau1/tau3) + D21*(tau/tau3))*r2^3
	factor12 = muE*D31*(tau^2 - tau1^2)*(tau1/tau3)
	factor13 = 6*r2^3 + muE*(tau^2 - tau3^2)
	Q1 = (1/D0)*((factor11 + factor12)/factor13 - D11)

	factor31 = (6*(D13*(tau3/tau1) - D23*(tau/tau1))*r2^3
	factor32 = muE*D13*(tau^2 - tau3^2)*(tau3/tau1)
	factor33 = 6*r2^3 + muE*(tau^2 - tau3^2)
	Q3 = (1/D0)*((factor31 + factor32)/factor33 - D33)

	#Calculate range vectors of observed object from center of Earth
	vec_r1 = R[1] + Q1*Q[1]
	vec_r2 = R[2] + Q2*Q[2]
	vec_r3 = R[3] + Q3*Q[3]

	#Calculate the lagrange coefficients (Note: These are approximations)
	f1 = 1 - (1/2)*(muE/r2^3)*tau1^2
	f3 = 1 - (1/2)*(muE/r2^3)*tau3^2
	g1 = tau1 - (1/6)*(muE/r2^3)*tau1^3
	g3 = tau3 - (1/6)*(muE/r2^3)*tau3^3

	#Calculate the velocity vector at t2
	vec_v2 = (1/(f1*g3 - f3*g1))*(-f3*vec_r1 + f1*vec_r3)

	#Compile vec_r2 and v2 into state vector that can be used to calculate
	#the orbit of the object
	SV = np.array([vec_r2,vec_v2])	

	#If improved_orbit flag is set to true, run state vectors through the
	#improved Gaussian Estimate. Will return a more accurate state vector
	if improved_orbit == True:
		D = np.array([D11,D12,D13],[D21,D22,D23],[D31,D32,D33])
		return Improve_Gauss_Estimate(SV,D0,D)	
	else:
		return SV
				

def Improve_Gauss_Estimate(SV,D0,D,R,Q,i):
	'''Applies iterative improvements to State Vector determined in the
	   previous application of Gauss's Method by finding exact values for
	   lagrange coefficients, as detailed in pg. 243 of Curtis'''
	# SV: State vectors (SV[0] = vec_r2, SV[1] = vec_v2)
	# D0: Scalar triple product Q1(Q2xQ3) 
	# D: Array containing scalar triple products of position and slant vectors
	# R: Array containing position vectors from t1, t2, and t3
	# Q: Array containing slant vectors from t1, t2, t3
	# i: number of times to iterate the formula
	
	#Unpack state vector
	vec_r2 = SV[0]
	vec_v2 = SV[1]	


	#Unpack scalar triple product array
	D11 = D[0][0]
	D12 = D[0][1]
	D13 = D[0][2]
	D21 = D[1][0]
	D22 = D[1][1]
	D23 = D[1][2]
	D31 = D[2][0]
	D32 = D[2][1]
	D33 = D[2][2]

	#Calculate magnitudes of vector r2 and vector v2
	r2 = np.sqrt(vec_r2.dot(vec_r2))
	v2 = np.sqrt(vec_v2.dot(vec_v2.dot(vec_v2))

	#Calculate the reciprocal of the semi-major axis (alpha)
	alpha = 2/r2 - v2^2/muE

	#Calculate the radial components of vec_v2 (vr2)
	vr2 = np.dot(vec_v2,vec_r2/r2)
	
	#Solve for the univseral anomalies X1 and X3 for time intervals
	#tau1 and tau3, respectively
	X1 = Solve_Kep_for_Univ_Anom(tau1,r2,vr2,alpha)
	X3 = Solve_Kep_for_Univ_Anom(tau3,r2,vr2,alpha)

	#Solve for explicit values of f1, f3, g1, and g3, coefficients of lagrange
	#equation
	S1 = Stumpff(alpha*X1^2)
	S3 = Stumpff(alpha*X3^2)

	f1 = 1- (X1^2/r2)*S1.StumpffC()
	g1 = tau1 - (1/np.sqrt(muE))*X1^3*S1.StumpffS()

	f3 = 1 - (X3^2/r2)*S3.StumpffC()
	g3 = tau3 - (1/np.sqrt(muE))*X3^3*S3.StumpffS()

	#Solve for coefficients for vec_r1 and vec_r3 in linear sum of the two
	#to form vec_r2
	c1 = g3/(f1*g3 - f3*g1)
	c3 = -g1/(f1*g3 - f3*g1)

	#Calculate updated values for the slant ranges
	Q1 = (1/D0)*(-D11 + (1/c1)*D21 - (c3/c1)*D31)
	Q2 = (1/D0)*(-c1*D12 + D22 - c3*D32)
	Q3 = (1/D0)*((-c1/c3)*D13 + (1/c3)*D23 - D33)

	#Calculate updated vectors from center of earth to object at t1, t2, and t3
	vec_r1 = R1 + Q1*Q[1]
	vec_r2 = R2 + Q2*Q[2]
	vec_r3 = R3 + Q3*Q[3]

	#Calculate updated velocity vector at t2
	vec_v2 = 1/(f1*g3 - f3*g1)*(-f3*r1 + f1*r3)

	#Iteratively improve the estimate to desired precision, and return final
	#state vector
	SV = np.array([vec_r2,vec_v2])
	if i == 1:
		return SV 
	else:
		Q = np.array([Q1,Q2,Q3])
		Improve_Gauss_Estimate(SV,D0,D,R,Q,i-1)


		

	
	 


