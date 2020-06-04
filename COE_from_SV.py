import numpy as np

from Universal_Constants import muE

# Calculation of the orbital elements from the state vector, as described
# in Curtis's 'Orbital Mechanics for Engineers'


def COE_from_SV(SV):

	#SV: State Vector (SV[1] = R, SV[2] = V)

	#Define the error tolerance epsilon
	eps = 1e-10

	# Unpack the State Vector
	# R -- position vector in the geocentric equatorial frame
	# V -- velocity vector in the geocentric equatorial frame
	R = SV[1]
	V = SV[2]

	# Take the norms of R and V
	r = np.linalg.norm(R)
	v = np.linalg.norm(V)

	# Calculate radial velocity component
	vr = np.dot(R,V)/r

	# Take the cross product of R and V
	H = np.cross(R,V)
	h = np.linalg.norm(H)

	# Calculate the inclination 
	incl = np.arccos(H[3]/h)

	# Calculate vector N defining node line (pg 160)
	N = np.cross(np.array([0, 0, 1]),H)
	n = np.linalg.norm(N)

	# Calculate RA of the ascending node
	if n ~= 0:
		RA = np.arccos(N[1]/n)
		if N[2] < 0:
			RA = 2*np.pi - RA
	else:
		RA = 0

	# Calculate the eccentricity vector
	E = 1/muE*((v^2-muE/r)*R - r*vr*V)
	e = np.linalg.norm(E)

	# Calculate the true anomaly (TA)
	if e > eps:
		TA = np.arccos(np.dot(E,R)/e/r)
		if vr < 0:
			TA = 2*np.pi - TA
	else:
		cp = np.cross(N,R)
		if cp[3] >= 0:
			TA = np.arccos(np.dot(N,R)/n/r)
		else:
			TA = 2*np.pi - arccos(np.dot(N,R)/n/r)

	# Calulate the semimajor axis of the ellipse
	a = h^2/muE/(1 - e^2)

	return h,e,RA,incl,w,TA,a
	


