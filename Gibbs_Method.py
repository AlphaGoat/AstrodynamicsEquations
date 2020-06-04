import numpy as np

from Universal_Constants import muE

# Gibbs Method of Preliminary Orbit Determination, as derived in Curtis's 
# 'Orbital Mechanics for Engineering Students' pgs. 194-198

# Solves for the orbit of an object based on three successive observations

# R1, R2, R3 are the position vectors for the object at t1, t2, t3, respectively# R1, R2, and R3 are represented by numpy arrays 


def Gibbs_Solver(r1_vec,r2_vec,r3_vec):

	# r_vec: position vectors from center of Earth to observed object

	# Define tolerance for determination of coplanar vectors

	# Take the norm of the three position vectors
	r1_mag = np.linalg.norm(r1_vec)
	r2_mag = np.linalg.norm(r2_vec)
	r3_mag = np.linalg.norm(r3_vec)

	# Take the cross products of the three position vectors
	c12 = np.cross(r1_vec,r2_vec)
	c23 = np.cross(r2_vec,r3_vec)
	c31 = np.cross(r3_vec,r1_vec)
	
	# Check if R1, R2, and R3 are coplanar; if not, set error flag
	if np.absolute(np.dot(r1_vec,c23)/r1_mag/np.linalg.norm(c23)) > tol
		ierr = 1

	# Relational variables N,D, and S, as used in the Gibbs derivation used in 
	# Curtis's text
	N = r1_mag*c23 + r2_mag*c31 + r3_mag*c12
	D = c12 + c23 +c31
	S = r1_vec*(r2_vec-r3_vec) + R2(r3_vec-r1_vec) + R3(r1_vec-r2_vec)

	# Take the norms of N and D
	n = np.linalg.norm(N)
	d = np.linalg.norm(D)

	# Solve for the velocity at position 2
	v2 = np.sqrt(muE/(n*d))*(np.cross(D,r2_vec)/r2_mag + S)

	#Return state vector at position 2
	SV = np.array([r2_vec,v2])

	return ie, SV	
	


