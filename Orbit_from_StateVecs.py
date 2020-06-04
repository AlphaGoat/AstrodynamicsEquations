import numpy as np

from Universal_Constants import muE

#Determines the orbital elements of an object from the state vector,
#as detailed in algorithm 4.1 of Curtis's 'Orbital Mechanics for Engineers',
#pg. 159-161

def Orbit_from_StateVecs(r,v):
	# r: position vector
	# v: velocity vector

	#Calculate the distance
	distance = np.sqrt(r.dot(r))

	#Calculate the speed	
	speed = np.sqrt(v.dot(v))

	#Calculate the radial velocity
	rad_velocity = np.dot(r,v)/distance

	#Calculate the specific angular momentum
	ang_moment = np.cross(r,v)

	#Calculate the magnitude of the specific angular momentum
	mag_ang_moment = np.sqrt(ang_moment.dot(ang_moment))

	#Calculate the inclination
	inc = np.arccos(ang_moment[3]/mag_ang_moment)

	#Calculate vector for the node line
	N = np.cross([0,0,1],ang_moment)

	#Calculate the magnitude of N
	mag_N = n.sqrt(N.dot(N))

	#Calculate the right ascension of the ascending node (RA)
	if N[2] >= 0:
		RA = np.arccos(N[1]/mag_N)
	else:
		RA = 2*np.pi - np.arccos(N[1]/mag_N) 

	#Calculate the eccentricity vector (e_vec)
	e_vec = (1/muE)*((speed^2-(muE/distance))*r - distance*rad_velocity*v)

	#Calculate the eccentricity (e)
	e = np.sqrt(e_vec.dot(e_vec))

	#Calculate the argument of perigee (omega)
	if e_vec[3] >= 0:
		omega = arccos(np.dot(N,e_vec)/(mag_N*e))
	else:
		omega = 2*np.pi - arccos(np.dot(N,e_vec)/(mag_N*e))

	#Calculate the True Anomaly (ta)
	if rad_velocity >= 0:
		ta = np.arccos(np.dot(e_vec,r)/(e*distance))
	else:
		ta = 360 - np.arccos(np.dot(e_vec,r)/(e*distance))

	return mag_ang_moment,inc,RA,e,omega,ta

