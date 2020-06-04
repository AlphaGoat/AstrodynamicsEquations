import numpy as np

from Universal_Constants import oblate, Re, wE

# Determines the state vectors of an object based on range and 
# angle measurements
# As detailed in pg 228-232 of Curtis's 'Orbital Mechanics for Engineers'

def StateVecs_from_AngRange(H,lat,lst,qoppa,A,a,dt_qoppa,dt_A,dt_a):
	# H: elevation of the observing site
	# lat: latitude of the observing site
	# lst: local sidereal time of the observing time
	# qoppa: range to the target [m]
	# A: azimuth to target [rad]
	# a: angular elevation to target

	#Calculate geocentric position vector (R)
	factor = Re/np.sqrt(1-(2*oblate-oblate^2)*np.sin(lat)^2) + H
	R = np.array([factor*np.cos(lat)*np.cos(lst),
				  factor*np.cos(lat)*np.sin(lst),
				  factor*np.sin(lst)])

	#Calculate the topocentric declination
	dec = np.arcsin(np.cos(lat)*np.cos(A)*np.cos(a)+np.sin(lat)*np.sin(a))
	
	#Calculate the topocentric right ascension (ra)
	if A > 0 and A < np.pi:
		h = 2*pi - np.arccos((np.cos(lat)*np.sin(a)-\
								np.sin(lat)*np.cos(A)*np.cos(a))/np.cos(dec))
	else:
		h = np.arccos((np.cos(lat)*np.sin(a)-\
								np.sin(lat)*np.cos(A)*np.cos(a))/np.cos(dec))
	ra = lst - h

	#Calculate the direction cosine unit vector 
	del_qoppa = np.array([np.cos(dec)*np.cos(ra),
						  np.cos(dec)*np.sin(ra),
						  np.sin(dec)])

	#Calculate the geocentric position vector
	r = R + qoppa*del_qoppa

	#Calculate the inertial velocity of the observing site (dt_R)
	sigma = np.array([0,0,wE])
	dt_R = 	np.cross(sigma,R)

	#Calculate the declination rate (dt_dec)
	factor1 = -dt_A*np.cos(lat)*np.sin(A)*np.cos(a)
	factor2 = dt_a*(np.sin(lat)*np.cos(a) - np.cos(lat)*np.cos(A)*np.sin(a))
	dt_dec = (1/np.cos(dec))*(factor1 + factor2)

	#Calculate the right ascension rate (dt_ra)
	factor3 = dt_A*np.cos(A)*np.cos(a) - dt_a*np.sin(A)*np.sin(a) + \
					dt_dec*np.sin(A)*np.cos(a)*np.tan(dec)
	factor4 = np.cos(lat)*np.sin(a) - np.sin(lat)*np.cos(A)*np.cos(a)
	dt_ra = wE + (factor3/factor4)

	#Calculate the direction cosine rate vector (dt_del_qoppa)
	Ifactor = -dt_ra*np.sin(ra)*np.cos(dec) - dt_dec*np.cos(ra)*np.sin(dec)
	Jfactor = dt_ra*np.cos(ra)*np.cos(dec) - dt_dec*np.sin(ra)*np.sin(dec)
	Kfactor = dt_dec*np.cos(dec)
	dt_del_qoppa = np.array([Ifactor,Jfactor,Kfactor])

	#Calculate the geocentric velocity vector 
	v = dt_R + dt_qoppa*del_qoppa + qoppa*dt_del_qoppa

	#Return State Vector
	SV = np.array([r,v])

	return SV

	
	

