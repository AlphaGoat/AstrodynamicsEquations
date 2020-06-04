import numpy as np
import Geocentric_Equatorial_Coordinates as GEC #not implemented yet

# Class implementing Topocentric Horizon Coordinate System for Angles-Only
# Orbit determination

# See Pages 223-228 from Curtis's Orbital Mechancis for Engineering Students

deg2rad = np.pi/180

class THC(object):
	'''Implements the Topocentric Horizon Coordinate System'''
	def __init__(self,q,a,Az,theta,psi,deg_rad_flag):

		# q: distance from the observer to the observed object
		# a: elevation
		# Az: azimuth
		# theta: angle between Greenwich Meridian and longitude of site
		# psi: angle between Equator and Latitude of site
		# deg_rad_flag: flag detailing whether psi and theta were input as
		#               radians or degrees (based on input lat/lon) 

		#Distance from the observer to the observed object
		self.q = q 

		# Relative Position Vector (qoppa):
		self.a = a
		self.az = az
		self.qoppa = np.array([np.sin(Az)*np.cos(a),
							   np.cos(Az)*np.cos(a),
							   np.cos(a)])
		
		#Angle between Greenwich Meridian and the longitude of the site
		if deg_rad_flag = 'degrees':
			self.theta = theta * deg2rad
		elif deg_rad_flag = 'radians':
			self.theta = theta
		else:
			Astro_Error.Flag_Error('Error: deg_rad_flag not 
									initialized properly')

		#Angle between Equator and Latitude of the site
		if deg_rad_flag = 'degrees':
			self.psi = psi * deg2rad
		elif deg_rad_flag = 'radians':
			self.psi = psi
	
		#Defining Transform matrix converting GEC to topocentric base vectors 	
		Q_Xx = np.array([-np.sin(self.theta),np.cos(self.theta),0],
					    [-np.sin(self.psi)*np.cos(self.theta),
						 -np.sin(self.psi)*np.sin(self.theta),np.cos(self.psi)],
					    [np.cos(self.psi)*np.cos(self.theta),np.cos(self.psi)*/
					     np.sin(self.theta),np.sin(self.psi)])

		#Reverse Transform of Q, converting THC to GEC
		Q_xX = np.array([-np.sin(self.theta),-np.sin(self.psi)*/
						 np.cos(self.theta),np.cos(self.psi)*/
						 np.cos(self.theta)],
						[np.cos(self.theta),-np.sin(self.psi)*/
						 np.sin(self.theta),np.cos(self.psi)*/
						 np.sin(self.theta)]
						[0,np.cos(self.psi),np.sin(self.psi)])

		#Defining topocentric base vectors ijk from IJK base vectors of 
		#geocentric equatorial frame	 	
		self.base = np.dot(Q,GEC.base)

		
	def convert_to_TEC(self,v)
		'''Returns vector converted into Topocentric Equatorial Coordinates''' 	
		
		return np.dot(self.Q_xX,v)
	
			
		
