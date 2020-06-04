import numpy as np

#Importing Universal Constants
#oblate: oblateness of the Earth
#Re: radius of the Earth
from Univ_Consts import oblate, Re

class TEC(object):
	'''Implements the Topocentric Equatorial Coordinate System'''
	def __init__(self,r,lst,lon,lat,H=0,qoppa):
		# R: GEC position vector of observed object
		self.r = r	
		# lst: local sidereal time
		self.lst = lst
		# lon: longitude of observing station
		self.lon = lon
		# lat: latitude of observing station
		self.lat = lat
		# H: elevation of observing station above ellipsoid surface
		self.H = H

	def position_vector(self):
		'''Determines position vector R for the observer'''
		factor = (Re/np.sqrt(1-(2*oblate-oblate^2)*np.sin(self.lat)^2)
					+self.H)*np.cos(self.lat)
		R = np.array([factor*np.cos(self.lat)*np.cos(self.lon),
					  factor*np.cos(self.lat)*np.sin(self.lon),
					  factor*np.sin(self.lat)])

		return R
 
	def qoppa_vector(self):
		'''Determines the relative position vector qoppa, aka the relative
		   distance between the observer and the observed object'''
		qoppa = self.r - self.position_vector()

		return qoppa

	def topocentric_declension(self):
		'''Determines the Topocentric Declension from the calculated qoppa
		   vector'''
		qoppa = qoppa_vector()
		mag_qoppa = np.sqrt(qoppa.dot(qoppa))
		unit_qoppa = qoppa/mag_qoppa
		declension = np.arcsin(unit_qoppa[3])

		return declension

	def topocentric_right_ascension(self):
		'''Determines the Topocentric Right Ascension from the calculated 
		   qoppa vector and the declension'''
		qoppa = qoppa_vector()
		declension = topocentric_declension()
		#determine whether right ascension lies between 0-pi or pi-2pi
		if qoppa[2]/np.cos(declension) < 0:
			return 	np.pi/2 - np.arcsin(qoppa[2]/np.cos(declension))
		else 
			return np.arccos(qoppa[1]/np.cos(declension))	

	
