import numpy as np

from Universal_Constants import oblate, Re

# File defining the differnet coordinate systems used to determine orbit, as
# well as the transforms between them

class Coordinates(object):

	def __init__(self,*args,**kwargs):
	

class TEC(Coordinates):
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
            return  np.pi/2 - np.arcsin(qoppa[2]/np.cos(declension))
        else
            return np.arccos(qoppa[1]/np.cos(declension))



class THC(Coordinates):
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
                                                                                		#Defining topocentric base vectors ijk from IJK base vectors of 
		#geocentric equatorial frame
		self.base = np.dot(Q,GEC.base)

	def convert_to_TEC(self,v)
		'''Returns vector converted into Topocentric Equatorial Coordinates'''

		return np.dot(self.!_xX,vd)
