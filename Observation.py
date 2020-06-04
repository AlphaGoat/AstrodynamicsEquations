import numpy as np
import time
import requests
from timezonefinder import TimezoneFinder
from pytz import timezone, utc
from pytz.exceptions import UnknownTimeZoneError

import Orbit_Builder

# File detailing information of observation site class/ Individual Observations

class Observation_Site():

	def __init__(self,site_name,lat,lon,H,coord_sys,
					obs_equipment,*args,**kwargs):

		self.site_name = site_name
		# lat: latitude of site
		self.lat = lat
		# lon: longitude of site
		self.lon = lon
		# H: elevation of the site above sea level
		self.H = H
		# coord_sys: flag detailing coordinate system being used
		self.coor_sys = coord_sys
		# obs_equipment : list of detection equipment being used
		#                 (e.g., radar array, telescope, etc.)
		self.obs_equipment = obs_equipment

		# Retrieve the timezone of the site as well as utc offset
		self.timezone, self.utc_offset	= self.UTC_Offset()

	def UTC_Offset(self):
		'''Retrieves the time zone and the Universal Time offset
		   for the Observation Site'''

        # IF google timezone API fails, use timezonfinder library instead
		tzf = TimezoneFinder()
		# Determine if import of optimized algorithms was succesful
		if tzf.using_numba() == True:

			# Test timezone functions to see which provides an answer
			# with input longitude/latitude
			if tzf.timezone_at(lng=self.lon,lat=self.lat) != None:

				self.time_zone = tzf.timezone_at(lng=self.lon,lat=self.lat)

			elif tzf.certain_timezone_at(lng=self.lon,
									lat=self.lat) != None:

				self.time_zone = tzf.certain_timezone_at(lng=self.lon,
									lat=self.lat) != None

			else:
				self.time_zone = closest_timezone_at(lng=self.lon,
									lat=self.lat)

			# Retrieve the date now
			time_zone_now = datetime.datetime.now(
										pytz.timezone(self.time_zone)

			# Determine UTC offset based on current date and the
			# time zone of the site
			self.utc_offset	= time_zone_now.utcoffset.totalseconds()/60/60

			return self.timezone, self.utc_offset

		else:
			print('Error: Unable to determine timezone of input
					latitude and longitude')
			return None, None

class Individual_Observation(Observation_Site):

	# Keep dictionary to keep track of current objects beings observed
	current_observations = {}

	def __init__(self,OID,ra,dec,az,el,t,
				 r=None,qoppa=None,SV=None,
				 *args,**kwargs):
		# OID: Object Identification
		self.OID = OID
		# ra: right ascension
		self.ra = ra
		# dec: declination
		self.dec = dec
		# az: azimuth
		self.az = az
		# el: elevation
		self.el = el
		# t: time of observation (in yyyy-dd-mmT00:00)
		self.t = t
		# r: range to target (from Earth center)
		self.r = r
		# qoppa: range to target (from obs site)
		self.qoppa = qoppa
		# SV: State Vector
		self.SV = SV

		# Inherit attributes from parent class (Observation_Site)
		super(Individual_Observation).__init__(lat,lon,coord_sys)

		self.lat = lat
		self.lon = lon
		self.coord_sys = coord_sys

	def Initiate_Orbit_Build(self):
		'''If an object has been spotted for the first time, and has not
		been seen by any other site, build an orbit class for the object.
		Otherwise, contribute information about object gathered by the
		observation into the already constructed orbit'''
		if self.OID not in Orbit.instances:
			current_observations[self.OID] = Orbit(self.OID,num_obs=1)
		else:

	def Local_Sidereal_Time(self):
		'''Calculate local sidereal time, given input longitude'''

		#Determine the number of days between the time of observation and J2000
		time_of_obs = np.datetime64(self.t)
		J2000 = np.datetime64('2000-01-01T02:00')
		diff = time_of_obs - J2000
		diff_days = diff/np.timedelta64(1,'D')

		#Retrieve the current universal time at time of observation
		UT = float(time_of_obs.hour) + float(time_of_obs.minute)/60 + \
						self.utc_offset

		#Calculate and return the local sidereal time
		LST = (100.46 + (0.985647 * diff_days) + self.lon + \
						(15 * UT)) % 360

		if LST < 0: LST = LST + 360

		return LST

	def Construct_Orbit_Info_Dict(self):
		'''Compile information gathered about orbit into dictionary, to
		   be sent to the Orbit instance being built'''


	def RAnDec_to_AzEl(self):
		'''If values for the azimuth and elevation of the observation are
		   not provided, calculate for it from provided ra and dec values'''

		# Retrieve the local sidereal time for the site
		lst = self.Local_Sidereal_Time()
		# Calculate the hour angle (HA), given the provided ra and lst
		ha = lst - self.ra
		# Convert hour angle and latitude into radians for calculations
		ha = ha * (np.pi/180)
		dec = self.dec * (np.pi/180)
		# Calculate elevation
		self.el = np.sin(dec)*np.sin(selflat) + np.cos(dec)*np.cos(ha)* \
							np.cos(self.lat)
		# Calculate azimuth
		cosA = (np.sin(dec) - np.sin(self.el)*np.sin(self.lat)) / (
							np.cos(self.el)*np.cos(self.lat))
		A = np.arccos(A)
		if np.sin(HA) >= 0:
			self.az = A
		else:
			self.az = 360 - A



