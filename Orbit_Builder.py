import numpy as np


class Orbit():
	'''Orbit class for keeping track of information needed to build an orbit,
		ultimately determining the six orbital elements once the orbit has
		been built'''
	#Instances: Dictionary used to keep track of (a) how many instances of the
	#           orbit class have been created and (b) if the instance for a 
	#           particular observed object has been created or not
	instances = []

	def __init__(self,OID,H=None,i=None,ra=None,
					  e=None,omega=None,ta=None,
					  obs=None,num_obs=0,
				      *args,**kwargs):

		# H: moment of angular momentum
		self.H = H
		# i: inclination
		self.i = i
		# ra: right ascension
		self.ra = ra
		# e: eccentricity
		self.e = e
		# omega: argument of perigee
		self.omega = omega
		# ta: true anomaly
		self.ta = ta
		# obs: dict containing information on number of observation, info
		#      gathered during obs, observation sites, etc. First thing looked
		#      at when building orbit
		self.obs = obs
		# num_obs: number of observations performed
		self.num_obs = num_obs

		# Load the instances list with object identification
		Orbit.instances.append(self)

		# initialize dictionary containing information of each observation
		# (keyword for each site)
		self.orbit_info = {}
		# initialize dictionary containing number each type of observation
		# occurs 
		self.orbit_numbers = {}
		
	def Unload_Obs(self,obs_dict):
		'''Every time a new observation occurs, unload observational 
		   data into appropiate dictionary'''
		# Retrieve site name and initialize an empty dicitonary tied to
		# that key to load the rest of the data
		try:
			site_name = obs_dict['site name']
			# Check if earlier observations from site have been loaded into
			# the dictionary. If so, add another key specifying the number
			# of the observation
			if site_name in self.orbit_info:
				obs_index = len(self.orbit_info[site_name])
			else:	
				self.orbit_info[site_name] = {}
				obs_index = 0
			self.orbit_info[site_name][obs_index] = {}
			# coordinates is simply an identifier, ignore 
			coordinates = None
		except KeyError:
			# If no site name is provided (shame on the user!!!!), use
			# latitude/longitude coordinates instead
			try:
				site_name = None
				longitude = str(obs_dict['longitude'])
				latitude = str(obs_dict['latitude'])
				coordinates = latitude+','+longitude]
				# Check if the coordinates provided have been loaded into 
				# the dictionary. If so, count the number of obs and make 
				# this entry +=1
				if coordinates in self.orbit_info:
					obs_index = len(self.orbit_info[coordinates])
				else:
					self.orbit_info[coordinates] = {}
					obs_index = 0
				self.orbit_info[coordinates][obs_index] = {}
			except KeyError:
				# If no latitude and longitude coordinates are provided,
				# something is seriously fudged. Discard this observation
				# and print a notification

				print('Error: No coordinates provided for observation
						site.\nEnsure that coordinates and site name
						are input correctly.\nObservation discarded')

				return
			
		# If everything has run so far, see what data the input 
		# dictionary has and input it into the orbit dictionaries

		# Let program know whether to list observational data under
		# site name or coordinates
		if site_name != None:
			obs_ident = site_name
		else:
			obs_ident = coordinates
			
		try:
			# Request time of observation
			time = obs_dict['time']
			self.orbit_info[obs_ident][obs_index]['time'] = time
		except KeyError:
			pass

		try:
			# Request range data
			mag_qoppa = obs_dict['range']
			self.orbit_info[obs_ident][obs_index]['range'] = mag_qoppa
		except KeyError:
			pass

		try:
			# Request azimuth
			az = obs_dict['azimuth']
			self.orbit_info[obs_ident][obs_index]['azimuth'] = az
		except KeyError:
			pass

		try:
			# Request elevation
			el = obs_dict['elevation']
			self.orbit_info[obs_ident][obs_index]['elevation'] = el
		except KeyError:
			pass

	def	Orbit_Build(self):
		'''Takes information gathered from observations of object and 
		   determines whether or not an orbit can be built for the object'''
		num_obs = 0
		num_ranges = 0
		num_az = 0
		num_el = 0

		ranges = []
		azimuthes = []
		elevations = []
		# Retrieve the number of each type of observation made	
		for coordinates in self.orbit_info:
 			# Compile number of total observations
			num_obs += len(self.orbit_info[coordinates])		
			for obs_index in self.orbit_info[coordinates]: 	
			
				# Compile number of range estimates 
				try:
					ranges.append(
						self.orbit_info[coordinates][obs_index]['range'])
					num_ranges += 1	
				except KeyError:
					pass
		
				# Compile number of azimuth estimates
				try:
					azimuthes.append(
						self.orbit_info[coordinates][obs_index]['azimuth'])
					num_az += 1
				except KeyError:
					pass

				# Compile number of elevation estimates
				try:
					elevations.append(
						self.orbit_info[coordinates][obs_index]['elevation'])
					num_el += 1
				except KeyError:
					pass

		# Now that the number of each data type has been collected,
		# the best method to solve for the orbit can be determined.
		# If not enough information is present to calculate the 
		# state vectors of the orbit, return and wait for additional
		# observations

		if num_ranges == 3:
			
						
				

			
					
				
				 		 
						 
						
		
		
	
		

		

		
