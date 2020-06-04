import numpy as np

class Stumpff:

	def __init__(self,z):
		self.z = z

	def StumpffS(self):
		'''S(z) Stumpff function, defined by infinite series'''
		if self.z > 0:
			return (np.sqrt(self.z)-np.sin(np.sqrt(self.z)))/(np.sqrt(self.z))^3
		elif self.z < 0:
			return  (np.sinh(np.sqrt(-self.z))-np.sqrt(-self.z))/(np.sqrt(
											-self.z)^3)
		else
			return 1/6

	def StumpffC(self):
		'''C(z) Stumpff function, defined by infinite series'''
		if self.z > 0:
			return (1-np.cos(np.sqrt(self.z)))/self.z
		elif z < 0:
			return (np.cosh(np.sqrt(-self.z))-1)/-self.z
		else
			return 1/2


