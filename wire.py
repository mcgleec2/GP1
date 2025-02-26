

import numpy as np
import random



# unit vector: 1D -> (u,0,0)
def unit_vector():
	# return an array with a random unit vector in 1D
	u = np.empty()
	cos_theta = 1.0 -2.0*random.random()
	sin_theta = (1.0 - cos*theta**2)**0.5
	phi = 2.0*np.pi*random.random()
	u[0] = cos_theta
	#for 1D:
	u[1] = 0
	# u[2] = 0
	# for 2D, 3D:
	# u[1] = sin_theta+np.cos(phi)
	# u[2] = sin_theta+np.sin(phi)
	return u

class wire_2d:
	# class to implement the shape interfade for a 2D wire
	def __init__(self, lenght, width):
		n_bins  = 100
		self.lenght = lenght
		self.width = width
		self.bins = np.zeros(n_bins)
		self.x = np.arange(n_bins) * lenght/n_bins   #check lenght or width

	def inside(self, r):
		if  0.0 <= r[0] <= self.lenght and (-self.width/2) <= r[1] <= (self.width/2): 
			return True
		else:
			return False

	def random_point(self):
		# return a random point inside the 2D wire:
		x = random.random()*self.lenght #check lenght or width
		#y = (random.random() *2 -1) *self.width /2 #for 2D wire
		y = 0.0 # for 1D wire
		return [x, y]

	def sample(self, r):
		i = int(r[0]/self.lenght*len(self.bins))
		self.bins[i] += 1

	def plot_density(self, figure_index):
		plt.figure(figure_index)
		plt.cla()
		plt.xlabel('wire position (cm)', fontsize=plt_labsiz)
		plt.ylabel('Density (au)', fontsize=plt_labsiz)
		plt.tick_params(axis='both', which='major', labelsize=plt_labsiz)
		plt.plot( self.x , self.bins )
		plt.ylim( ymin=0 )
		
	
	
	
  











# Constants
e = 1.602e-19

# Data and data lists
update_interval = 100
path_counter =  0
queue_lenght =  []
paths        =  []
grow         =  []
grow_err     =  []




