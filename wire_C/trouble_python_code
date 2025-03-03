'''
*System description:
Electron movement in 1D inside a 1D wire. In the code you can find some
sections that have been prepared to run the simulation for a 1D movement
of the elctrons inside a 2D wire

Also, even though the electron-electron interaction is being evaluated by 
couloumbs law, the position of the electrons is not updated. This will be 
new positions will be implemeted after a visual output of the simulation is obtained
'''

'''
References:
[1] This code was written by using as a reference the critical_cylinder.py 
code written by Miles Turner, PHY1063-Computational Phydics, DCU, 2024
'''

import sys
import argparse

sys.argv = sys.argv[:1]

parser = argparse.ArgumentParser()
parser.add_argument( '--mayavi', help='Use mayavi',
		    default=False,action='store_true')
#parser.add_argument( '--material', action = 'store', type = str
#                     , help = 'Name of material.',default='u235')
parser.add_argument( '--shape', action = 'store', type = str
                     ,help = 'Name of shape.',default='wire_1D')
parser.add_argument( '--length', action = 'store', type = float
                     , help = 'Characteristic length.',default=5.745868e-2)
parser.add_argument( '--width', action = 'store', type = float
                     , help = 'Characteristic width.',default=5.745868e-2/2)
parser.add_argument( '--paths', action = 'store', type = int
                     , help = 'Trajectories to sample',default=1000000) # whrn the que name is too small the code wont run for the default val but less
parser.add_argument( '--queue', action = 'store', type = int
                     , help = 'Queue length',default=10000)   # this is the only value we are changing
args    = parser.parse_args()


import numpy as np
import random
from   collections        import deque
import matplotlib.pyplot  as     plt
from   matplotlib         import rcParams
from mayavi           import mlab
from   scipy.optimize     import curve_fit
import sys

plt_labsiz = 20

#~~~~~~~ create a random unit vector (direction of initial 
# movemement of the electrons) in 1D -> (u,0,0), ~~~~~~~~~#	 
def unit_vector():
	# the component x is obtained as if a 2D or 3D total vector was obtained 
	u = np.empty()
	cos_theta = 1.0 -2.0*random.random()
	sin_theta = (1.0 - cos*theta**2)**0.5
	phi = 2.0*np.pi*random.random()
	u[0] = cos_theta
	#for 1D:
	u[1], u[2] = 0, 0
	# for 2D, 3D:
	# u[1] = sin_theta+np.cos(phi)
	# u[2] = sin_theta+np.sin(phi)
	return u


# ~~~~~~~ Fit an linear function to supplied data, and return the
# best fit coefficients and one standard deviation error. ~~~~~~~#
# class copied from  the critical_cylinder.py code [1]
class linear:  
    @staticmethod
    def f( x , p0 , p1 ):
        return p0*(x*p1 + 1.0)

    def fit( self , x , y ):
        p0 = [ y[0] , 0.0 ]
        p,cov=curve_fit(self.f,x,y,p0=p0)
        return p, np.sqrt( np.diag( cov ) )



 #~~~~~~~ class to implement the shape interface for a 2D wire ~~~~~~~#
class wire_1D:
	def __init__(self, lenght, width):
		n_bins  = 100                               # number of bin for elctron sampling inside the wire
		self.lenght = lenght                        # wire lenght
		self.width = width                          # wire width
		self.bins = np.zeros(n_bins)                # create the bins
		self.x = np.arange(n_bins) * lenght/n_bins  # set the x-axis position of bins

	# check if the array r is inside the wire, keep running if True, stop if False
	def inside(self, r):
		if  0.0 <= r[0] <= self.lenght and (-self.width/2) <= r[1] <= (self.width/2): 
			return True
		else:
			return False

	# return a random point as the initial position of the electrons inside the wire
	def random_point(self):
		x = random.random()*self.lenght #check lenght or width
		y = 0.0 # for 1D wire
		#y = (random.random() *2 -1) *self.width /2 #for 2D wire
		return [x, y]

	# Save the position r[0] to their respect bin
	def shape(self, r):
		i = int(r[0]/self.lenght*len(self.bins))  # check the bin that r[0] correspond to
		self.bins[i] += 1                         # increment the count of such bin

	# Plot the elctron density disctribution across the wire
	def plot_density(self, figure_index):
		plt.figure(figure_index)
		plt.cla() #clear previous plot
		plt.xlabel('wire position (cm)', fontsize=plt_labsiz)
		plt.ylabel('Density (au)', fontsize=plt_labsiz)
		plt.tick_params(axis='both', which='major', labelsize=plt_labsiz)
		plt.plot( self.x , self.bins )
		plt.ylim( ymin=0 )


#~~~~~~~ Checks that the shape is described in the code ~~~~~~#
# can facilitate the code in the future if ran for a 2D wire
def choose_shape( name , length, width ):
    if name == 'wire_1D':
        return wire_1D( length, width )
    else:
        sys.stderr.write( 'Unknown shape: {}\n'.format( name ) )
        sys.exit( 1 )


#~~~~~~~ Class that contanins all the definitions needed to evaluate the 
# movement of the electrons due to randomness and electrons-electrons
# interaction ~~~~~~~#
class ElectronElectron:
    def __init__(self, wire, max_paths, queue_length):
        self.wire = wire
        self.max_paths = max_paths
        self.queue_length = queue_length
        self.electrons = []

    def scatter(self, position, queue, electrons):
	# Constants:
        k_e=8.99e9        # Coulomb's constant (nm^2/C^2 )
        e_charge=1.6e-19  # Electrons charge

        # move the elctron in the random direction given by the def unit_vector
        move_x, move_y, move_z = unit_vector()
        new_position = (position[0] + move_x * random.uniform(0, 1),
                    position[1] + move_y * random.uniform(0, 1))
	    
	# Evaluate Coulomb's repulsion of the elcectrons and the ones sourranding it
	for electron in electrons:
          if electron != position: # avoid electrons which are in the same location, or self-interaction
              distance = np.linalg.norm(np.array(new_position) - np.array(electron))  # distance between electrons
              	   # if electrons are too close,  calculate repulsion, else: no interaction
	       if distance < 1e-10: 
		  force = k_e * (e_charge ** 2) / distance ** 2
		  # modify the direction because of coulombs law
		  move_x += force * (new_position[0] - electron[0])
		  move_y += force * (new_position[1] - electron[1])
		       
		  #update the position of the electron after Coulomb interaction
	    	  #new_position[0] += force_x 
		  #new_position[1] += forfe_y
		       
	ended = False
	# make sure electrons are still inside of the wire
        if not self.wire.inside(new_position):
            ended = True # dont save the new position
        else:
            queue.appendleft(new_position)  # Store the new position
        return ended	




#~~~~~~~run the simulation~~~~~~~~~#	
def run(self):
#shape is chosen with parameters
        shape = choose_shape(args.shape, args.length, args.width)
        initial_queue_length = args.queue
        maximum_paths = args.paths
#determine the number of points and put them into queue
        queue = deque()
        for i in range(initial_queue_length):
            queue.appendleft( shape.random_point() )

        plt.ion()
        update_interval = 10000  # Gather data after this number of trajectories
        path_counter    = 0      # Number of trajectories examined
        queue_length    = []     # Number of test particles in queue (for plotting)
        paths           = []     # Number of trajectories examined (for plotting)
        grow            = []     # Growth rate (for plotting)
        grow_err        = []
#loop is run while electrons are still in the queue until the max number of paths is reached
        while len(queue) > 0 and path_counter < maximum_paths:

            r = queue.pop()
#scatter function to demonstrate electron's movement through the wire
            while not self.scatter(r, queue, self.electrons):
                shape.sample(r)
#checks the queue length and paths, plots the electrons in queue vs the number of paths
            if path_counter%update_interval == 0:
                queue_length.append( len( queue ) )
                paths.append( len(queue_length)*update_interval )
                plt.figure( 1 )
                plt.cla()
                plt.tick_params(axis='both', which='major', labelsize=plt_labsiz)
                plt.plot( paths , np.array( queue_length ) )
#uses expodential plot to show the growth rate and errors                
		if len(paths) > 2:
                    p , s = exponential().fit( paths, np.array( queue_length ) )
                    grow.append( p[1] )
                    grow_err.append( s[1] )
                    if True:
                        plt.plot( paths , p[0]*np.exp( p[1]*np.array( paths ) ) )
                        plt.plot( paths , p[0]*np.exp( (p[1]+s[1])*np.array( paths ) ) )
                        plt.plot( paths , p[0]*np.exp( (p[1]-s[1])*np.array( paths ) ) )
                    else:
                        plt.plot( paths , p[0]*(p[1]*np.array( paths ) + 1.0) )
                        plt.plot( paths , p[0]*((p[1]+s[1])*np.array( paths ) + 1.0) )
                        plt.plot( paths , p[0]*((p[1]-s[1])*np.array( paths ) + 1.0) )
                    plt.xlabel( r'Paths' , fontsize=plt_labsiz)
                    plt.ylabel( r'Queued Test Particles' , fontsize=plt_labsiz)
                    plt.figure( 2 )
                    plt.cla()
                    plt.errorbar( paths[2:] , grow , yerr=grow_err , fmt='o' )
                    plt.xlabel( r'Paths' , fontsize=plt_labsiz)
                    plt.ylabel( r'Growth Rate' , fontsize=plt_labsiz)
                    plt.tick_params(axis='both', which='major', labelsize=plt_labsiz)
                plt.pause( 0.001 )
            path_counter += 1
        shape.plot_density( 3 )
#indicate when queue is empty or max number of paths reached and the loop stops
        if len( queue ) == 0:
            sys.stderr.write( 'Stopped.  Queue empty.\n' )
        else:
            sys.stderr.write( 'Stopped.  Maximum paths reached.\n' )
        sys.stderr.flush()  





wire = wire_1D(args.length, args.width)
simulation = ElectronElectron(wire, args.paths, args.queue, )

# Run the simulation
simulation.run()


plt.ioff()
plt.show()
