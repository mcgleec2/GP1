import sys
import argparse

sys.argv = sys.argv[:1]

parser = argparse.ArgumentParser()
parser.add_argument( '--mayavi', help='Use mayavi',
		    default=False,action='store_true')
#parser.add_argument( '--material', action = 'store', type = str
#                     , help = 'Name of material.',default='u235')
parser.add_argument( '--shape', action = 'store', type = str
                     ,help = 'Name of shape.',default='wire_2d')
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

class exponential:
    """Fit an exponential to supplied data, and return the
       best fit coefficients and one standard deviation error.
    """
    @staticmethod
    def f( x , p0 , p1 ):
        return p0*np.exp( x*p1 )

    def fit( self , x , y ):
        p0 = [ y[0] , 0.0 ]
        p,cov=curve_fit(self.f,x,y,p0=p0)
        return p, np.sqrt( np.diag( cov ) )

class linear:
    """Fit an linear function to supplied data, and return the
       best fit coefficients and one standard deviation error.
    """
    @staticmethod
    def f( x , p0 , p1 ):
        return p0*(x*p1 + 1.0)

    def fit( self , x , y ):
        p0 = [ y[0] , 0.0 ]
        p,cov=curve_fit(self.f,x,y,p0=p0)
        return p, np.sqrt( np.diag( cov ) )

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
		
def choose_shape( name , length, width ):
    if name == 'wire_2d':
        return wire_2d( length, width )
    else:
        sys.stderr.write( 'Unknown shape: {}\n'.format( name ) )
        sys.exit( 1 )


class ElectronElectron:
    def __init__(self, wire, max_paths, queue_length):
        self.wire = wire
        self.max_paths = max_paths
        self.queue_length = queue_length
        self.electrons = []

    def scatter(self, position, queue, electrons):
        k_e=8.99e9
        e_charge=1.6e-19

        """Move the electron in a random direction in 2D."""
        move_x, move_y, move_z = unit_vector()
        new_position = (position[0] + move_x * random.uniform(0, 1),
                    position[1] + move_y * random.uniform(0, 1))

        for electron in electrons:
          if electron != position:
              distance = np.linalg.norm(np.array(new_position) - np.array(electron))  # distance between electrons
              if distance < 1e-10:  # if eelectrons are too close:
                  # Repulsion force using coulomb
                  force = k_e * (e_charge ** 2) / distance ** 2
                  # modify the direction because of coulombs law
                  move_x += force * (new_position[0] - electron[0])
                  move_y += force * (new_position[1] - electron[1])


        ended = False
        if not self.wire.inside(new_position):
            ended = True
        else:
            queue.appendleft(new_position)  # Keep the new position in the queue for future simulation
        return ended	
	
	
def run(self):

        shape = choose_shape(args.shape, args.length, args.width)
        initial_queue_length = args.queue
        maximum_paths = args.paths

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

        while len(queue) > 0 and path_counter < maximum_paths:

            r = queue.pop()

            while not self.scatter(r, queue, self.electrons):
                shape.sample(r)

            if path_counter%update_interval == 0:
                queue_length.append( len( queue ) )
                paths.append( len(queue_length)*update_interval )
                plt.figure( 1 )
                plt.cla()
                plt.tick_params(axis='both', which='major', labelsize=plt_labsiz)
                plt.plot( paths , np.array( queue_length ) )
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

        if len( queue ) == 0:
            sys.stderr.write( 'Stopped.  Queue empty.\n' )
        else:
            sys.stderr.write( 'Stopped.  Maximum paths reached.\n' )
        sys.stderr.flush()  







wire = wire_2d(args.length, args.width)
simulation = ElectronElectron(wire, args.paths, args.queue, )

# Run the simulation
simulation.run()


plt.ioff()
plt.show()




# Constants
e = 1.602e-19

# Data and data lists
update_interval = 100
path_counter =  0
queue_lenght =  []
paths        =  []
grow         =  []
grow_err     =  []




