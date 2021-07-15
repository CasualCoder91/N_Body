from point import Point
import numpy as np
from config import eps_magnitude

def generate_velocity_and_index(points_t0, points_t1):
	""" 
    Extracting the stars(=points) at t1 from the image online gives information about position and magnitude but not velocity or id of the star.
	This function finds the corresponding stars from the previous timestep and sets the id of the future timestep 
	as well as the velocity of the stars from the previous timestep accordingly.

    :returns: updated point arrays
    """
	max_dist_pos = 0
	for i in range(len(points_t0)-1): # loop through all points at timestep 0
		minDist = -1
		future_point_index = -1

		for j in range(len(points_t1)-1): #compare to all points at timestep 1
			current_dist = points_t0[i].get_distance(points_t1[j])
			if (currentDist < minDist or minDist == -1) and abs(1 - points_t1[j].magnitude / points_t0[i].magnitude) < eps_magnitude:
				minDist = currentDist
				future_point_index = j
			if maxDistPos < currentDist: 
				maxDistPos = currentDist
		points_t0[i].velocity[0] = points_t1[future_point_index].position[0] - points_t0[i].position[0]; #tecnically division by dt needed but dt is equal for all points
		points_t0[i].velocity[1] = points_t1[future_point_index].position[1] - points_t0[i].position[1];
		points_t1[future_point_index].id = points_t0[i].id
	return points_t0, points_t1

def estimate_MKs():
	x = [19.8,18.5,17.7,15,11,10,7.3,6.0,5.4,5.1,4.7,4.3,3.92,3.38,2.75,2.68,2.18,2.05,1.98,1.86,1.93,1.88,1.83,1.77,1.81,1.75,1.61,1.50,1.46,1.44,1.38,1.33,1.25,1.21,1.18,1.13,1.08,1.06,1.03,1.00,0.99,0.985,0.98,0.97,0.95,0.94,0.90,0.88,0.86,0.82,0.78,0.73,0.70,0.69,0.64,0.62,0.59,0.57,0.54,0.50,0.47,0.44,0.40,0.37,0.27,0.23,0.184,0.162,0.123,0.102,0.093,0.090,0.088,0.085,0.080,0.079,0.078,0.077,0.076,0.075,0.01]
	y = [-3.20,-3.073,-2.942,-2.587,-2.126,-1.848,-1.198,-0.956,-0.708,-0.553,-0.433,-0.192,-0.075,0.254,0.621,0.648,0.949,1.07,1.172,1.447,1.587,1.607,1.655,1.702,1.694,1.756,1.836,1.941,2.045,2.116,2.188,2.291,2.50,2.579,2.76,2.915,2.947,3.043,3.120,3.236,3.282,3.319,3.345,3.409,3.495,3.532,3.693,3.827,3.889,3.938,4.100,4.247,4.397,4.56,4.81,4.95,5.01,5.15,5.36,5.64,5.75,5.98,6.18,6.55,7.10,7.36,7.93,8.20,8.80,9.22,9.50,9.70,9.81,9.92,10.30,10.40,10.50,10.55,10.77,11.00,20.00 ]

	coef = np.polyfit(x[:6],y[:6],1)
	poly1d_fn = np.poly1d(coef) 

	plt.plot(x[:6],y[:6], 'yo', x[:6], poly1d_fn(x[:6]), '--k')

	estimate = poly1d_fn(50)
	print(estimate)

	plt.show()
	return