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
