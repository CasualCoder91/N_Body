import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from PIL.Image import core as Image

def plot_cube(cube_definition):
    cube_definition_array = [
        np.array(list(item))
        for item in cube_definition
    ]

    points = []
    points += cube_definition_array
    vectors = [
        cube_definition_array[1] - cube_definition_array[0],
        cube_definition_array[2] - cube_definition_array[0],
        cube_definition_array[3] - cube_definition_array[0]
    ]

    points += [cube_definition_array[0] + vectors[0] + vectors[1]]
    points += [cube_definition_array[0] + vectors[0] + vectors[2]]
    points += [cube_definition_array[0] + vectors[1] + vectors[2]]
    points += [cube_definition_array[0] + vectors[0] + vectors[1] + vectors[2]]

    points = np.array(points)

    edges = [
        [points[0], points[3], points[5], points[1]],
        [points[1], points[5], points[7], points[4]],
        [points[4], points[2], points[6], points[7]],
        [points[2], points[6], points[3], points[0]],
        [points[0], points[2], points[4], points[1]],
        [points[3], points[6], points[7], points[5]]
    ]

    faces = Poly3DCollection(edges, linewidths=1, edgecolors='k')
    faces.set_facecolor((0,0,1,0.1))

    ax.add_collection3d(faces)

    # Plot the points themselves to force the scaling of the axes
    ax.scatter(points[:,0], points[:,1], points[:,2], s=0)

    #ax.set_aspect('equal')

def plt_star_and_cube_file():
    cubes = Import=np.loadtxt('tree.dat',delimiter=',')
    stars = Import=np.loadtxt('stars.dat',delimiter=',')
    #print(cubes)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for line in cubes:
        cube_definition = [
            (line[0],line[4],line[2]), (line[0],line[1],line[2]), (line[3],line[4],line[2]), (line[0],line[4],line[5])
        ]
        plot_cube(cube_definition)
    ax.scatter(stars[:,0], stars[:,1], stars[:,2], s=50,c='red')
    plt.show()

def plot_star_series(n_plots,stepsize):
    for i in np.arange(0,n_plots,stepsize):
        name = 'stars'+str(i)
        stars = Import=np.loadtxt(name+'.dat',delimiter=',')
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(stars[:,0], stars[:,1], stars[:,2], s=50,c='red')
        plt.savefig(name+'.jpg')
        plt.close(fig)

def plot_star_series_with_label(n_plots,stepsize):
    for i in np.arange(0,n_plots,stepsize):
        name = 'stars'+str(i)
        stars = Import=np.loadtxt(name+'.dat',delimiter=',')
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(stars[:,0], stars[:,1], stars[:,2], s=50,c='red')
        ax.set_xlim([4.0,4.0])
        ax.set_ylim([4.0,4.0])
        ax.set_zlim([4.0,4.0])
        for x,y,z,label in zip(stars[:,0], stars[:,1], stars[:,2], stars[:,3]):
             ax.text(x, y, z, label)
        plt.savefig(name+'.jpg')
        plt.close(fig)

plot_star_series(100000,100)

#stars = Import=np.loadtxt('stars0.dat',delimiter=',')
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(stars[:,0], stars[:,1], stars[:,2], s=50,c='red')
#plt.show()
