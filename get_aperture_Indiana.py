import csv
import matplotlib.pyplot as plt
import csv
import numpy as np
from scipy.interpolate import griddata
from scipy.spatial import distance
'''
Andesite
files = 'surfaces_hang.csv' and 'surfaces_foot.csv'
x.append(float(row[3]))
y.append(float(row[1]))
z.append(float(row[2]))


'''

def load_data(filename):
    x = []
    y = []
    z = []

    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader)  # Skip the header row if present
        for row in csvreader:
            x.append(float(row[1]))
            y.append(float(row[0]))
            z.append(float(row[2]))
            #print("%e,%e,%e\n"%(float(row[1]),float(row[0]),float(row[2])))

    return np.array(x), np.array(y), np.array(z)

def calculate_z(x, y, z, x_query, y_query):
    points = np.column_stack((x, y))
    values = z
    z_query = griddata(points, values, (x_query, y_query), method='nearest')
    return z_query

def matrixZ(X,Y,x,y,z,output_folder,case_filename):
 
    Z = calculate_z(x, y, z, X, Y)
    #np.savetxt(output_folder + "Z_" + case_filename + ".csv", np.column_stack((X.ravel(), Y.ravel(), Z.ravel())), fmt="%.2g", delimiter=",")

    return Z

def calculate_distance(Zup, Zdown):
    distance = np.sqrt((Zup - Zdown) ** 2)
    return distance

a = 687
b = 271
output_folder = '/home/monique/fem-fenics/two-phase/fracture_results/'

upper_wall = output_folder + 'hanging_wall_IL2.csv'
lower_wall = output_folder + 'foot_wall_IL2.csv'

x = np.linspace(-0.0473847, 0.0509658, 4662)
y = np.linspace(-0.0175, 0.018, 1821)
X,Y = np.meshgrid(x,y)

xup, yup, zup = load_data(upper_wall)
xlow, ylow, zlow = load_data(lower_wall)

Zup = matrixZ(X,Y,xup,yup,zup,output_folder,"hanging_IL2")
Zdown = matrixZ(X,Y,xlow,ylow,zlow,output_folder,"foot_IL2")

aperture = calculate_distance(Zup, Zdown)
np.savetxt(output_folder + 'aperture_Indiana_IL2.csv', np.column_stack((X.ravel(), Y.ravel(), aperture.ravel())), fmt="%.8e", delimiter=",")  