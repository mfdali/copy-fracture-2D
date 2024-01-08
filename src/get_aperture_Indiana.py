import csv
import matplotlib.pyplot as plt
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

def load_data(filename, x_col, y_col, z_col):
    x = []
    y = []
    z = []

    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader)  # Skip the header row if present
        for row in csvreader:
            x.append(float(row[x_col]))
            y.append(float(row[y_col]))
            z.append(float(row[z_col]))
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

def il2b_params():

    a = 687
    b = 271
    input_folder = '/mnt/g/My Drive/Petrobras - Meio Fraturado/2023/fracture Indiana/'
    filename_up = 'hanging_wall_IL2.csv'
    filename_down = 'foot_wall_IL2.csv'
    x_col = 1
    y_col = 0
    z_col = 2
    x_steps = 4662
    y_steps = 1821
    return x_col, y_col, z_col, x_steps, y_steps, input_folder, filename_up, filename_down

def ilta_params():
    input_folder = '/mnt/g/My Drive/Petrobras - Meio Fraturado/2023/fracture Indiana/geometry/ILTA/'
    filename_up = 'ILTA_B_topology.csv'
    filename_down = 'ILTA_A_topology.csv'
    x_col = 1
    y_col = 3
    z_col = 2
    x_steps = 1466
    y_steps = 9596

    return x_col, y_col, z_col, x_steps, y_steps, input_folder, filename_up, filename_down


output_folder = '/mnt/g/My Drive/Petrobras - Meio Fraturado/2023/fracture Indiana/'

sample = 'ILTA'

if sample.lower() == 'ilta':
    x_col, y_col, z_col, x_steps, y_steps, input_folder, filename_up, filename_down = ilta_params()
if sample.lower() == 'il2b':    
    x_col, y_col, z_col, x_steps, y_steps, input_folder, filename_up, filename_down = il2b_params()

upper_wall = input_folder + filename_up
lower_wall = input_folder + filename_down

xup, yup, zup = load_data(upper_wall, x_col, y_col, z_col)
xlow, ylow, zlow = load_data(lower_wall, x_col, y_col, z_col)

xmin = min(np.min(xup),np.min(xlow))
xmax = max(np.max(xup),np.max(xlow))
ymin = min(np.min(yup),np.min(ylow))
ymax = max(np.max(yup),np.max(ylow))

num_points_in_line = len(yup)
min_distances_up = []

#for i in range(num_points_in_line):
i = 20
point1 = (xup[i], yup[i])  # Coordinates of the current point
min_distance = float("inf")  # Initialize with a large value
for j in range(num_points_in_line):
    if i != j:  # Skip the current point
        point2 = (xup[j], yup[j])  # Coordinates of another point
        dist = distance.euclidean(point1, point2)
        min_distance = min(min_distance, dist)
min_distances_up.append(min_distance)

x = np.linspace(xmin, xmax, x_steps)
y = np.linspace(ymin, ymax, y_steps)
X,Y = np.meshgrid(x,y)

Zup = matrixZ(X,Y,xup,yup,zup,output_folder,"hanging_IL2")
Zdown = matrixZ(X,Y,xlow,ylow,zlow,output_folder,"foot_IL2")

aperture = calculate_distance(Zup, Zdown)
np.savetxt(output_folder + 'aperture_Indiana_' + sample + '.txt', np.column_stack((X.ravel(), Y.ravel(), aperture.ravel())), fmt="%.8e", delimiter=",")  