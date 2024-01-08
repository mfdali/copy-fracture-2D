import csv
import numpy as np
import gmsh
import math

def load_data(filename):
    x = []
    y = []
    z = []

    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        #next(csvreader)  # Skip the header row if present
        for row in csvreader:
            x.append(float(row[0]))
            y.append(float(row[1]))
            z.append(float(row[2]))

    return np.array(x), np.array(y), np.array(z)

def indiana_properties(output_folder):
    a = 0.0983505
    b = 0.0355
    x_ = a # Lenght
    y_ = b # Width

    aperture = output_folder + 'aperture_Indiana_IL2.csv'
    project_name = 'stokes-ss-hvar-IL2_gmsh'
    case = ['0']
    h_ = [0]

    x, y, z = load_data(aperture) 
    xah = x+0.0473847
    yah = y+0.0175
    zh_ = z-np.min(z)-(6.59062e-05)#-(7.5e-5)
    
    return x_,y_,xah, yah, zh_, aperture, project_name, case, h_

file_path = "/mnt/g/My Drive/Fenicsx/gmsh/stl/"
x_,y_,xah, yah, zh_, aperture, project_name, case, h_ = indiana_properties(file_path)