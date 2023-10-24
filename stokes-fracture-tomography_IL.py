import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import griddata
from scipy.spatial import distance
from dolfin import *
from mshr import *
import random #pseudorandom
import numpy as np
#import matplotlib.pyplot as pl
from math import pi, sqrt, floor, ceil,pow
#import time
from datetime import datetime
#from secrets import randbelow #not cryptographically secure
import time
import matplotlib.pyplot as plt
from ufl.operators import min_value
from scipy.interpolate import LinearNDInterpolator

class Params:
  def __init__(self):
    self.muw = 0.001
    self.Swc = 0
    self.nw = 2
    self.muo = 0.001
    self.rhoo = 1000
    self.rhow = 1000
    self.kromax = 1
    self.krwmax = 1
    
    self.Sor = 0
    self.no = 2
    self.Cw = 0
    self.Co = 0
    self.kpm = 1E-11
    self.pin = 5*6894.76
    self.pout = 6894.76

class MeshHeight(UserExpression):

    def __init__(self, **kwargs):
        #random.seed(2 + MPI.rank(MPI.comm_world))
        random.seed(datetime.now())
        super().__init__(**kwargs)
        
    def set_h_values(self,X,Y,Z):
        self.X, self.Y, self.Z = X, Y, Z
        self.interpolator = LinearNDInterpolator((self.X, self.Y), self.Z)

    def eval(self, values, x):
        '''tol = 1E-14
        mask = (
            (x[0] >= self.X - tol - 0.00005)
            & (x[0] <= self.X + tol + 0.00005)
            & (x[1] >= self.Y - tol - 0.00005)
            & (x[1] <= self.Y + tol + 0.00005)
        )
        if mask.any():
            values[0] = griddata((self.X, self.Y), self.Z, (x[0], x[1]), method="linear")'''
        f_h = self.interpolator([x[0], x[1]])
        if f_h > 1E-10:
            values[0] = self.interpolator([x[0], x[1]])
        else:
            values[0] = 0

    def value_shape(self):
        return ()

#class LoadAperture():

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

def _calculate_z(x, y, z, x_query, y_query):
    points = np.column_stack((x, y))
    values = z
    z_query = griddata(points, values, (x_query, y_query), method='nearest')
    return z_query

def matrixZ(X,Y,x,y,z):
    Z = _calculate_z(x, y, z, X, Y)
    return Z

def calculate_distance(Zup, Zdown):
    distance = np.sqrt((Zup - Zdown) ** 2)
    return distance

#class NumModel():

def Dolfin_Mesh(res,a,b):

    # Domain size
    larea = a
    harea = b

    # Domain construction
    rectangle = Rectangle(Point(0., 0.), Point(larea, harea))
    domain = rectangle

    # Mesh generation
    mesh = generate_mesh(domain,res)

    # Define boundary condition
    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    left = AutoSubDomain(lambda x: near(x[0], 0.0))
    right = AutoSubDomain(lambda x: near(x[0], larea))
    bottom = AutoSubDomain(lambda x: near(x[1], 0.0))
    top = AutoSubDomain(lambda x: near(x[1], harea))

    # Define boundary markers
    left.mark(boundaries, 1)
    top.mark(boundaries, 2)
    right.mark(boundaries, 3)
    bottom.mark(boundaries, 4)
    '''
    # 2D plot
    plot(mesh)
    #plt.axis('off')
    plt.show()
    plt.savefig("mesh-savefig.png")'''
    return mesh, boundaries

def Max(a, b): return (a+b+abs(a-b))/Constant(2)
def Min(a, b): return (a+b-abs(a-b))/Constant(2)

'''def func(x,y,h):
    return h((x,y))'''

def thickness(xah, yah, zah,x_,y_,case_filename):

    from matplotlib import cm
    from matplotlib.ticker import LinearLocator

    eq = case_filename

    #grid_pos = grid_position()
    
    #Aperture
    h = MeshHeight()
    h.set_h_values(xah,yah,zah)
    
    tol = 0.001 # avoid hitting points outside the domain
    x = np.linspace(0, x_, 1001)
    y = np.linspace(0, y_, 1001)
    X,Y = np.meshgrid(x,y)
    #h_min = hmin
    #h_max = hmax
    
    Z = matrixZ(X,Y, xah, yah,zah)
    #pointsy = [(0, y_) for y_ in y] # 2D points
    #pointsx = [(x_,0) for x_ in x] # 2D points
    #w_line = np.array([inflow(point) for point in pointsx])
    #p_line = np.array([inflow(point) for point in pointsy])
    '''h_min = np.min(Z)
    h_max = np.max(Z)
    h_av = np.average(Z)'''

    #fig, axs = plt.subplots(2, 2, figsize=(9, 3))
    fig = plt.figure(figsize=(9, 4))

    # Plot the surface.

    ax = fig.add_subplot(1, 2, 1, projection='3d')
    surf = ax.plot_surface(X, Y, Z, cmap=cm.viridis, rstride=4, cstride=4, 
                        linewidth=0.2)#, antialiased=False)
    #ax.set_aspect('equal')
    ax.view_init(30, 70)
    # Plot projections of the contours for each dimension. By choosing offsets
    # that match the appropriate axes limits, the projected contours will sit on
    # the 'walls' of the graph
    #ax.set(xlim=(0.0, 0.2), ylim=(0.0, 0.05), zlim=(0.0,0.02), xlabel='X', ylabel='Y', zlabel='Z')

    #ax.contourf(X, Y, Z, zdir='z', offset=0, cmap='jet')
    #ax.contour(X, Y, Z, zdir='x', offset=X.min(), cmap='coolwarm')
    #ax.contourf(X, Y, Z, zdir='y', offset=0.05, cmap='coolwarm')

    fig.colorbar(surf, shrink=0.5, aspect=10, location='bottom',orientation='horizontal')

    #ax = fig.add_subplot(2, 2, 2)
    #ax.plot(x, w_line, 'k', linewidth=2) # magnify w
    #ax.plot(y, p_line, 'b--', linewidth=2)
    
    
    #plt.savefig(output_folder + case_filename+'-expression.png')
    
    ax = fig.add_subplot(1, 2, 2)
    CS = ax.contour(X, Y, Z)
    ax.clabel(CS, inline=True, fontsize=10)
    ax.grid(True)
    #ax.set_title('Simplest default with labels')
    #plt.savefig(output_folder + case_filename+'_-contour.png')
    
    #fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    
    '''
    # Customize the z axis.
    ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    # A StrMethodFormatter is used automatically
    ax.zaxis.set_major_formatter('{x:.02f}')
    '''
    # Add a color bar which maps values to colors.
    #fig.colorbar(surf, shrink=0.5, aspect=5)
    #plt.savefig(output_folder + case_filename+'-surf.png')
    fig.savefig(output_folder + case_filename + '-all.png',dpi=120)

    plt.close()
    
    return h,eq#,h_av,h_min,h_max


def Num_Solution(x_,y_,xah, yah, zah,case_filename,output_folder):
    
    ### Parameters ###
    res = 128
    start_time = time.time()

    # Load mesh
    mesh, boundaries = Dolfin_Mesh(res,x_,y_)

    #Material properties
    Prop = Params()
    mu = Constant(Prop.muw) #Water viscosity [Pa.s]
    mu_ = 7.5*mu #R.C. Givler, S.A. Altobelli, J. Fluid Mech. 258 (1994) 355.
    pin = Constant(Prop.pin) #Imposed pressure at the entrance [Pa]
    pout = Constant(Prop.pout) #Imposed pressure on output [Pa]
    thick,eq = thickness(xah, yah, zah, x_, y_, case_filename) #Varing the height

    ####### Numerical Solution #######

    #Function space over the mesh
    V = VectorElement('CG',mesh.ufl_cell(),2)
    Q = FiniteElement('CG',mesh.ufl_cell(),1)
    Element = V*Q
    W = FunctionSpace(mesh,Element)

    Y = FunctionSpace(mesh, "Lagrange", 1)
    h_test = interpolate(thick,Y)
    h_min = h_test.vector().min()
    h_max = h_test.vector().max()

    #Define variational problem
    (u,p) = TrialFunctions(W)
    (v,q) = TestFunctions(W)
                    
    ##############################################################
    # Parameters

    # Define expressions used in variational forms
    dp = pin-pout
    #g = Expression("b-(a/l)*y[0]" , degree=1 ,a = Constant(dp), b = Constant(pin), l = 0.027) #g = div(u)
    #u_in = Constant((0.0)) #Initial velocity in x [m/s]
    noslip = Constant((0.0,0.0)) #No-slip condition for velocity, u=0 at y=h
    f = Constant((0.0,0.0)) #External force
    n = FacetNormal(mesh) #Normal vector to mesh

    ##############################################################

    #Define Dirichlet boundary conditions
    bc1 = DirichletBC(W.sub(0).sub(1),Constant(0.0),boundaries,1)
    bc2 = DirichletBC(W.sub(0),noslip,boundaries,2)
    bc3 = DirichletBC(W.sub(0).sub(1),Constant(0.0),boundaries,3)
    bc4 = DirichletBC(W.sub(0),noslip,boundaries,4)

    bcs = [bc1,bc2,bc3,bc4]

    '''mesh coordinates'''
    D = mesh.topology().dim()
        
    # get vertex coordinates
    coords = mesh.coordinates()
    
    # get elements at BC1
    mks = bc3.markers()
    
    print(len(mks))

    # Find indices where the x-coordinate equals 'a'
    indices_inlet = [i for i, (x, y) in enumerate(coords) if near(x, 0)]
    indices_outlet = [i for i, (x, y) in enumerate(coords) if near(x, x_)]

    # Extract the corresponding coordinates
    coordinates_inlet = coords[indices_inlet]
    coordinates_outlet = coords[indices_outlet]
    
    h_inlet = [h_test((coordinates_inlet[i])) for i in range(len(coordinates_inlet))]
    h_outlet = [h_test((coordinates_outlet[i])) for i in range(len(coordinates_outlet))]

    # Convert the array to a DataFrame
    df1 = pd.DataFrame(h_inlet)
    df2 = pd.DataFrame(h_outlet)

    # Save the DataFrame to a CSV file
    df1.to_csv('h_inlet_andesite.csv', index=False, header=False)
    df2.to_csv('h_outlet_andesite.csv', index=False, header=False)

    #Define measures associated with the boundaries and holes
    ds = Measure('ds',domain=mesh, subdomain_data = boundaries)
    dx = Measure('dx',domain=mesh)

    ##############################################################
    # Compensate the drag at glass walls (3D)
    Drag_force = -(12/(thick**2))*mu*inner(u,v)*dx

    #Define variational form for Stokes
    #a = (mu_*inner(grad(u),grad(v))*dx(1) + (mu/k)*inner(u,v)*dx(0) - div(v)*p*dx(1) - div(v)*p*dx(0)-div(u)*q*dx(1) - div(u)*q*dx(0))
    #L = (inner(f,v)*dx(1) + inner(f,v)*dx(0) - pin*dot(v,n)*ds(1) - pout*dot(v,n)*ds(3))

    F = (mu*inner(grad(u),grad(v))*dx - div(v)*p*dx - div(u)*q*dx) + pin*dot(v,n)*ds(1) + pout*dot(v,n)*ds(3) - Drag_force - inner(f,v)*dx
    (a, L) = system(F)

    #Compute solution
    U = Function(W, name = "field")
    #solve(a,U.vector(),L)
    '''
    problem = LinearVariationalProblem(a, L, U, bcs)
    solver = LinearVariationalSolver(problem)
    solver.solve()
    '''
    A  = assemble(a)
    b = assemble(L)
    [bc.apply(A,b) for bc in bcs]
    #[bc.apply(U.vector()) for bc in bcs]
    solve(A,U.vector(),b,"mumps")

    ##############################################################

    #Get sub functions
    (u, p) = U.split()
    #ux, uy = u.split(deepcopy=True)

    # inlet flow
    form = -dot(n, u) * ds(1)
    inflow = assemble(form)

    # outlet flow
    form = dot(n, u) * ds(3)
    outflow = assemble(form)
    deviation = ((inflow-outflow)/inflow)*100
    print("Flow: %.10g\t %.10g\t deviation: %e\n" % (inflow,outflow,deviation))

    exec_time = time.time() - start_time

    print(exec_time)
    #Save data to file
    #data_file.write("%s,%s,%s,%g,%g,%g,%g,%2.10g,%2.10g,%g\n" % (mode,mesh_type,sm,triangles,exec_time,k,ps,inflow,outflow,k_medium))
    
    
    # Save solution in VTK format
    file = File(output_folder + case_filename + "-u.pvd")
    file << u
    file = File(output_folder + case_filename + "-p.pvd")
    file << p
    file = File(output_folder + case_filename + "-ah.pvd")
    file << h_test

    # Plot field and save figure
    plt.figure()
    plot(u, mode = "glyphs")
    '''
    plot(u,
    wireframe = True,              # use wireframe rendering
    scalarbar = False,             # hide the color mapping bar
    scale = 2.0,                    # scale the warping/glyphs
    title = "Fancy plot"           # Set your own title
    )'''
    
    plt.xlabel('u')
    plt.savefig(output_folder + case_filename + "-u.png")
    plt.close()

    print("\nmin: %.5e"%h_min)
    print("\nmax: %.5e"%h_max)

    a_H2 = assemble(thick*dx) 
    Area = float(assemble(1 * dx))
    int_H2 = a_H2/Area
    print("a_h_average2 = %e" % int_H2)

    try:
        a_H = assemble((1/(thick**3))*dx)
        Area = float(assemble(1 * dx))
        int_H = pow((Area/a_H),(1/3))
        print("int dA/hy^3 = %e" % a_H)
        print("a_h_average = %e" % int_H)
    except:
        a_H = 0
        int_H = 0

    #mesh size
    triangles = mesh.num_cells()
    print("%g\n" % triangles)
    info("Num DOFs {}".format(W.dim()))

    k_eq = (outflow * mu * x_)/(y_*dp)
    k_D = k_eq*1013250000000

    print("k_eq: %.5g\t k_darcy: %g\n" % (k_eq,k_D))

    
    print_to_file.write("%s,%s,%g,%e,%e,%.10g,%.10g,%e,%g,%g,%g,%g,%g\n" % (case_filename,eq,int_H2,h_min,h_max,inflow,outflow,deviation,k_eq,exec_time,Area,a_H,int_H))

'''Andesite'''
def andesite_properties(output_folder):
    a = 687
    b = 271
    x_ = a*0.0001 # Lenght
    y_ = b*0.0001 # Width

    aperture = output_folder + 'aperture.csv'
    project_name = 'stokes-ss-hvar-proj394_'
    case = ['1e-3','9.5e-4','9e-4','8.9e-4']
    h_ = [0,1e-4,1.5e-4,1.55e-4]

    x, y, z = load_data(aperture) 
    xah = x*0.0001
    yah = y*0.0001
    zah = ((z-np.min(z))*0.0001)-(4.94506e-05)#-(1.55e-4)
    return x_,y_,xah, yah, zah, aperture, project_name, case, h_

'''Indiana'''
def indiana_properties(output_folder):
    a = 0.0983505
    b = 0.0355
    x_ = a # Lenght
    y_ = b # Width

    aperture = output_folder + 'aperture_Indiana_IL2.csv'
    project_name = 'stokes-ss-hvar-IL2_'
    case = ['0','1','2']
    h_ = [0,5e-5,7.5e-5]

    x, y, z = load_data(aperture) 
    xah = x+0.0473847
    yah = y+0.0175
    zah = z-np.min(z)-(6.59062e-05)#-(7.5e-5)
    
    return x_,y_,xah, yah, zah, aperture, project_name, case, h_

output_folder = '/home/monique/fem-fenics/two-phase/fracture_results/'

#x_,y_,xah, yah, zah, aperture, project_name, case, h_ = indiana_properties(output_folder)
x_,y_,xah, yah, zah, aperture, project_name, case, h_ = andesite_properties(output_folder)

print_to_file = open(output_folder + project_name + ".txt", "a")
#print_to_file.write("Case,Equation,Av_Height,np.min,np.max,Inflow,Outflow,%Deviation,K_eq,Exec_time,Area,Height_3,Av_lub\n")

for step in range(len(h_)):
    case_filename = project_name + case[step]
    zh_ = zah-h_[step]
    Num_Solution(x_,y_,xah, yah,zh_,case_filename,output_folder)

print_to_file.close()
