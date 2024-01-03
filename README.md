# copy-fracture-2D
Get aperture of real fractured core sample

Calculate the gap between two surfaces from .stl fractured core samples

![Core sample](https://github.com/mfdali/copy-fracture-2D/blob/main/ILB-4-15_git.png?raw=true)

1. ## Extract csv file
   From a stl file get the surface points of each side of the core sample. Extract the csv file.

   ### Fracture Surface
   ![Surface A](https://github.com/mfdali/copy-fracture-2D/blob/main/ILB-4-15_A.png?raw=true)

2. ## Calculate distance between surfaces
   
   File: [fracture_aperture_from_core_sample_url.ipynb](https://github.com/mfdali/copy-fracture-2D/blob/main/fracture_aperture_from_core_sample_url.ipynb)
   
   Load the two csv files with the surface points, calculate and save the fracture aperture.
   
   ### Fracture Aperture
   ![Heatmap](https://github.com/mfdali/copy-fracture-2D/blob/main/heatmap.png?raw=true)

   ### Aperture distribution
   ![Histogram](https://github.com/mfdali/copy-fracture-2D/blob/main/histogram.png?raw=true)

3. ## 2.5D Numerical Simulation
   
   File: [stokes-fracture_core-sample.py](https://github.com/mfdali/copy-fracture-2D/blob/main/stokes-fracture_core-sample.py)
   
   Create 2D mesh and simulate single-phase flow through fracture using the distance between plates extracted in fracture_aperture_from_core_sample_url.ipynb

   ### Velocity field
   ![Velocity](https://github.com/mfdali/copy-fracture-2D/blob/main/stokes-ss-hvar-ILB_4_15_2e5-0-u.png?raw=true)

4. ## Mesh Test
   
   ![MeshTest](https://github.com/lmmp-puc-rio/copy-fracture-2D/blob/main/mesh_test_ILB_4_15.png?raw=true)
   
5. ## Manipulate fracture aperture

   File: [fracture_simulation_analysis.ipynb](ttps://github.com/lmmp-puc-rio/copy-fracture-2D/blob/main/fracture_simulation_analysis.ipynb)
   
   Manipulate the gap between the two surfaces of the fracture. Create another aperture based on the original fracture aperture in step 2.
   
   $new \ aperture = original \ a^2_h - \Delta height$

   Example:
   
   $\Delta h = 0.02 mm$

   | Case     | \Delta height (mm) | \bar{a}^2_h (mm) | max({a}^2_h) (mm) | k_simulation |
   | ----------- | ----------- | ----------- | ----------- | ----------- |
   | 0    | 0.02 | 0.180 | 2.17 |2322.58 |
   | 1    | 0.04 |  0.160 | 2.15 | 1787.41 |
   | 2 | 0.06 | 0.140 | 2.13 | 1323.27 |
   | 3 | 0.08 | 0.120 | 2.11 | 930.06 |
   | 4 | 0.10 | 0.100 | 2.09 | 607.93 |

   Fracture permeability equation
   $k_frac = \frac{\bar{a}^2_h}{12}$
   
   ### Permeability Graph: Simulation x Cubic law
   
   ![HvarGraph](ttps://github.com/lmmp-puc-rio/copy-fracture-2D/blob/main/mesh_test_ILB_4_15.png?raw=true)
