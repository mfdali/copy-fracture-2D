# copy-fracture-2D
Get aperture of real fractured core sample

Calculate the gap between two surfaces from .stl fractured core samples

![Core sample](https://github.com/mfdali/copy-fracture-2D/blob/main/ILB-4-15_git.png?raw=true)

1. ## Extract csv file
   From a stl file get the surface points of each side of the core sample. Extract the csv file.

   ### Fracture Surface
   ![Surface A](https://github.com/mfdali/copy-fracture-2D/blob/main/ILB-4-15_A.png?raw=true)

3. ## Calculate distance between surfaces
   
   File: [fracture_aperture_from_core_sample_url.ipynb](https://github.com/mfdali/copy-fracture-2D/blob/main/fracture_aperture_from_core_sample_url.ipynb)
   
   Load the two csv files with the surface points, calculate and save the fracture aperture.
   
   ### Fracture Aperture
   ![Heatmap](https://github.com/mfdali/copy-fracture-2D/blob/main/heatmap.png?raw=true)

   ### Aperture distribution
   ![Histogram](https://github.com/mfdali/copy-fracture-2D/blob/main/histogram.png?raw=true)

5. ## Numerical Simulation
   
   File: [stokes-fracture_core-sample.py](https://github.com/mfdali/copy-fracture-2D/blob/main/stokes-fracture_core-sample.py)
   
   Simulate single-phase flow through fracture using the distance between plates extracted in fracture_aperture_from_core_sample_url.ipynb

   ### Velocity field
   ![Velocity](https://github.com/mfdali/copy-fracture-2D/blob/main/stokes-ss-hvar-ILB_4_15_2e5-0-u.png?raw=true)
