import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt

class Grid: 
    def __init__(self, x=[], y=[], z=[]): 
        self.x, self.y, self.z = x, y, z
        self.has_grid_center = False
        
    def get_xyz_array(self): 
        return self.x, self.y, self.z 

    def get_xyz_set(self):
        return np.array([self.x, self.y, self.z]).transpose()

    def get_subgrid_xyzi(self, xmin, xmax, ymin, ymax, zmin, zmax): 
        ix = np.where((self.x > xmin) & (self.x < xmax) & 
                      (self.y > ymin) & (self.y < ymax) & 
                      (self.z > zmin) & (self.z < zmax))
        return self.x[ix], self.y[ix], self.z[ix], ix

    def get_grid_center(self): 
        if not self.has_grid_center: 
            self.initialize_grid_center()
        return self.x0, self.y0, self.z0 
        
    def initialize_grid_center(self): 
        self.x0 = self.x[np.argmin(abs(self.x - abs(1/2*(max(self.x)+min(self.x)))))]
        self.y0 = self.y[np.argmin(abs(self.y - abs(1/2*(max(self.y)+min(self.y)))))] 
        self.z0 = self.z[np.argmin(abs(self.z - abs(1/2*(max(self.z)+min(self.z)))))] 
        self.initialized_grid_center = True

    def scale_xyz(self, r0): 
        self.x = self.x / r0 
        self.y = self.y / r0 
        self.z = self.z / r0 
        self.initialize_grid_center()
        
        
class COMSOLGrid(Grid): 
    def __init__(self, df, sim_unit=1e-6,
                 header_x='% x', header_y='y', header_z='z', 
                 **kwargs): 
        super().__init__()
        self.x = np.array(df[header_x]) * sim_unit
        self.y = np.array(df[header_y]) * sim_unit
        self.z = np.array(df[header_z]) * sim_unit
        self.x_step, self.y_step, self.z_step = self.initialize_resolution()

    def initialize_resolution(self):  
        dx, dy, dz = 0, 0, 0
        for i in range(len(self.x)-1):
            dx0 = self.x[i+1] - self.x[i] 
            dy0 = self.y[i+1] - self.y[i] 
            dz0 = self.z[i+1] - self.z[i]
            dx = dx0 if dx0 > 0 else dx 
            dy = dy0 if dy0 > 0 else dy 
            dz = dz0 if dz0 > 0 else dz 
            if dx > 0 and dy > 0 and dz > 0: 
                break
        return dx, dy, dz

    def get_resolution(self): 
        return self.x_step, self.y_step, self.z_step

    def gen_subgrid(self, xmin, xmax, ymin, ymax, zmin, zmax): 
        x, y, z, _ = self.get_subgrid_xyzi(xmin, xmax, ymin, ymax, zmin, zmax) 
        return Grid(x, y, z)

    def gen_subcube(self, L_cube): 
        x0, y0, z0 = self.get_grid_center() 
        x, y, z, _ = self.get_subgrid_xyzi(x0-L_cube/2, x0+L_cube/2, 
                                           y0-L_cube/2, y0+L_cube/2, 
                                           z0-L_cube/2, z0+L_cube/2)
        return CubeGrid(x, y, z, L_cube)
            

class CubeGrid(Grid): 
    def __init__(self, x, y, z, L_cube): 
        super().__init__(x, y, z) 
        self.L_cube = L_cube
    
class CustomizedGrid(Grid): 
    def __init__(self, xmin, xmax, xstep, ymin, ymax, ystep, zmin, zmax, zstep): 
        super().__init__() 
        self.x = np.arange(xmin, xmax, xstep) 
        self.y = np.arange(ymin, ymax, ystep) 
        self.z = np.arange(zmin, zmax, zstep) 

class CustomCubeGrid(CustomizedGrid):
    def __init__(self, rmax, rmin, rstep): 
        super().__init__(rmin, rmax, rstep, rmin, rmax, rstep, rmin, rmax, rstep)
        self.L = rmax - rmin 