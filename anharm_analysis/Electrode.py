import pandas as pd
import numpy as np 
from .Grid import COMSOLGrid

class SimulatedElectrode: 
    def __init__(self, name, file, sim_unit=1e-6): 
        self.V, self.Ex, self.Ey, self.Ez, self.sim_grid = [], [], [], [], [] 

    def load_from_file(self): 
        raise NotImplementedError 

    def get_V_in_cube(self, V0=1, L_cube=None): 
        if L_cube is None: 
            return V0 * self.V 
        x0, y0, z0 = self.sim_grid.get_grid_center() 
        V_idx = self.sim_grid.get_subgrid_xyzi(x0-L_cube/2, x0+L_cube/2, 
                                               y0-L_cube/2, y0+L_cube/2, 
                                               z0-L_cube/2, z0+L_cube/2)[-1] 
        return V0 * self.V[V_idx]

    def get_V_in_ROI(self,V0=1, ROI_grid=None):
        if ROI_grid is None: 
            return V0 * self.V 
        x0, y0, z0 = ROI_grid.get_grid_center() 
        L_cube = ROI_grid.L_cube
        V_idx = self.sim_grid.get_subgrid_xyzi(x0-L_cube/2, x0+L_cube/2, 
                                               y0-L_cube/2, y0+L_cube/2, 
                                               z0-L_cube/2, z0+L_cube/2)[-1] 
        return V0 * self.V[V_idx]   

    def get_subcube_Vxyz(self, V0=1, L_cube=None): 
        if L_cube is None: 
            return V0 * self.V 
        x0, y0, z0 = self.sim_grid.get_grid_center() 
        X, Y, Z, V_idx = self.sim_grid.get_subgrid_xyzi(x0-L_cube/2, x0+L_cube/2, 
                                               y0-L_cube/2, y0+L_cube/2, 
                                               z0-L_cube/2, z0+L_cube/2)
        return V0 * self.V[V_idx], X, Y, Z

    def set_sim_grid(self, grid): 
        self.sim_grid = grid


class COMSOLElectrode(SimulatedElectrode): 
    def __init__(self, name, file, sim_unit=1e-6, 
                 header_x='% x', header_y='y', header_z='z', electrode_list=None, quantities = ['V', 'Ex', 'Ey', 'Ez'], **kwargs): 
        self.electrode_list = electrode_list
        self.quantities = quantities
        self.V, self.Ex, self.Ey, self.Ez, self.sim_grid = self.load_from_file(name, file, 
                                                                               unit=sim_unit, 
                                                                               header_x=header_x, 
                                                                               header_y=header_y, 
                                                                               header_z=header_z, **kwargs) 

    def load_from_file(self, electrode, file, excitation_prefix='V', sim_prefix='esbe.', skiprows=8, 
                       unit=1e-6, header_x='% x', header_y='y', header_z='z', **kwargs): 
        df = pd.read_csv(file, skiprows=skiprows) 
        if electrode != 'default':
            all_headers = [i for i in df.keys() if f'{excitation_prefix}{electrode}=1' in i] 
            V_headers = [i for i in all_headers if f'{sim_prefix}V' in i] 
            Ex_headers = [i for i in all_headers if f'{sim_prefix}Ex' in i]
            Ey_headers = [i for i in all_headers if f'{sim_prefix}Ey' in i] 
            Ez_headers = [i for i in all_headers if f'{sim_prefix}Ez' in i] 
            assert len(V_headers)==1, "Check output results or output headers, found 0 or 1+ potential data"
        else: 
            V_headers = [i for i in df.keys() if f'{sim_prefix}V' in i] 
            Ex_headers = [i for i in df.keys() if f'{sim_prefix}Ex' in i]
            Ey_headers = [i for i in df.keys() if f'{sim_prefix}Ey' in i] 
            Ez_headers = [i for i in df.keys() if f'{sim_prefix}Ez' in i] 
        V = np.array(df[V_headers[0]]) 
        Ex = np.array([]) if len(Ex_headers) != 1 else np.array(df[Ex_headers[0]]) 
        Ey = np.array([]) if len(Ey_headers) != 1 else np.array(df[Ey_headers[0]]) 
        Ez = np.array([]) if len(Ez_headers) != 1 else np.array(df[Ez_headers[0]]) 
        grid = COMSOLGrid(df, unit, header_x, header_y, header_z) 
        return V, Ex, Ey, Ez, grid
class COMSOLElectrodeAdvanced(SimulatedElectrode): 
    def __init__(self, name, file, sim_unit=1e-6, 
                 header_x='% x', header_y='y', header_z='z', electrode_list=None, quantities = ['V', 'Ex', 'Ey', 'Ez'], **kwargs): 
        self.electrode_list = electrode_list
        self.quantities = quantities
        self.V, self.Ex, self.Ey, self.Ez, self.sim_grid = self.load_from_file(name, file, 
                                                                               unit=sim_unit, 
                                                                               header_x=header_x, 
                                                                               header_y=header_y, 
                                                                               header_z=header_z, **kwargs) 
    # Change the name based search for each electrode to a robust way by index. Now the electrode list will determine the order of the electrodes, so you need to make sure it matches the order of the electrodes in the COMSOL output file. 
    def load_from_file(self, electrode, file, excitation_prefix='V', sim_prefix='', skiprows=8, 
                       unit=1e-6, header_x='% x', header_y='y', header_z='z', **kwargs): 
    
        # Read CSV. Note: we don't use comment='%' because the header line itself starts with '% x'
        df = pd.read_csv(file, skiprows=skiprows) 
        
        # 1. Clean up coordinate headers (COMSOL leaves a '% ' in front of x)
        df.columns = [col.strip() for col in df.columns]
        
        # 2. Identify the spatial coordinate columns
        coord_cols = [header_x.strip(), header_y.strip(), header_z.strip()]
        data_cols = [col for col in df.columns if col not in coord_cols]
        
        # 3. Determine how many quantities are exported per electrode sweep
        num_quantities_per_sweep = len(self.quantities)  # e.g., 4 if ['V', 'Ex', 'Ey', 'Ez'] or 1 if ['V']
        
        if electrode != 'default':
            try:
                electrode_idx = self.electrode_list.index(electrode)
            except AttributeError:
                raise ValueError(f"Could not map {electrode}. Pass the ordered electrode list to map indices.")
                
            # Calculate exactly which chunk of columns belongs to this electrode sweep
            start_col_idx = electrode_idx * num_quantities_per_sweep
            electrode_cols = data_cols[start_col_idx : start_col_idx + num_quantities_per_sweep]
        else: 
            # Fallback for 'default' layout or single-sweep files
            electrode_cols = data_cols[:num_quantities_per_sweep]

        # 4. Map the columns strictly by position according to self.quantities
        # We build a dictionary mapping the quantity name to its corresponding data array
        data_dict = {}
        for idx, quantity in enumerate(self.quantities):
            if idx < len(electrode_cols):
                data_dict[quantity] = np.array(df[electrode_cols[idx]])
            else:
                data_dict[quantity] = np.array([])  # Fallback if file has fewer columns than expected

        # Extract arrays based on what was requested, defaulting to empty arrays if missing
        V = data_dict.get('V', np.array([]))
        Ex = data_dict.get('Ex', np.array([]))
        Ey = data_dict.get('Ey', np.array([]))
        Ez = data_dict.get('Ez', np.array([]))
        
        grid = COMSOLGrid(df, unit, header_x.strip(), header_y.strip(), header_z.strip()) 
        return V, Ex, Ey, Ez, grid

class COMSOLElectrodeRF(SimulatedElectrode): 
    def __init__(self, file, sim_unit=1e-6, 
                 header_x='% x', header_y='y', header_z='z', quantities=None, **kwargs): 
        
        # Default to exactly 3 quantities for the RF fields
        self.quantities = quantities if quantities is not None else ['Ex', 'Ey', 'Ez']
        
        # Load from file using simple positional mapping
        self.Ex, self.Ey, self.Ez, self.sim_grid = self.load_from_file(
            file, 
            unit=sim_unit, 
            header_x=header_x, 
            header_y=header_y, 
            header_z=header_z, 
            **kwargs
        ) 

    def load_from_file(self, file, skiprows=8, unit=1e-6, 
                       header_x='% x', header_y='y', header_z='z', **kwargs): 
    
        # Read CSV. We don't use comment='%' because the header line starts with '% x'
        df = pd.read_csv(file, skiprows=skiprows) 
        
        # 1. Clean up coordinate headers (COMSOL leaves a '% ' in front of x)
        df.columns = [col.strip() for col in df.columns]
        
        # 2. Identify the spatial coordinate columns
        coord_cols = [header_x.strip(), header_y.strip(), header_z.strip()]
        data_cols = [col for col in df.columns if col not in coord_cols]
        
        # 3. Determine how many quantities we expect to grab
        num_quantities = len(self.quantities)  # Typically 3 (Ex, Ey, Ez)
        
        # Grab exactly the first 'num_quantities' columns of physical data
        electrode_cols = data_cols[:num_quantities]

        # 4. Map columns strictly by position according to self.quantities
        data_dict = {}
        for idx, quantity in enumerate(self.quantities):
            if idx < len(electrode_cols):
                data_dict[quantity] = np.array(df[electrode_cols[idx]])
            else:
                data_dict[quantity] = np.array([])  # Fallback if file has fewer columns than expected

        # Extract arrays cleanly
        Ex = data_dict.get('Ex', np.array([]))
        Ey = data_dict.get('Ey', np.array([]))
        Ez = data_dict.get('Ez', np.array([]))
        
        grid = COMSOLGrid(df, unit, header_x.strip(), header_y.strip(), header_z.strip()) 
        return Ex, Ey, Ez, grid
    
class ANSYSElectrode(SimulatedElectrode): 
    def __init__(self, name, file, sim_unit=1e-6, 
                 header_x='x', header_y='y', header_z='z', **kwargs): 
    
        self.V, self.sim_grid = self.load_from_file(name, file, 
                                                    unit=sim_unit, 
                                                    header_x=header_x, 
                                                    header_y=header_y, 
                                                    header_z=header_z, **kwargs) 

    def load_from_file(self, name, file, skiprows=2, 
                       unit=1e-6, header_x='x', header_y='y', header_z='z', header_V='V',
                       **kwargs): 
        df = pd.read_csv(f'{file}_{name}.fld', sep=r' ', skiprows=skiprows, names=['x', 'y', 'z', '_', 'V'])
        V = np.array(df[header_x]), np.array(df[header_x]), np.array(df[header_x]), np.array(df[header_V])
        grid = COMSOLGrid(df, unit, header_x, header_y, header_z) 
        return V, grid