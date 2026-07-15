import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import csv
from .utils import eval_spherical_harmonics_by_term, get_r0_from_unit, get_Cj_fit, eval_spherical_harmonics, get_Cj_list, \
                   plot_all_potentials, plot_Mj, plot_V_DC, plot_potential_contours, \
                   plot_cutline_fits, compute_a, find_freq_shift, plot_potential_contours_slices
from .Grid import COMSOLGrid 
from .Electrode import COMSOLElectrode, ANSYSElectrode, COMSOLElectrodeRF, COMSOLElectrodeAdvanced
                   

class SimulatedTrap: 
    def __init__(self, result_file, electrodes, unit='um', L_ROI=50, quantities = ['V', 'Ex', 'Ey', 'Ez'], **kwargs): 
        self.constructed_V_total = False
        self.electrodes = {}
        self.r0 = get_r0_from_unit(unit)
        self.unit = unit
        self.electrode_list = electrodes
        self.quantities = quantities

    def get_V_matrix_ROI(self): 
        V = [] 
        for ei in self.electrodes: 
            # self.electrodes[ei].V_ROI = self.electrodes[ei].get_V_in_cube(L_cube=self.ROI_grid.L_cube)
            self.electrodes[ei].V_ROI = self.electrodes[ei].get_V_in_ROI(ROI_grid=self.ROI_grid)
            V.append(self.electrodes[ei].V_ROI) 
        return np.array(V).T
    def get_V_matrix(self):
        V = [] 
        for ei in self.electrodes: 
            V.append(self.electrodes[ei].V) 
        return np.array(V).T

    def get_electrode_voltages(self, C=0, Ey=0, Ez=0, Ex=0, 
                               U3=0, U4=0, U2=-1, U5=0, U1=0, **kwargs): 
        Cj_ideal = get_Cj_list(C=C, Ey=Ey, Ez=Ez, Ex=Ex, 
                               U3=U3, U4=U4, U2=U2, U5=U5, U1=U1, **kwargs)
        x, y, z = self.ROI_grid.get_xyz_array() 
        V_ideal = eval_spherical_harmonics(Cj_ideal, x, y, z) 
        V_DC = np.linalg.lstsq(self.V_matrix_ROI, V_ideal, rcond=None)[0] 
        return V_DC 
    
    def get_electrode_voltages_C_ignored_full(self, C=0, Ey=0, Ez=0, Ex=0, 
                               U3=0, U4=0, U2=-1, U5=0, U1=0, **kwargs): 
        Cj_ideal = get_Cj_list(C=C, Ey=Ey, Ez=Ez, Ex=Ex, 
                               U3=U3, U4=U4, U2=U2, U5=U5, U1=U1, **kwargs)
        x, y, z = self.sim_grid.get_xyz_array() 
        V_ideal = eval_spherical_harmonics(Cj_ideal, x, y, z) 
        
        # 1. Flatten V_ideal if it isn't already 1D to match the rows of the matrix
        V_ideal_flat = V_ideal.ravel()

        # 2. Append a column of ones to the matrix to represent the constant offset
        #    This creates a temporary matrix: [V_matrix_ROI | 1]
        ones_column = np.ones((self.V_matrix.shape[0], 1))
        A_augmented = np.hstack([self.V_matrix, ones_column])
        
        # 3. Perform the least-squares fit on the augmented matrix
        res = np.linalg.lstsq(A_augmented, V_ideal_flat, rcond=None)[0] 
        
        # 4. Slice off the last element (which is the redundant constant offset)
        #    and return just the electrode voltages
        V_DC = res[:-1]
        
        return V_DC
    
    def get_electrode_voltages_C_ignored(self, C=0, Ey=0, Ez=0, Ex=0, 
                               U3=0, U4=0, U2=-1, U5=0, U1=0, **kwargs): 
        Cj_ideal = get_Cj_list(C=C, Ey=Ey, Ez=Ez, Ex=Ex, 
                               U3=U3, U4=U4, U2=U2, U5=U5, U1=U1, **kwargs)
        x, y, z = self.ROI_grid.get_xyz_array() 
        V_ideal = eval_spherical_harmonics(Cj_ideal, x, y, z) 
        
        # 1. Flatten V_ideal if it isn't already 1D to match the rows of the matrix
        V_ideal_flat = V_ideal.ravel()

        # 2. Append a column of ones to the matrix to represent the constant offset
        #    This creates a temporary matrix: [V_matrix_ROI | 1]
        ones_column = np.ones((self.V_matrix_ROI.shape[0], 1))
        A_augmented = np.hstack([self.V_matrix_ROI, ones_column])
        
        # 3. Perform the least-squares fit on the augmented matrix
        res = np.linalg.lstsq(A_augmented, V_ideal_flat, rcond=None)[0] 
        
        # 4. Slice off the last element (which is the redundant constant offset)
        #    and return just the electrode voltages
        V_DC = res[:-1]
        
        return V_DC

    def construct_V_total(self, C=0, Ey=0, Ez=0, Ex=0, 
                          U3=0, U4=0, U2=-1, U5=0, U1=0, original = False, **kwargs): 
        if original:
            V_DC = np.ones(len(self.electrodes))
            V_total = np.dot(self.V_matrix_ROI, V_DC) 
            self.V_DC = V_DC
            self.V_total = V_total 
            self.constructed_V_total = True
        else:
            V_DC = self.get_electrode_voltages(C=C, Ey=Ey, Ez=Ez, Ex=Ex, 
                                            U3=U3, U4=U4, U2=U2, U5=U5, U1=U1, **kwargs)
            V_total = np.dot(self.V_matrix_ROI, V_DC) 
            self.V_DC = V_DC
            self.V_total = V_total 
            self.constructed_V_total = True
    def construct_V_total_full_range(self, C=0, Ey=0, Ez=0, Ex=0, 
                          U3=0, U4=0, U2=-1, U5=0, U1=0, original = False, **kwargs):
        if original:
            V_DC = np.ones(len(self.electrodes))
            V_total = np.dot(self.V_matrix, V_DC) 
            self.V_DC_full = V_DC
            self.V_total_full = V_total 
            self.constructed_V_total_full = True
        else:
            V_DC = self.get_electrode_voltages_C_ignored_full(C=C, Ey=Ey, Ez=Ez, Ex=Ex, 
                                            U3=U3, U4=U4, U2=U2, U5=U5, U1=U1, **kwargs)
            V_total = np.dot(self.V_matrix, V_DC) 
            self.V_DC_full = V_DC
            self.V_total_full = V_total 
            self.constructed_V_total_full = True

    def expand_spherical_harmonics(self, order=2): 
        assert self.constructed_V_total, "Construct the desired total potential" + \
                                         "using self.construct_V_total(...) first"
        x, y, z = self.ROI_grid.get_xyz_array() 
        if self.RF_file is not None:
            x = x - self.RF_x0
            y = y - self.RF_y0
            z = z - self.sim_grid.z0 # assume linear trap with axial as z axis.
        else:
            x = x - self.sim_grid.x0
            y = y - self.sim_grid.y0
            z = z - self.sim_grid.z0
        self.Cj_fit = get_Cj_fit(self.V_total, x, y, z, order=order)
        self.V_fit = eval_spherical_harmonics(self.Cj_fit, x, y, z)
    def expand_spherical_harmonics_electrode(self, ei, order =2):
        x, y, z = self.ROI_grid.get_xyz_array() 
        if self.RF_file is not None:
            x = x - self.RF_x0
            y = y - self.RF_y0
            z = z - self.sim_grid.z0 # assume linear trap with axial as z axis.
        else:
            x = x - self.sim_grid.x0
            y = y - self.sim_grid.y0
            z = z - self.sim_grid.z0
        self.electrodes[ei].Cj_fit = get_Cj_fit(self.electrodes[ei].V_ROI, x, y, z, order=order)
        self.electrodes[ei].V_fit = eval_spherical_harmonics(self.electrodes[ei].Cj_fit, x, y, z)
    def expand_spherical_harmonics_all_electrodes(self, order =2):
        for ei in self.electrodes:
            self.expand_spherical_harmonics_electrode(ei, order=order)


    def plot_V_fit(self, m=1, plot_scale=1): 
        assert hasattr(self, "Cj_fit"), "Construct spherical harmonics expansion" + \
                                        "using self.expand_spherical_harmonics(...) first"
        x, y, z = self.ROI_grid.get_xyz_array() 
        plot_all_potentials((x, y, z), self.V_total, self.V_fit, m=1, plot_scale=plot_scale, unit=self.unit)

    def plot_Mj(self, Mj_threshold=0.01, logy=True, title='', save_fig=False): 
        Cj_fit = abs(self.Cj_fit) if logy else self.Cj_fit
        ylabel = f'$|M_j|$ (1/{self.unit}$^l$)' if logy else f'$M_j$ (1/{self.unit}$^l$)'
        plot_Mj(Cj_fit, Mj_threshold=Mj_threshold, title=title, 
                ylabel=ylabel, save_fig=save_fig, logy=logy)

    def plot_V_DC(self):
        plot_V_DC(self.V_DC, self.electrodes.keys(), U2=self.Cj_fit[6], unit=self.unit)

    def plot_potential_contours(self):  
        x, y, z = self.ROI_grid.get_xyz_array() 
        plot_potential_contours(self.V_total, x, y, z)

    def plot_potential_contours_slices(self):  
        x, y, z = self.ROI_grid.get_xyz_array() 
        plot_potential_contours_slices(self.V_total, x, y, z, unit=self.unit)

    def plot_potential_contours_full(self):
        x, y, z = self.sim_grid.get_xyz_array()
        plot_potential_contours_slices(self.V_total_full, x, y, z, unit=self.unit)

    def plot_cutline_fits(self):
        x, y, z = self.ROI_grid.get_xyz_array()
        self.cutline_fit_coeff = plot_cutline_fits(self.V_total, x, y, z, *self.ROI_grid.get_grid_center(), unit=self.unit)

    def plot_estimated_frequency_shift(self, Amin=0, Amax=100, logx=True): 
        if not hasattr(self, "cutline_fit_coeff"): 
            self.plot_cutline_fits()
        fig, ax = plt.subplots() 
        A = np.linspace(Amin, Amax, 1000)
        a_x = compute_a(self.cutline_fit_coeff['x'])
        a_y = compute_a(self.cutline_fit_coeff['y']) 
        a_z = compute_a(self.cutline_fit_coeff['z'])
        x_shift = find_freq_shift(A, a_x)
        y_shift = find_freq_shift(A, a_y)
        z_shift = find_freq_shift(A, a_z)
        ax.plot(A, x_shift, label='x') 
        ax.plot(A, y_shift, label='y') 
        ax.plot(A, z_shift, label='z') 
        ax.set_xlabel(f'Amplitude ({self.unit})')
        ax.set_ylabel(r'$|\Delta \omega / \omega|$')
        ax.grid() 
        ax.set_yscale('log')
        if logx: 
            ax.set_xscale('log')
        ax.legend()
        plt.tight_layout()
        plt.show()
    def get_RF_null_pos(self):
      idx = np.argmin(np.abs(self.RF_field.Ex)**2 + np.abs(self.RF_field.Ey)**2 + np.abs(self.RF_field.Ez)**2)
      self.RF_x0, self.RF_y0, self.RF_z0 = self.sim_grid_RF.x[idx], self.sim_grid_RF.y[idx], self.sim_grid_RF.z[idx]

    def get_cfile(self, output_file, fused_electrodes = None):
        
        self.expand_spherical_harmonics_all_electrodes()
        Cj_fit_all = []
        for ei in self.electrode_list:
            Cj_fit_all.append(self.electrodes[ei].Cj_fit[1:])  # Exclude the constant term (C0)
        self.cfile = np.linalg.pinv(np.column_stack(Cj_fit_all))

        # Define headers
        headers = ['', 'Ey', 'Ez', 'Ex', 'U3', 'U4', 'U2', 'U5', 'U1']

        with open(output_file, mode='w', newline='') as f:
            writer = csv.writer(f)
            
            # Write the header row
            writer.writerow(headers)
            
            # Write each electrode row alongside its corresponding cfile row data
            for label, data_row in zip(self.electrode_list, self.cfile):
                # [label] + list(data_row) puts the electrode name in the first empty column
                writer.writerow([label] + list(data_row))
        return self.cfile
    
    def get_cfile(self, output_file, fused_electrodes = None):
        self.expand_spherical_harmonics_all_electrodes()
        
        # 1. Gather all individual electrode fit vectors (excluding C0)
        Cj_fit_all = []
        for ei in self.electrode_list:
            Cj_fit_all.append(self.electrodes[ei].Cj_fit[1:])
        
        # C_matrix shape: (num_features, num_electrodes)
        C_matrix = np.column_stack(Cj_fit_all)
        
        # 2. Process fused_electrodes into a standard list of lists
        groups = []
        if fused_electrodes is not None:
            # If a single list of strings is provided: ['tr1', 'tl1'] -> [['tr1', 'tl1']]
            if isinstance(fused_electrodes, list) and len(fused_electrodes) > 0:
                if isinstance(fused_electrodes[0], str):
                    groups = [fused_electrodes]
                else:
                    groups = fused_electrodes

        # Map electrode labels to their current column indices in self.electrode_list
        label_to_idx = {label: i for i, label in enumerate(self.electrode_list)}
        
        # Track which indices have been grouped/fused to avoid duplicates
        fused_indices_set = set()
        valid_groups = []
        for gp in groups:
            # Filter group to make sure the electrodes actually exist in self.electrode_list
            valid_gp_indices = [label_to_idx[name] for name in gp if name in label_to_idx]
            if len(valid_gp_indices) > 1:
                valid_groups.append(valid_gp_indices)
                fused_indices_set.update(valid_gp_indices)

        # 3. Construct the reduced fitting matrix
        # Columns of non-fused electrodes stay as they are.
        # Columns of fused groups are summed together.
        reduced_columns = []
        # This list will keep track of which group or original index each column in our reduced matrix maps to
        column_mapping = [] 
        
        # First, add the unfused electrodes
        for idx in range(len(self.electrode_list)):
            if idx not in fused_indices_set:
                reduced_columns.append(C_matrix[:, idx])
                column_mapping.append([idx]) # Single-item list represents its own index
                
        # Next, add the fused groups as single merged columns
        for gp_indices in valid_groups:
            merged_col = np.sum(C_matrix[:, gp_indices], axis=1)
            reduced_columns.append(merged_col)
            column_mapping.append(gp_indices)

        # Stack them back to a matrix
        reduced_C_matrix = np.column_stack(reduced_columns)
        
        # 4. Compute the pseudoinverse of the reduced system
        # pinv_reduced shape: (num_reduced_electrodes, num_features)
        pinv_reduced = np.linalg.pinv(reduced_C_matrix)
        
        # 5. Map the reduced pseudoinverse back to the original full-sized electrode layout
        # Initialize cfile to zeros: shape (num_electrodes, num_features)
        self.cfile = np.zeros((len(self.electrode_list), pinv_reduced.shape[1]))
        
        for reduced_idx, original_indices in enumerate(column_mapping):
            row_data = pinv_reduced[reduced_idx, :]
            for orig_idx in original_indices:
                # All physical electrodes in a fused group get the exact same row data (voltage coefficient mapping)
                self.cfile[orig_idx, :] = row_data

        # Define headers
        headers = ['', 'Ey', 'Ez', 'Ex', 'U3', 'U4', 'U2', 'U5', 'U1']

        # 6. Write out results to the CSV
        with open(output_file, mode='w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            
            for label, data_row in zip(self.electrode_list, self.cfile):
                writer.writerow([label] + list(data_row))
                
        return self.cfile

    
    def get_cfile_from_matrix(self, output_file):
        # this code is for some case where the trap is very anharmonic.
        voltages_all = []

        # 1. Define the parameters you want to iterate over in order
        param_names = ['Ey', 'Ez', 'Ex', 'U3', 'U4', 'U2', 'U5', 'U1']

        # 2. Iterate through the parameters
        for hot_param in param_names:
            # Dynamically build a dictionary where the current param is 1, others are 0
            # (Note: U2 defaults to -1 in your signature, so we explicitly set it to 0 here)
            one_hot_args = {name: (1 if name == hot_param else 0) for name in param_names}
            
            # 3. Call your method passing the dictionary as keyword arguments
            voltages = self.get_electrode_voltages(C=0, **one_hot_args)
            voltages_all.append(voltages)
        self.cfile = np.column_stack(voltages_all)  # Stack the results into a matrix
        
        # Define headers
        headers = ['', 'Ey', 'Ez', 'Ex', 'U3', 'U4', 'U2', 'U5', 'U1']

        with open(output_file, mode='w', newline='') as f:
            writer = csv.writer(f)
            
            # Write the header row
            writer.writerow(headers)
            
            # Write each electrode row alongside its corresponding cfile row data
            for label, data_row in zip(self.electrode_list, self.cfile):
                # [label] + list(data_row) puts the electrode name in the first empty column
                writer.writerow([label] + list(data_row))

        return self.cfile
    def reconstructROI(self, L_ROI):
        if self.RF_file is not None:
            self.ROI_grid = self.sim_grid.gen_subcube(L_cube=L_ROI, center=(self.RF_x0, self.RF_y0, self.RF_z0))
        else:
            self.ROI_grid = self.sim_grid.gen_subcube(L_cube=L_ROI) 


class COMSOLTrap(SimulatedTrap): 
    def __init__(self, result_file, electrodes, unit='um', L_ROI=float('inf'), skiprows=8, **kwargs): 
        super().__init__(result_file, electrodes, unit, L_ROI)
        self.sim_grid = COMSOLGrid(pd.read_csv(result_file, skiprows=skiprows), **kwargs)  # Grid used in COMSOL simulation
        self.sim_grid.scale_xyz(self.r0) 
        self.ROI_grid = self.sim_grid.gen_subcube(L_cube=L_ROI) 
        if len(electrodes) == 0:
            self.electrodes['default'] = COMSOLElectrode('default', result_file, **kwargs)
        else:
            for ei in electrodes: 
                self.electrodes[ei] = COMSOLElectrode(ei, result_file, **kwargs) 
                self.electrodes[ei].set_sim_grid(self.sim_grid)
        self.V_matrix_ROI = self.get_V_matrix_ROI() 

class COMSOLTrapAdvanced(SimulatedTrap):
    def __init__(self, result_file, electrodes, quantities = ['V', 'Ex', 'Ey', 'Ez'], unit='um', L_ROI=float('inf'), skiprows=8, RF_file = None, **kwargs): 
        super().__init__(result_file, electrodes, unit, L_ROI)
        self.sim_grid = COMSOLGrid(pd.read_csv(result_file, skiprows=skiprows), **kwargs)  # Grid used in COMSOL simulation
        self.sim_grid.scale_xyz(self.r0) 
        self.RF_file = RF_file
        if RF_file is not None:
            self.sim_grid_RF= COMSOLGrid(pd.read_csv(RF_file, skiprows=skiprows), **kwargs)
            self.sim_grid_RF.scale_xyz(self.r0)
            self.RF_field = COMSOLElectrodeRF(RF_file, quantities=['Ex', 'Ey', 'Ez'], **kwargs)
            self.get_RF_null_pos()
            self.ROI_grid = self.sim_grid.gen_subcube(L_cube=L_ROI, center=(self.RF_x0, self.RF_y0, self.RF_z0))
        else:
            self.ROI_grid = self.sim_grid.gen_subcube(L_cube=L_ROI) 
        if len(electrodes) == 0:
            self.electrodes['default'] = COMSOLElectrodeAdvanced('default', result_file, **kwargs)
        else:
            for ei in electrodes: 
                self.electrodes[ei] = COMSOLElectrodeAdvanced(ei, result_file, electrode_list = electrodes, quantities = quantities, **kwargs) 
                self.electrodes[ei].set_sim_grid(self.sim_grid)
        self.V_matrix_ROI = self.get_V_matrix_ROI() 
        self.V_matrix = self.get_V_matrix()
        
    

class ANSYSTrap(SimulatedTrap): 
    def __init__(self, result_file, electrodes, unit='um', L_ROI=float('inf'), skiprows=2, **kwargs): 
        super().__init__(result_file, electrodes, unit, L_ROI)
        self.sim_grid = COMSOLGrid(pd.read_csv(result_file, skiprows=skiprows), **kwargs)  # Grid used in COMSOL simulation
        self.sim_grid.scale_xyz(self.r0) 
        self.ROI_grid = self.sim_grid.gen_subcube(L_cube=L_ROI) 
        for ei in electrodes: 
            self.electrodes[ei] = ANSYSElectrode(ei, result_file, **kwargs) 
            self.electrodes[ei].set_sim_grid(self.sim_grid)
        
        self.V_matrix_ROI = self.get_V_matrix_ROI() 