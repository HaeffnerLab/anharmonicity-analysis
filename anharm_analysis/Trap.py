import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from .utils import get_r0_from_unit, get_Cj_fit, eval_spherical_harmonics, get_Cj_list, \
                   plot_all_potentials, plot_Mj, plot_V_DC, plot_potential_contours, \
                   plot_cutline_fits, compute_a, find_freq_shift
from .Grid import COMSOLGrid 
from .Electrode import COMSOLElectrode
                   

class SimulatedTrap: 
    def __init__(self, result_file, electrodes, unit='um', L_ROI=50): 
        self.constructed_V_total = False
        self.electrodes = {}
        self.r0 = get_r0_from_unit(unit)
        self.unit = unit

    def get_V_matrix_ROI(self): 
        V = [] 
        for ei in self.electrodes: 
            V.append(self.electrodes[ei].get_V_in_cube(L_cube=self.ROI_grid.L_cube)) 
        return np.array(V).T

    def get_electrode_voltages(self, C=0, Ey=0, Ez=0, Ex=0, 
                               U3=0, U4=0, U2=-1, U5=0, U1=0, **kwargs): 
        Cj_ideal = get_Cj_list(C=C, Ey=Ey, Ez=Ez, Ex=Ex, 
                               U3=U3, U4=U4, U2=U2, U5=U5, U1=U1, **kwargs)
        x, y, z = self.ROI_grid.get_xyz_array() 
        V_ideal = eval_spherical_harmonics(Cj_ideal, x, y, z) 
        V_DC = np.linalg.lstsq(self.V_matrix_ROI, V_ideal, rcond=None)[0] 
        return V_DC 

    def construct_V_total(self, C=0, Ey=0, Ez=0, Ex=0, 
                          U3=0, U4=0, U2=-1, U5=0, U1=0, **kwargs): 
        V_DC = self.get_electrode_voltages(C=C, Ey=Ey, Ez=Ez, Ex=Ex, 
                                           U3=U3, U4=U4, U2=U2, U5=U5, U1=U1, **kwargs)
        V_total = np.dot(self.V_matrix_ROI, V_DC) 
        self.V_DC = V_DC
        self.V_total = V_total 
        self.constructed_V_total = True

    def expand_spherical_harmonics(self, order=2): 
        assert self.constructed_V_total, "Construct the desired total potential" + \
                                         "using self.construct_V_total(...) first"
        x, y, z = self.ROI_grid.get_xyz_array() 
        self.Cj_fit = get_Cj_fit(self.V_total, x, y, z, order=order)
        self.V_fit = eval_spherical_harmonics(self.Cj_fit, x, y, z)

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

class COMSOLTrap(SimulatedTrap): 
    def __init__(self, result_file, electrodes, unit='um', L_ROI=50, skiprows=8, **kwargs): 
        super().__init__(result_file, electrodes, unit, L_ROI)
        self.sim_grid = COMSOLGrid(pd.read_csv(result_file, skiprows=skiprows), **kwargs)  # Grid used in COMSOL simulation
        self.sim_grid.scale_xyz(self.r0) 
        self.ROI_grid = self.sim_grid.gen_subcube(L_cube=L_ROI) 
        for ei in electrodes: 
            self.electrodes[ei] = COMSOLElectrode(ei, result_file, **kwargs) 
            self.electrodes[ei].set_sim_grid(self.sim_grid)
        self.V_matrix_ROI = self.get_V_matrix_ROI() 