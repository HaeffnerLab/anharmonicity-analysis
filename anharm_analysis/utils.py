import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from juliacall import Main as jl
from scipy.interpolate import griddata
import os 

if not os.path.exists("../.setup_log"): 
    jl.seval("using Pkg")
    jl.Pkg.add("SphericalHarmonicExpansions")
    os.mkdir("../.setup_log")
    
jl.seval("using SphericalHarmonicExpansions") 


def big_plt_font(): 
    plt.rcParams.update({'font.size': 14,
                     'lines.markersize': 12,
                     'lines.linewidth': 2.5,
                     'xtick.labelsize': 15,
                     'ytick.labelsize': 15,
                     'errorbar.capsize': 2}) 
    
def plot_3D_potential(coord, Phi, size=1, n=100, scale=1,
                   cmap='viridis', title=r'$\Phi(x,y,z)$',
                   ax=None):
    x, y, z = coord
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Create scatter plot
    sc = ax.scatter(x[::n]*scale, y[::n]*scale, z[::n]*scale, 
                    s=size, c=Phi[::n], 
                    marker='.', cmap=cmap)
    
    # Add color bar
    cb = plt.colorbar(sc, ax=ax, shrink=0.5, aspect=10)
    cb.set_label('Potential (V)')
    
    # Labels and display
    ax.set_xlabel('x (um)')
    ax.set_ylabel('y (um)')
    ax.set_zlabel('z (um)')
    
    plt.show()

def plot_potential(coord, Phi, size=1, n=100, scale=1, unit='um',
                   cmap='viridis', title=r'$\Phi(x,y,z)$',
                   ax=None):
    x, y, z = coord
    if ax is None:
        fig = plt.figure(figsize=(6,6))
        ax = plt.axes(projection='3d')
    im = ax.scatter(x[::n]*scale, y[::n]*scale, z[::n]*scale, 
                    s=size, c=Phi[::n], 
                    marker='.', cmap=cmap)
    ax.set_xlabel(f'x ({unit})')
    ax.set_ylabel(f'y ({unit})')
    ax.set_zlabel(f'z ({unit})')
    ax.set_title(title)
    if ax is None:
        fig.colorbar(im, ax=ax, shrink=0.8)
        plt.tight_layout()
        plt.show()
    else:
        return im
        
def plot_all_potentials(plot_coord_fit, V, Phi, m, unit='um', plot_scale=1):
    """
    Plot the original, fitted, and residual potential.
    """
    fig1 = plt.figure(figsize=plt.figaspect(1/3))
    ax1 = fig1.add_subplot(1, 3, 1, projection='3d')
    ax2 = fig1.add_subplot(1, 3, 2, projection='3d')
    ax3 = fig1.add_subplot(1, 3, 3, projection='3d')
    
    im1 = plot_potential(plot_coord_fit, V.flatten(), size=10, n=m, scale=plot_scale, 
                   unit=unit, title=r'Constructed $\Phi(x,y,z)$', ax=ax1)
    im2 = plot_potential(plot_coord_fit, Phi.flatten(), size=10, n=m, scale=plot_scale, 
                   unit=unit, title=r'Fit $\Phi(x,y,z)$', ax=ax2)
    im3 = plot_potential(plot_coord_fit, abs(Phi.flatten()-V.flatten()), size=10, n=m, scale=plot_scale, 
                   unit=unit, title=r'Residual $\Delta\Phi(x,y,z)$', cmap='Reds', ax=ax3)
    fig1.colorbar(im1, ax=ax1, shrink=0.8)
    fig1.colorbar(im2, ax=ax2, shrink=0.8)
    fig1.colorbar(im3, ax=ax3, shrink=0.8)
    plt.tight_layout()
    plt.show()

def add_value_labels(ax, spacing=0.1, threshold=0.01, decimal=2):
    """Add labels to the end of each bar in a bar chart.

    Arguments:
        ax (matplotlib.axes.Axes): The matplotlib object containing the axes
            of the plot to annotate.
        spacing (int): The distance between the labels and the bars.
    """

    # For each bar: Place a label
    for rect in ax.patches:
        # Get X and Y placement of label from rect.
        y_value = rect.get_height()
        if abs(y_value) < threshold:
            continue
        x_value = rect.get_x() + rect.get_width() / 2

        # Number of points between bar and label. Change to your liking.
        space = spacing
        # Vertical alignment for positive values
        va = 'bottom'

        # If value of bar is negative: Place label below bar
        if y_value < 0:
            # Invert space to place label below
            space *= -1
            # Vertically align label at top
            va = 'top'

        # Use Y value as label and format number with one decimal place
        if decimal == 3:
            label = "{:.3f}".format(y_value)
        else:
            label = "{:.2f}".format(y_value)

        # Create annotation
        ax.annotate(
            label,                      # Use `label` as label
            (x_value, y_value),         # Place label at end of the bar
            xytext=(0, space),          # Vertically shift label by `space`
            textcoords="offset points", # Interpret `xytext` as offset in points
            ha='center',                # Horizontally center label
            va=va)                      # Vertically align label differently for
                                        # positive and negative values.
        
def plot_Mj(Mj, mutipole_names=['C', 'Ey', 'Ez', 'Ex', 'U3', 'U4', 'U2', 'U5', 'U1'], 
            Mj_threshold=0.01, xlabel=r'$j$', ylabel=r'$M_j$',
            title='', save_fig=False, logy=False):
    
    fig, ax = plt.subplots(figsize=(max(6,0.3*len(Mj)), 4))
    ax.bar(list(range(1,len(Mj)+1)), Mj.flatten())
    add_value_labels(ax, threshold=Mj_threshold)
    #ax.axvline(np.argmax(abs(Mj))+1, label='j ='+str(np.argmax(abs(Mj))+1), 
    #            linestyle='--', color='r', alpha=0.7)
    tick_name = list(mutipole_names)
    tick_name += list(range(len(tick_name)+1, len(Mj)+1))
    df = pd.DataFrame({'Mj': [f'{float(i):.3e}' for i in Mj]}, index=tick_name).transpose()
    try:
        display(df)
    except:
        print(df)
    ax.set_xticks(range(1,len(Mj)+1), tick_name, rotation = -90)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if logy: 
        ax.set_yscale('log')
    ax.grid()
    ax.set_title(title)
    plt.tight_layout()
    if save_fig:
        plt.savefig(f'{title}.pdf')
    plt.show() 

def plot_V_DC(V_DC, electrode_names, U2=-1, unit='um'):
    fig, ax = plt.subplots(figsize=(max(6,0.3*len(V_DC)), 4))
    ax.bar(list(range(1,len(V_DC)+1)), V_DC.flatten())
    add_value_labels(ax, threshold=0, decimal=3)
    tick_name = list(electrode_names)
    ax.set_xticks(range(1,len(electrode_names)+1), tick_name)
    ax.set_xlabel('Electrode')
    ax.set_ylabel(f'Voltage (V) [at $U_2={U2:.1f}$ V/{unit}$^2$]')
    ax.grid()
    plt.tight_layout()
    plt.show() 

def plot_potential_contours(V, x, y, z, unit='um'):
    """
    Plot contour plots of the potential V in the xy, yz, and xz planes.

    Parameters:
        V (ndarray): 1D array of strictly positive potential values.
        x, y, z (ndarray): 1D coordinate arrays, same shape as V.
    """
    assert V.ndim == x.ndim == y.ndim == z.ndim == 1, "All inputs must be 1D arrays"
    assert len(V) == len(x) == len(y) == len(z), "All arrays must be the same length"

    grid_res = 100
    cmap = 'seismic'

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)

    planes = [
        ('xy', x, y, f'x ({unit})', f'y ({unit})'),
        ('yz', y, z, f'y ({unit})', f'z ({unit})'),
        ('xz', x, z, f'x ({unit})', f'z ({unit})'),
    ]

    for ax, (label, coord1, coord2, xlabel, ylabel) in zip(axes, planes):
        xi = np.linspace(min(coord1), max(coord1), grid_res)
        yi = np.linspace(min(coord2), max(coord2), grid_res)
        X, Y = np.meshgrid(xi, yi)
        points = np.column_stack((coord1, coord2))
        Z = griddata(points, V, (X, Y), method='linear')

        contour = ax.contourf(X, Y, Z, levels=50, cmap=cmap)
        ax.set_title(f'{label}-plane')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    cbar = fig.colorbar(contour, ax=axes.ravel().tolist(), orientation='vertical')
    cbar.set_label('Potential (V)')

    plt.show()

def plot_cutline_fits(V, x, y, z, x0=0.0, y0=0.0, z0=0.0, tol=1e-6, unit='um'):
    """
    Plot 1D cutlines through the point (x0, y0, z0) along x, y, z directions.
    Fits quadratic and 8th-order polynomials to each and plots them.

    Parameters:
        V, x, y, z (ndarray): 1D arrays of potential and positions.
        x0, y0, z0 (float): The reference point through which cutlines pass.
        tol (float): Tolerance to select axis-aligned points.

    Returns:
        dict: Dictionary containing 8th-order polynomial coefficients for each direction.
    """
    directions = {
        'x': (x, np.where((np.abs(y - y0) < tol) & (np.abs(z - z0) < tol)), x0, f'x ({unit})'),
        'y': (y, np.where((np.abs(x - x0) < tol) & (np.abs(z - z0) < tol)), y0, f'y ({unit})'),
        'z': (z, np.where((np.abs(x - x0) < tol) & (np.abs(y - y0) < tol)), z0, f'z ({unit})'),
    }

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    coeffs = {}

    for ax, (label, (coord, mask, c0, coord_label)) in zip(axes, directions.items()):
        coord_cut = coord[mask]
        V_cut = V[mask]

        if len(coord_cut) < 10:
            ax.set_title(f'Insufficient data for {label}-cut')
            continue

        # Sort for clean plotting
        sorted_idx = np.argsort(coord_cut)
        coord_cut = coord_cut[sorted_idx]
        V_cut = V_cut[sorted_idx]

        # Plot original data
        ax.plot(coord_cut, V_cut, 'o', label='Data')
        coord_plot = np.linspace(min(coord_cut), max(coord_cut), 1000)

        # Fit quadratic
        p2 = np.polynomial.polynomial.Polynomial.fit(coord_cut, V_cut, 2)
        V_fit2 = p2(coord_plot)   #p.polyval(p2, coord_cut)
        ax.plot(coord_plot, V_fit2, 'r--', label='Quadratic Fit')

        # Fit 8th-order polynomial
        p8 = np.polynomial.polynomial.Polynomial.fit(coord_cut, V_cut, 8)
        V_fit8 = p8(coord_plot) #np.polyval(p8, coord_cut)
        ax.plot(coord_plot, V_fit8, 'k-', label='8th-Order Fit')

        coeffs[label] = p8.convert().coef
        coeffs[f'{label}2'] = p2.convert().coef

        ax.set_title(f'Cutline along {label}-axis through trap center')
        ax.set_xlabel(coord_label)
        ax.set_ylabel('Potential (V)')
        ax.legend()
        ax.grid()

    plt.tight_layout()
    plt.show()
    
    return coeffs
    
def eval_spherical_harmonics(C, x, y, z):
    jl.seval(f'C = {list(C)}')
    jl.seval("c = SphericalHarmonicCoefficients(C)")
    jl.seval("@polyvar x y z")
    jl.seval("f = sphericalHarmonicsExpansion(c, x, y, z)")
    fastf = jl.seval("(x,y,z) -> fastfunc(f)([x,y,z])")
    result = jl.broadcast(fastf, x, y, z) 
    return np.array(result)

def eval_spherical_harmonics_by_term(x, y, z, order=2): 
    N = (order+1)**2 
    V_list = []
    for i in range(N): 
        C = np.zeros(N) 
        C[i] = 1
        V_list.append(eval_spherical_harmonics(C, x, y, z)) 
    return np.array(V_list).T
    
def get_Cj_list(C=0, Ey=0, Ez=0, Ex=0, U3=0, U4=0, U2=-1, U5=0, U1=0, **kwargs):
    multipole_coeffs = [C, Ey, Ez, Ex, U3, U4, U2, U5, U1]
    L = 2
    for mj in kwargs: 
        term = int(mj[1:])  
        while (L+1)**2 < term: 
            L += 1 
    
    N_terms = (L+1)**2 
    C = np.zeros(N_terms)
    for i in range(N_terms): 
        if i < 9: 
            C[i] = multipole_coeffs[i] 
        elif f'm{i}' in kwargs: 
            C[i] = kwargs[f'm{i}'] 
    return C

def get_Cj_fit(V, x, y, z, order=2): 
    V_spherical_harmonics_matrix = eval_spherical_harmonics_by_term(x, y, z, order=order) 
    return np.linalg.lstsq(V_spherical_harmonics_matrix, V, rcond=None)[0]

def get_r0_from_unit(unit): 
    unit_dict = {'pm': 1e-12, 'nm': 1e-9, 'um': 1e-6, 'mm': 1e-3, 'cm': 1e-2, 'm': 1}
    return unit_dict[unit]

def compute_a(C0):
    """
    Compute the coefficients a_i (i=0...len(C)-1) given C_i values.

    C should be a sequence (list, tuple, etc.) such that C[i] == C_i.
    Returns a list a of the same length, with
      a[2] = -15*(C3)**2/16 + 3*C4/4
      a[3] = C3 * a[2]
      a[4] = ...
      …
      a[7] = ...
    and a[i]=0 for any i<2 or i>7 (or if C is too short).
    """
    C = C0 / C0[2]
    N = len(C)
    # initialize all a[i]=0
    a = np.zeros(N)

    # a2
    if N > 4:
        a[2] = -15 * C[3]**2 / 16 +  3 * C[4] / 4

    # a3 = C3 * a2
    if N > 3:
        a[3] = C[3] * a[2]

    # a4
    if N > 6:
        a[4] = (
            -2565 * C[3]**4  / 1024
            +  645 * C[3]**2 * C[4] / 128
            -   21 * C[4]**2     /  64
            -  105 * C[3] * C[5] /  32
            +   15 * C[6]        /  16
        )

    # a5
    if N > 6:
        a[5] = (
            -2565 * C[3]**5   / 512
            +  765 * C[3]**3 * C[4]  /  64
            -   69 * C[3] * C[4]**2  /  32
            -   15 * C[3]**2 * C[5]  /   2
            +    3 * C[4] * C[5]     /   4
            +   15 * C[3] * C[6]     /   8
        )
        # (equivalently you can use the simplified form
        #  a[5] = (C[5] - 2*C[3]*C[4]) * a[2] + 2*C[3] * a[4])

    # a6
    if N > 8:
        a[6] = (
            -205845 * C[3]**6        /  16384
            + 159795 * C[3]**4 * C[4] /   4096
            -  21039 * C[3]**2 * C[4]**2 / 1024
            +     81 * C[4]**3        /   256
            -  13545 * C[3]**3 * C[5] /   512
            +   1995 * C[3] * C[4] * C[5] / 128
            -    315 * C[5]**2       /   128
            +   3015 * C[3]**2 * C[6] /   256
            -     57 * C[4] * C[6]   /    64
            -    315 * C[3] * C[7]   /    64
            +     35 * C[8]          /    32
        )

    # a7
    if N > 7:
        a[7] = (
            3 * C[3] * a[6]
            + (-3*C[3]**3 - 4*C[3]*C[4] + 2*C[5]) * a[4]
            + ( 3*C[3]**3*C[4]
              + 4*C[3]*C[4]**2
              - 2*C[4]*C[5]
              - 3*C[3]*C[6]
              +   C[7]
              ) * a[2]
        )

    return a 

def find_freq_shift(A, a): 
    total_shift = np.zeros(np.shape(A))
    for k in range(2, len(a)): 
        total_shift += a[k] * A**k
    return np.abs(total_shift)