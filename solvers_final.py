"""Python file containing the solvers for project 2."""

# import relevant libraries
import numpy as np
import matplotlib.pyplot as plt
import math

# import the parameters from the dictionaries file
from dictionaries import *

# Define the analytical solvers for Task C
# Define the simple functions first


# for the a coefficient
def a_analytic(epsilon):
    """Computes the a coefficient for calculation of the analytic solution of the ocean gyre simulation from Mushgrave (1985).
    
    Inputs:
    epsilon - constant
    
    Outputs:
    value of a coefficient.
    """ 
    return (-1 - np.sqrt(1 + (2*np.pi*epsilon)**2)) / (2*epsilon)

# for the b coefficient
def b_analytic(epsilon):
    """Computes the b coefficient for calculation of the analytic solution of the ocean gyre simulation from Mushgrave (1985).    
    
    Inputs:
    epsilon - constant
    
    Outputs:
    value of b coefficient.
    """
    return (-1 + np.sqrt(1 + (2*np.pi*epsilon)**2)) / (2*epsilon)


# for epsilon
def epsilon_analytic(gamma, L, beta):
    """Computes the value of episilon for the calculation of the analytic solution for the ocean gyre simulation from Mushgrave (1985).
    
    Inputs:
    gamma - linear drag coefficient (s^-1)
    L - dimensions of computational domain (m)
    beta - constant used for B-plane approximation (m^-1 s^-1)
    
    Outputs:
    value of epsilon.
    """
    return gamma / (L*beta)

# now define the functions f_1 and f_2 which use the simple functions above

# for the f_1 function
def f1(x, a, b):
    """f_1 function for calculating the analytical solution of the ocean gyre using methods from Mushgrave (1985).
    
    Inputs:
    x - the value of the x domain
    a - the a coefficient
    b - the b coefficient
    
    Outputs:
    value of f1 (at value of x).
    """
    
    # compute the numerator
    numerator = (math.exp(a) - 1) * b * math.exp(b * x) + (1 - math.exp(b)) * math.exp(a * x)
    
    # compute the denominator
    denominator = math.exp(b) - math.exp(a)
    
    return np.pi * (1 + numerator / denominator)

# for the f_2 function   
def f2(x, a, b):
    """f_2 function for calculating the analytical solution of the ocean gyre using methods from Mushgrave (1985).

    Inputs:
    x - the value of the x domain
    a - the a coefficient
    b - the b coefficient
    
    Outputs:
    value of f2 (at value of x).
    """
    
    # compute the numerator
    numerator = (math.exp(a) - 1)*b*math.exp(b*x) + (1 - math.exp(b))*a*math.exp(a*x)
    
    # compute the denominator
    denominator = math.exp(b) - math.exp(a)
    
    return numerator/denominator


# define the function for the analytic solution
def analytic_solution(params_analytic):
    """Analytic solver for the SWEs using equations (3), (4) and (5) from project brief specifying
    the solutions at (x, y) using methods from Mushgrave (1985).

    Inputs:
    params_analytic - dictionary containing the constants to be used

    Outputs:
    x - array of x values
    y - array of y values
    u - analytic solution for fluid motion in x-direction
    v - analytic solution for fluid motion in y-direction
    eta - analytic solution for the deviation of water surface from its initial level.
    """

     # extract the constants
    f0, beta, g, gamma, rho, H, tau0, L, eta0 = [params_analytic[key] for key in
                                                 ['f0', 'beta', 'g', 'gamma', 'rho', 'H', 'tau0', 'L', 'eta0']]

    # extract the grid characteristics
    x_points, y_points = params_analytic['x_points'], params_analytic['y_points']
    gridbox_size = params_analytic['gridbox_size']

    # set up the x and y domains
    x = np.linspace(0, x_points*gridbox_size, x_points+1)
    y = np.linspace(0, y_points*gridbox_size, y_points+1)

    # set up the arrays for the solution
    u, v, eta = np.zeros((x_points, y_points)), np.zeros((x_points, y_points)), np.zeros((x_points, y_points))

    # compute the constants
    epsilon = epsilon_analytic(gamma, L, beta)
    a, b = a_analytic(epsilon), b_analytic(epsilon)
    tau_coeff = tau0 / (np.pi*gamma*rho*H)

    # simplify the functions
    pi, sin, cos = np.pi, np.sin, np.cos

    # compute the solutions
    # note the indexing - numpy arrays are indexed as [y, x]
    for j in range(y_points):
        for i in range(x_points):
            
            # compute the solutions
            u[j, i] = -tau_coeff * f1(x[i]/L, a, b) * cos(pi*y[j]/L)
            
            v[j, i] = tau_coeff * f2(x[i]/L, a, b) * sin(pi*y[j]/L)
            
            eta[j, i] = eta0 + tau_coeff * (f0*L/g) * (
                gamma/(f0*pi) * f2(x[i]/L, a, b) * cos(pi*y[j]/L)
                + 1/pi * f1(x[i]/L, a, b) * (
                    sin(pi*y[j]/L) * (1 + beta*y[j]/f0)
                    + beta*L/(f0*pi) * cos(pi*y[j]/L)
                )
            )

    return u, v, eta, x, y


# now we set up the functions for the numerical solution
# -----------------------------------------------------#

# we need to compute coriolis and map this onto both the u-grid and the v-grid
def coriolis(y_ugrid, y_vgrid, x_points, f0, beta):
    """Function for computing the coriolis parameter at the u-grid and v-grid points using distance from origin.
    
    Inputs:
    y_ugrid - the y-coordinates of the u-grid points
    y_vgrid - the y-coordinates of the v-grid points
    x_points - the number of x-grid points
    f0 - the reference coriolis parameter
    beta - the gradient of the coriolis parameter
    
    Outputs:
    f_u - the coriolis parameter at the u-grid points
    f_v - the coriolis parameter at the v-grid points
    """

    # calcuate the coriolis parameter at the u-grid points mapped onto y
    f_u = f0 + (beta*y_ugrid)
    
    # calculate the coriolis parameter at the v-grid points mapped onto y
    f_v = f0 + (beta*y_vgrid)

    # repeat each element of f_u and f_v onto the u-grid and v-grid respectively
    f_u = np.tile(f_u, (x_points + 1, 1)).T
    f_v = np.tile(f_v, (x_points, 1)).T

    # return the coriolis parameter at the u-grid and v-grid points
    return f_u, f_v


# # implerment the function for computing the zonal wind stress at the u-grid points
def zonal_wind_stress(y_ugrid, x_points, L, tau0):
    """Function for computing the zonal wind stress at the u-grid points using distance from origin.
    
    Inputs:
    y_ugrid - the y-coordinates of the u-grid points
    x_points - the number of x-grid points
    L - the length of the domain
    tau0 - the maximum wind stress
    
    Outputs:
    tau - the zonal wind stress at the u-grid points
    """

    # repeat the y-coordinates of u-grid points along the x-axis
    # then reshape the array to be the same shape as the u-grid
    y = np.repeat(y_ugrid, x_points+1).reshape(len(y_ugrid), x_points+1)

    # compute the zonal wind stress at the u-grid points using the tau formula
    tau = tau0 * -np.cos(np.pi * y / L)

    # return the zonal wind stress at the u-grid points
    return tau

# we need to calculate the gradient of eta
# mapped onto the x and y points of our arakawa c-grid

# define a function for computing the gradient of eta mapped onto the x and y points of our arakawa c-grid
def eta_gradient(eta, x_points, y_points, dx, dy):
    """Function for computing the gradient of eta on the arakawa c-grid for both the u grid and the v grid.
    
    Inputs:
    eta - the free surface height
    x_points - the number of x-grid points
    y_points - the number of y-grid points
    dx - the grid spacing in the x-direction
    dy - the grid spacing in the y-direction
    
    Outputs:
    deta_dx - the gradient of eta in the x-direction
    deta_dy - the gradient of eta in the y-direction
    """
    # create arrays of zeros for the gradient of eta in the x-direction and y-direction
    array_zeros_x = np.zeros((y_points, 1))
    array_zeros_y = np.zeros((1, x_points))

    # concatenate the zeros to the start/end of the eta array for the x-direction
    eta_zeros_end_x = np.concatenate((eta, array_zeros_x), axis=1)
    eta_zeros_start_x = np.concatenate((array_zeros_x, eta), axis=1)

    # concatenate the zeros to the start/end of the eta array for the y-direction
    eta_zeros_end_y = np.concatenate((eta, array_zeros_y), axis=0)
    eta_zeros_start_y = np.concatenate((array_zeros_y, eta), axis=0)

    # calculate the gradient of eta in the x-direction
    deta_dx = (eta_zeros_end_x - eta_zeros_start_x)/(dx)
    
    # calculate the gradient of eta in the y-direction
    deta_dy = (eta_zeros_end_y - eta_zeros_start_y)/(dy)

    # return the gradient of eta in the x-direction and y-direction
    return deta_dx, deta_dy


# For the forward-backward time scheme, we must set up a function(s) which calculates \\
# the average of the four adjacent V-grid points for each U-grid point (and vice versa).
# define a function for mapping the v-grid points onto the u-grid points
def v_to_ugrid_mapping(v, y_points):
    """Function for mapping the v-grid values onto the u-grid.
    
    Inputs:
    v - the v-grid values
    y_points - the number of y-grid points
    
    Outputs:
    index of v mapped onto the u-grid set up for forward-backward time scheme
    """
    # create an array for appending zeros to the v-grid values
    # at the boundaries of the domain
    v_boundary = np.zeros((y_points, 1))

    # specify the index of v_j_i and v_j_iminus1
    v_j_i, v_jminus1_i = v[1:, :], v[:-1, :]

    # append zeros to the v-grid values at the boundaries of the domain
    # and concatenate the arrays
    v_same = np.concatenate((v_jminus1_i, v_boundary), axis=1)
    v_shift_up_row = np.concatenate((v_boundary, v_jminus1_i), axis=1)
    v_shift_down_row = np.concatenate((v_j_i, v_boundary), axis=1)
    v_shift_down_row_right_column = np.concatenate((v_boundary, v_j_i), axis=1)

    # take the average of the values mapped onto the u-grid points
    v_to_ugrid = (v_same + v_shift_up_row + v_shift_down_row + v_shift_down_row_right_column)/4

    # return the index of v mapped onto the u-grid set up for forward-backward time scheme
    return v_to_ugrid

# define a function for mapping the u-grid points onto the v-grid points
def u_to_vgrid_mapping(u, x_points):
    """Function for mapping the u-grid values onto the v-grid.
    
    Inputs:
    u - the u-grid values
    x_points - the number of x-grid points
    
    Outputs:
    index of u mapped onto the v-grid set up for forward-backward time scheme
    """

    # create an array for appending zeros to the u-grid values
    # at the boundaries of the domain
    u_boundary = np.zeros((1, x_points))

    # specify the index of u_j_i and u_j_iminus1
    u_j_i, u_j_iminus1 = u[:, 1:], u[:, :-1]

    # specify the index of u mapped onto the v-grid points
    # set up for the forward-backward time scheme
    u_same = np.concatenate((u_j_iminus1, u_boundary), axis=0)
    u_shift_right_column = np.concatenate((u_boundary, u_j_iminus1), axis=0)
    u_shift_left_column = np.concatenate((u_j_i, u_boundary), axis=0)
    u_shift_left_column_up_row = np.concatenate((u_boundary, u_j_i), axis=0)

    # take the average of the values mapped onto the v-grid points
    u_to_vgrid = (u_same + u_shift_right_column + u_shift_left_column + u_shift_left_column_up_row)/4

    # return the index of u mapped onto the v-grid set up for forward-backward time scheme
    return u_to_vgrid


# for Task E we calculate the energy of the analytic and numerical solutions

# we want to compute the energy to test the stability of the model
def energy(u, v, eta, dx, dy, rho, H, g):
    """Function for computing the energy of the model.
    
    Inputs:
    u - the zonal velocity
    v - the meridional velocity
    eta - the free surface height
    dx - the grid spacing in the x-direction
    dy - the grid spacing in the y-direction
    rho - the density of the fluid
    H - the depth of the fluid
    g - the gravitational acceleration
    
    Outputs:
    energy - the total energy of the model over the domain.
    """
    # calculate the kinetic energy
    kinetic_energy = H*(u**2 + v**2)
    
    # calculate the potential energy
    potential_energy = g*eta**2
    
    # calculate the total energy per grid point
    energy = 1/2*rho*(kinetic_energy + potential_energy)
    
    # calculate the total energy over the whole domain
    total_energy = np.sum(energy)*dx*dy

    # return the total energy
    return total_energy


# define a function which computes the energy in the FB scheme
def compute_energy(energy_array, energy_difference, energy_analytic, u, v, eta, u_analytic, v_analytic, eta_analytic, dx, dy, rho, H, g):
    """Function for computing the energy of the model.
    
    Inputs:
    energy_array - the array of energy values
    energy_difference - the array of energy difference values
    energy_analytic - the array of analytic energy values
    u - the zonal velocity
    v - the meridional velocity
    eta - the free surface height
    u_analytic - the analytic zonal velocity
    v_analytic - the analytic meridional velocity
    eta_analytic - the analytic free surface height
    dx - the grid spacing in the x-direction
    dy - the grid spacing in the y-direction
    rho - the density of the fluid
    H - the depth of the fluid
    g - the gravitational acceleration
    
    Outputs:
    energy_array - the array of energy values
    energy_difference - the array of energy difference values
    energy_analytic - the array of analytic energy values
    """
    # interpolate u and v onto the eta grid
    u_eta = (u[:, :-1] + u[:, 1:])/2
    v_eta = (v[:-1, :] + v[1:, :])/2

    # calculate the difference between the analytic and numerical solutions
    u_diff = u_analytic - u_eta
    v_diff = v_analytic - v_eta
    eta_diff = eta_analytic - eta

    # calculate the energy at the current time step
    # and add it to the array of energy values
    energy_array = np.append(energy_array, energy(u_eta, v_eta, eta, dx, dy, rho, H, g))
    energy_difference = np.append(energy_difference, energy(u_diff, v_diff, eta_diff, dx, dy, rho, H, g))

    # calculate the analytic energy at the current time step
    # create an array of the correct length for the analytic energy
    energy_analytic = np.append(energy_analytic, energy(u_analytic, v_analytic, eta_analytic, dx, dy, rho, H, g))

    # return the energy arrays
    return energy_array, energy_difference, energy_analytic, u_diff, v_diff, eta_diff

# now implement the forward-backward time scheme
# first u before v
# then v before u
# ahead of this we set up some simple functions for the FB function

# define a function which resets the boundary conditions for u
def reset_u_boundary_conditions(u, x_points):
    """Function for resetting the boundary conditions for u."""

    # set the boundary conditions for u
    u[:, 0] = 0
    u[:, x_points] = 0

    # return u
    return u

# define a function which resets the boundary conditions for v
def reset_v_boundary_conditions(v, y_points):
    """Function for resetting the boundary conditions for v."""

    # set the boundary conditions for v
    v[0, :] = 0
    v[y_points, :] = 0

    # return v
    return v

# define a function which computes u before v in the FB scheme
def compute_u_before_v(u, v, eta, x_points, y_points, dx, dy, dt, rho, H, g, gamma, tau_zonal, tau_meridional, coriolis_u, coriolis_v, v_to_ugrid_mapping, u_to_vgrid_mapping):
    """Function for computing u before v in the forward-backward time scheme."""

    # index the u and v arrays to the correct grid points
    u_j_i, u_j_iminus1 = u[:, 1:], u[:, :-1]
    v_j_i, v_jminus1_i = v[1:, :], v[:-1, :]

    # compute eta at the next time step
    eta_next = eta - H*dt*((u_j_i - u_j_iminus1)/dx + (v_j_i - v_jminus1_i)/dy)

    # compute the gradient of eta at the next time step
    deta_dx, deta_dy = eta_gradient(eta_next, x_points, y_points, dx, dy)

    # compute the zonal velocity at the next time step
    u_next = u + coriolis_u*dt*v_to_ugrid_mapping(v, y_points) - g*dt*deta_dx - gamma*dt*u + (tau_zonal/(rho*H))*dt

    # reset the boundary conditions for u
    u_next = reset_u_boundary_conditions(u_next, x_points)

    # compute the meridional velocity at the next time step
    # notice the use of u_next here
    v_next = v - coriolis_v*dt*u_to_vgrid_mapping(u_next, x_points) - g*dt*deta_dy - gamma*dt*v + (tau_meridional/(rho*H))*dt

    # reset the boundary conditions for v
    v_next = reset_v_boundary_conditions(v_next, y_points)

    # return the next time step for eta, u, and v
    return eta_next, u_next, v_next

# define a function which computes v before u in the FB scheme
def compute_v_before_u(u, v, eta, x_points, y_points, dx, dy, dt, rho, H, g, gamma, tau_zonal, tau_meridional, coriolis_u, coriolis_v, v_to_ugrid_mapping, u_to_vgrid_mapping):
    """Function for computing v before u in the forward-backward time scheme."""

    # index the u and v arrays to the correct grid points
    u_j_i, u_j_iminus1 = u[:, 1:], u[:, :-1]
    v_j_i, v_jminus1_i = v[1:, :], v[:-1, :]

    # compute eta at the next time step
    eta_next = eta - H*dt*((u_j_i - u_j_iminus1)/dx + (v_j_i - v_jminus1_i)/dy)

    # compute the gradient of eta at the next time step
    deta_dx, deta_dy = eta_gradient(eta_next, x_points, y_points, dx, dy)

    # compute the meridional velocity at the next time step
    v_next = v - coriolis_v*dt*u_to_vgrid_mapping(u, x_points) - g*dt*deta_dy - gamma*dt*v + (tau_meridional/(rho*H))*dt

    # reset the boundary conditions for v
    v_next = reset_v_boundary_conditions(v_next, y_points)

    # compute the zonal velocity at the next time step
    # notice the use of v_next here
    u_next = u + coriolis_u*dt*v_to_ugrid_mapping(v_next, y_points) - g*dt*deta_dx - gamma*dt*u + (tau_zonal/(rho*H))*dt

    # reset the boundary conditions for u
    u_next = reset_u_boundary_conditions(u_next, x_points)

    # return the next time step for eta, u, and v
    return eta_next, u_next, v_next

# now we implement the forward-backward time scheme
# ------------------------------------------------#

# implement the forward-backward time scheme
def forward_backward_time_scheme(params):
    """Function for solving the shallow water equations using the forward-backward time scheme.
    
    Inputs:
    params - a dictionary containing the parameters for the model
    
    Outputs:
    u - the zonal velocity
    v - the meridional velocity
    eta - the free surface height
    """
    # extract the constants from the dictionary
    f0, g, rho, H, tau0, beta, gamma = params['f0'], params['g'], params['rho'], params['H'], params['tau0'], params['beta'], params['gamma']
    
    # extract the spatial and temporal parameters from the dictionary
    L, dx, dy, x_points, y_points, dt, nt = params['L'], params['dx'], params['dy'], params['x_points'], params['y_points'], params['dt'], params['nt']
    
    # resolution flags
    use_higher_resolution, use_highest_resolution, task = params['use_higher_resolution'], params['use_highest_resolution'], params['task']

    # define the x and y arrays
    x_plotting, y_plotting, x_plotting_eta, y_u, y_v = np.arange(x_points + 1)*dx, np.arange(y_points + 1)*dy, np.arange(x_points)*dx, np.arange(y_points)*dy, np.arange(y_points + 1)*dy

    # define the 2D arrays for the zonal velocity, meridional velocity and free surface height
    u, v, eta = np.zeros((y_points, x_points + 1)), np.zeros((y_points + 1, x_points)), np.zeros((y_points, x_points))

    # define the zonal and meridional wind stress
    tau_zonal, tau_meridional = zonal_wind_stress(y_u, x_points, L, tau0), params['tau_meridional']

    # compute the values for the coriolis parameter
    coriolis_u, coriolis_v = coriolis(y_u, y_v, x_points, f0, beta)

    # determine which analytical solution to use
    if use_higher_resolution == 'True' or use_highest_resolution == 'True':
        u_analytic, v_analytic, eta_analytic, x, y = analytic_solution(params_analytic_higher_res if use_higher_resolution == 'True' else params_analytic_highest_res)
    else:
        u_analytic, v_analytic, eta_analytic, x, y = analytic_solution(params_analytic)

    # set up empty arrays for energy
    energy_array, energy_difference, energy_analytic = np.zeros(0), np.zeros(0), np.zeros(0)
    # ditto for time
    time_array = np.empty(0)

    # compute the analytic energy
    #energy_analytic = energy(u_analytic, v_analytic, eta_analytic, dx, dy, rho, H, g)
    #print('the value of energy analytic outside the loop is: ' + str(energy_analytic))

    # loop over time with intervals of 2 up to nt-2 for the forward-backward time scheme
    for i in range(0, nt - 2, 2):

        # compute the energy using the function
        energy_array, energy_difference, energy_analytic, u_diff, v_diff, eta_diff = compute_energy(energy_array, energy_difference, energy_analytic, u, v, eta, u_analytic, v_analytic, eta_analytic, dx, dy, rho, H, g)
        
        # set up the time array
        time_array = np.append(time_array, [i*dt])
        
        # compute the energy by calling the energy function
        # energy_array, energy_difference, energy_analytic, u_diff, v_diff, eta_diff = compute_energy(energy_array, energy_difference, energy_analytic, u, v, eta, u_analytic, v_analytic, eta_analytic, dx, dy, rho, H, g)

        # first compute u before v
        eta_next, u_next, v_next = compute_u_before_v(u, v, eta, x_points, y_points, dx, dy, dt, rho, H, g, gamma, tau_zonal, tau_meridional, coriolis_u, coriolis_v, v_to_ugrid_mapping, u_to_vgrid_mapping)

        # then compute v before u
        eta_next2, u_next2, v_next2 = compute_v_before_u(u_next, v_next, eta_next, x_points, y_points, dx, dy, dt, rho, H, g, gamma, tau_zonal, tau_meridional, coriolis_u, coriolis_v, v_to_ugrid_mapping, u_to_vgrid_mapping)

        # update the u, v and eta arrays with the value at the next next time step
        u = u_next2.copy()
        v = v_next2.copy()
        eta = eta_next2.copy()

    # compute the value of eta0 for use in the analytic solution above
    print('eta0 = ', eta[int(y_points/2), 0])

    # create the plots for each task
    if task == 'D1' or task == 'D2':
        # use the plotting function to plot the results
        plot_variables_D1_D2(u, v, eta, x_plotting, y_plotting, x_plotting_eta, y_points, params)

    elif task == 'D3':
        # plot the results as three seperate plots, with a colour bar for each
        # create a 2D contour plot of eta
        plot_differences_D3(u_diff, v_diff, eta_diff, x_plotting, y_plotting, x_plotting_eta, y_points, params)

    elif task == 'E':
        # plot the energy as a function of time
        plot_energy(energy_analytic, energy_array, energy_difference, time_array, params)
        # we want to return the energy difference array
        return energy_difference, energy_analytic, time_array
    

# now we define the plotting functions
# -----------------------------------#

# now set up a function to plot the results of the analytic solution (for Task C)
def plotting_taskC(params_analytic):
    """Function for plotting the results of the analytic solution for the ocean gyre simulation.

    Inputs:
    params_analytic - the dictionary containing the constants to be used (in this case 'params_analytic')

    Outputs:
    None
    """

    # compute the analytic solution
    u, v, eta, x, y = analytic_solution(params_analytic)

    # set up the subplots with a shared x and y axis
    fig, axs = plt.subplots(ncols=3, sharex=True, sharey=True, figsize=(15,5))

    # set the titles of the subplots
    titles = ['eta analytic', 'u analytic', 'v analytic']

    # iterate over the subplots and plot the results
    for ax, title, data, fig_name in zip(axs, titles, [eta, u, v], ['eta_fig_name', 'u_fig_name', 'v_fig_name']):
        
        ax.set_title(title)
        ax.set_xlabel('x (km)')
        ax.set_ylabel('y (km)')
        # use pcolormesh to plot the data
        ax.pcolormesh(x/1000, y/1000, data, cmap='jet')
        ax.set_aspect('equal')
        # save the figure using the appropriate file name from the dictionary
        plt.savefig(params_analytic[fig_name] + '.png', dpi=300)

    # show the plots
    plt.show()

# define a plotting function for Task D1 and D2
def plot_variables_D1_D2(u, v, eta, x_plotting, y_plotting, x_plotting_eta, y_points, params):
    """
    Plots u against x for the southern edge of the basin, v against y for the western edge of the basin,
    eta against x for the middle of the gyre, and a 2D contour plot of eta.

    Parameters:
    u (ndarray): Array of zonal velocity values.
    v (ndarray): Array of meridional velocity values.
    eta (ndarray): Array of surface displacement values.
    x_plotting (ndarray): Array of x values for plotting.
    y_plotting (ndarray): Array of y values for plotting.
    x_plotting_eta (ndarray): Array of x values for plotting eta.
    y_points (int): Number of points in the y direction.
    params (dict): Dictionary containing plot parameters.

    Returns:
    None
    """

    # plot u against x for for the southern edge of the basin
    fig1, ax1 = plt.subplots(figsize=(6, 6))
    ax1.plot(x_plotting/1000, u[0, :], label='numerical u')
    # set the x label
    ax1.set_xlabel('x (km)')
    # set the y label
    ax1.set_ylabel('u (m/s)')
    # set the title
    ax1.set_title('Zonal velocities at the southern edge of the basin')
    # save the plot
    fig1.savefig(params['u_fig_name'] + '.png')

    # plot v against y for for the western edge of the basin
    fig2, ax2 = plt.subplots(figsize=(6, 6))
    ax2.plot(y_plotting/1000, v[:, 1], label='numerical v')
    # set the x label
    ax2.set_xlabel('y (km)')
    # set the y label
    ax2.set_ylabel('v (m/s)')
    # set the title
    ax2.set_title('Meridional velocities at the western edge of the basin')
    # save the plot
    fig2.savefig(params['v_fig_name'] + '.png')

    # plot eta against x for for the middle of the gyre
    fig3, ax3 = plt.subplots(figsize=(6, 6))
    ax3.plot(x_plotting_eta/1000, eta[int(y_points/2), :], label='numerical eta')
    # set the x label
    ax3.set_xlabel('x (km)')
    # set the y label
    ax3.set_ylabel('eta (m)')
    # set the title
    ax3.set_title('Surface displacement at the middle of the gyre')
    # save the plot
    fig3.savefig(params['eta_fig_name'] + '.png')

    # create a 2D contour plot of eta
    fig4, ax4 = plt.subplots(figsize=(10, 8))
    # plot the contour
    contour = ax4.pcolormesh(x_plotting/1000, y_plotting/1000, eta[:,:], cmap='jet')
    # set the x label
    ax4.set_xlabel('x (km)')
    # set the y label
    ax4.set_ylabel('y (km)')
    # add a colour bar
    fig4.colorbar(contour)
    # set the title
    ax4.set_title('Surface displacement contour field')
    # save the plot
    fig4.savefig(params['eta_contour_fig_name'] + '.png')

    # specify the plots to show
    plt.show()


# create a function for the plotting in task D3
def plot_differences_D3(u_diff, v_diff, eta_diff, x_plotting, y_plotting, x_plotting_eta, y_points, params):
    """Plots the differences between the analytical and numerical solutions for u, v, and eta."""

    # plot the difference between the analytical and numerical solutions for u
    fig1, ax1 = plt.subplots(figsize=(10, 8))
    # plot the contour
    contour = ax1.pcolormesh(x_plotting/1000, y_plotting/1000, u_diff[:,:], cmap='jet')
    # set the x label
    ax1.set_xlabel('x (km)')
    # set the y label
    ax1.set_ylabel('y (km)')
    # add a colour bar
    fig1.colorbar(contour)
    # set the title
    ax1.set_title('Zonal velocity (u) difference')
    # save the plot
    fig1.savefig(params['u_fig_name'] + '.png')

    # plot the difference between the analytical and numerical solutions for v
    fig2, ax2 = plt.subplots(figsize=(10, 8))
    # plot the contour
    contour = ax2.pcolormesh(x_plotting/1000, y_plotting/1000, v_diff[:,:], cmap='jet')
    # set the x label
    ax2.set_xlabel('x (km)')
    # set the y label
    ax2.set_ylabel('y (km)')
    # add a colour bar
    fig2.colorbar(contour)
    # set the title
    ax2.set_title('Meridional velocity (v) difference')
    # save the plot
    fig2.savefig(params['v_fig_name'] + '.png')

    # plot the difference between the analytical and numerical solutions for eta
    fig3, ax3 = plt.subplots(figsize=(10, 8))
    # plot the contour
    contour = ax3.pcolormesh(x_plotting/1000, y_plotting/1000, eta_diff[:,:], cmap='jet')
    # set the x label
    ax3.set_xlabel('x (km)')
    # set the y label
    ax3.set_ylabel('y (km)')
    # add a colour bar
    fig3.colorbar(contour)
    # set the title
    ax3.set_title('Surface displacement (eta) difference')
    # save the plot
    fig3.savefig(params['eta_fig_name'] + '.png')

    # specify the plots to show
    plt.show()

# function for plotting energy
def plot_energy(energy_analytic, energy_array, energy_difference, time_array, params):
    """Plots the energy against time for the analytical and numerical solutions."""

    # print the values of the energy
    print('The energy for the analytical solution is: ' + str(energy_analytic))
    print('The energy for the steady state numerical solution is: ' + str(energy_array[-1]))

    # plot the energy against time for the analytical and numerical solutions
    fig1, ax1 = plt.subplots(figsize=(6, 6))
    ax1.plot(time_array/86400, energy_analytic, label='analytical energy', color='red')
    ax1.plot(time_array/86400, energy_array, label='numerical energy', color='blue')
    # set the x label
    ax1.set_xlabel('time (days)')
    # set the y label
    ax1.set_ylabel('energy (J)')
    # set the title
    ax1.set_title('Energy against time for numerical and analytical solutions')
    # include legend
    ax1.legend()
    # save the plot
    fig1.savefig(params['energy_fig_name'] + '.png')

    # plot the difference between the analytical and numerical solutions for energy against time
    fig2, ax2 = plt.subplots(figsize=(6, 6))
    ax2.plot(time_array/86400, energy_difference, label='energy difference', color='green')
    # set the x label
    ax2.set_xlabel('time (days)')
    # set the y label
    ax2.set_ylabel('energy difference (J)')
    # set the title
    ax2.set_title('Energy difference between the numerical and analytical solutions')
    # include legend
    ax2.legend()
    # save the plot
    fig2.savefig(params['energy_difference_fig_name'] + '.png')

    # specify the plots to show
    plt.show()

# function for plotting the energy differences for the different grid sizes
def plot_energy_difference(params_100, params_50, params_10):
    """Plots the energy differences for the different grid sizes on a single plot."""

    # call the scheme for the 100m grid size
    energy_difference_100, energy_analytic, time_array = forward_backward_time_scheme(params_100)

    # call the scheme for the 50m grid size
    energy_difference_50, energy_analytic_50, time_array_50 = forward_backward_time_scheme(params_50)

    # call the scheme for the 10m grid size
    energy_difference_10, energy_analytic_10, time_array_10 = forward_backward_time_scheme(params_10)

    # plot the energy difference against time for the different grid sizes
    fig1, ax1 = plt.subplots(figsize=(6, 6))
    ax1.plot(time_array/86400, energy_difference_100, label='100km grid size', color='red')
    ax1.plot(time_array_50/86400, energy_difference_50, label='50km grid size', color='blue')
    ax1.plot(time_array_10/86400, energy_difference_10, label='10km grid size', color='green')
    # set the x label
    ax1.set_xlabel('time (days)')
    # set the y label
    ax1.set_ylabel('energy difference (J)')
    # set the title
    ax1.set_title('Energy difference between the numerical and analytical solutions for different grid sizes')
    # include legend
    ax1.legend()
    # save the plot
    fig1.savefig('energy_difference_for_diff_grids.png')
    # specify the plots to show
    plt.show()