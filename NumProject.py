'''
Numerical solver for Schroedinger's equation.
'''
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from cycler import cycler
from scipy import integrate

rc('text', usetex=True)
beta = 64  # beta = 2 m a^2 V0 / h^2


def find_eigenenergies(x_axis, accuracy, n_max):
    del_x = x_axis * 2
    potential_arr = np.array([potential(x) for x in x_axis])
    min_energy = 1 / (beta * del_x**2) + potential_arr
    energy = np.min(min_energy)

    n = 0
    energies = []
    while n <= n_max:
        energy = find_eigenenergy(x_axis, accuracy, energy, n)
        energies.append(energy)
        n += 1

    print(energies)
    return energies


def find_eigenenergy(x_axis, accuracy, energy, state):
    order = math.floor(math.log10(energy))
    energy = round(energy, -order)
    del_E = 10**order
    factor = 1.0

    while math.floor(math.log10(del_E)) > -accuracy:
        wave = shooting(energy, x_axis)
        print('{:.20f}, nodes: {:.1f}'.format(energy, nodes(wave)))
        if nodes(wave) > state:
            newfactor = -1.0
        else:
            newfactor = 1.0

        if newfactor != factor:
            del_E = 0.5 * del_E
        factor = newfactor
        energy = energy + factor*del_E
        # plt.figure()
        # plt.plot(x_axis, wave)
        # plt.show()

    return energy

def potential(x):
    if x > 0.5 or x < -0.5:
        return 1

    return 0

def nodes(x):
    neighbor = np.roll(x, 1)
    neighbor[0] = 0
    return len(np.where(x * neighbor < 0)[0])

def solver(wave, x, energy):
    dpdx = [wave[1], -beta * (energy - potential(x)) * wave[0]]
    return dpdx

def shooting(energy, x_axis):
    initial = [0.01, 0.001]
    solution = integrate.odeint(solver, initial, x_axis, args=(energy,))
    return solution[:, 0]

def areaunder(wave, x_axis, xmax=0):
    '''
    Function that returns the area under the graph of a mathematical
    function passed as an argument. The area is calculated from x=0
    to x=xmax
    '''
    if xmax == 0:
        xmax = x_axis[-1]

    # Calculating the probability distribution of the wave
    y_val = probability(wave)
    area = 0
    for i in range(0, len(y_val) - 1):
        if x_axis[i] < xmax:
            # The middle value between to y data points is used
            area += (y_val[i] + (y_val[i + 1] - y_val[i]) / 2) * \
                (x_axis[i + 1] - x_axis[i])
    area = area * 2
    print("Area " + str(area))
    return area


def probability(wave):
    '''
    Function returning probability distribution for a given real wave function
    '''
    return wave * wave


def normalization(wave, x_axis):
    '''
    Function that returns normalized wave function by calculating overall probability
    and dividing the wave function by this normalization constant
    '''
    constant = 1 / math.sqrt(areaunder(wave, x_axis))
    # List comprehension is used to divide all data points by the
    # constant = 1/K from Equation 5
    return wave * constant


def ploteigenfcs(x_axis, lines, title, file_name='', save=False):
    '''
    Plots three different functions onto one figure. Currently the color used are
    black, blue, red, magneto, cyan, and yellow
    '''
    # Get line to plot for potential energy V(x)
    energy_line = np.array([potential(x) for x in x_axis])

    plt.figure()
    plt.rc('lines', linewidth=1.5)
    plt.rc('axes', prop_cycle=cycler(
        'color', ['k', 'b', 'r', 'g', 'm', 'c', 'y']))
    lines.update({'Potential': energy_line})
    for desc, wave in lines.items():
        plt.plot(x_axis, wave, label=desc)

    ymax = 0
    ymin = 0
    for wave in lines.values():
        ymax = np.maximum(np.amax(wave, axis=0), ymax)
        ymin = np.minimum(np.amin(wave, axis=0), ymin)
    plt.ylim(ymin - 0.1, ymax + 0.1)
    plt.xlim(-3, 3)
    plt.xlabel(r'Radius $u=r/a_{0}$')
    plt.ylabel(r'Eigenfunction $\psi(u)$')
    plt.title(title, ha='center', fontsize=14)
    plt.legend()
    plt.grid()

    if save:
        plt.savefig('files/' + str(file_name) + '.eps', format='eps',
                    dpi=1000)  # Optional line: saves the plot as .eps file


def numerical_slv():
    '''
    Executes the numerical solver script and generates plots
    '''
    x_axis = np.linspace(-3.0, 3.0, 12000)

    # Find first 3 bound eigenenergies
    energies = find_eigenenergies(x_axis, 18, 0) # I found that the precision needs to be at least 11

    eigenfc = shooting(energies[0], x_axis)
    normalized = normalization(eigenfc, x_axis) + energies[0]
#    eigenfc_exc = shooting(energies[1], x_axis)
#    normalized_exc = normalization(eigenfc_exc, x_axis) + energies[1]
#    eigenfc_exc2 = shooting(energies[2], x_axis)
#    normalized_exc2 = normalization(eigenfc_exc2, x_axis) + energies[2]
    lines = {
        'Ground State': normalized
#        'First Excited': normalized_exc,
#        'Second Excited': normalized_exc2
    }

    ploteigenfcs(x_axis, lines,
                 r'\textbf{Eigenfunctions vs Radial Distance}')
    plt.show()

def test():
    x_axis = np.linspace(-3.0, 3.0, 12000)
    energies = find_eigenenergies(x_axis, 11, 0)

    eigenfc = shooting(energies[0], x_axis)

    plt.figure()
    plt.plot(x_axis, eigenfc)
    plt.show()

if __name__ == "__main__":
    test()
