'''
Numerical solver for Schroedinger's equation.
'''
import math
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)

BETA = 64
ENERGY = 0.0980245145  # Ground state energy
ENERGY2 = 0.38272399  # First excited state energy
ENERGY3 = 0.807899  # Second excited state energy
ENERGY_UNB = 2
PTS = 30000
DELTA_U = 4 / PTS


def potential(position):
    '''
    A function used to return a potential. When u < 1.0, then potential is 0.
    When u = 1.0, then potential is 1/2. And when u > 1.0, then potential is 1
    '''
    if position < 1 / 2:
        return 0
    elif position == 1 / 2:
        return 1 / 2
    return 1


def shooting(initwave, initderiv, energy, noofpoints):
    '''
    First function to be called when generating a new eigenfunction. Arguments are
    used to specify initial conditions, energy epsilon, and number of points
    (spacing)
    '''
    # Initializing 3d array to store wave, derivatives, and position
    wave = [[initwave], [initderiv], [0]]
    # Calculate delta x based on number of points requested
    delta = 2.5 / noofpoints
    return generate(wave, energy, delta)


def generate(wave, energy, delta):
    '''
    Function that is actually performing all the calculations for the new
    eigenfunction
    '''
    while wave[2][-1] <= 2.5:
        # Equation 2 in code form
        newwave = wave[0][-1] + delta * wave[1][-1]
        # Equation 1 in code form
        newder = wave[1][-1] - delta * BETA * \
            (energy - potential(wave[2][-1])) * wave[0][-1]
        wave[0].append(newwave)
        wave[1].append(newder)
        # Calculating new u and adding it to the position array
        wave[2].append(wave[2][-1] + delta)
    return wave


def areaunder(y_val, x_val, xmax):
    '''
    Function that returns the area under the graph of a mathematical
    function passed as an argument. The area is calculated from x=0
    to x=xmax
    '''
    area = 0
    for i in range(0, len(y_val) - 1):
        if x_val[i] < xmax:
            # The middle value between to y data points is used
            area += (y_val[i] + (y_val[i + 1] - y_val[i]) / 2) * \
                (x_val[i + 1] - x_val[i])
    area = area * 2
    print("Area " + str(area))
    return area


def probability(wave):
    '''
    Function returning probability distribution for a given real wave function
    '''
    return [x**2 for x in wave]


def normalization(wave):
    '''
    Function that returns normalized wave function by calculating overall probability
    and dividing the wave function by this normalization constant
    '''
    # Calculating the probability distribution of the wave
    probabilitydstr = probability(wave[0])
    constant = 1 / math.sqrt(areaunder(probabilitydstr, wave[2], 4))
    # List comprehension is used to divide all data points by the
    # constant = 1/K from Equation 5
    return [x * constant for x in wave[0]]


def ploteigenfc(y_val, x_val, title, file_name, plots_no):
    '''
    Function that plots the specified eigenfunction and saves eps file in files directory.
    Specify title included in the plots and file name
    '''
    plt.figure(plots_no)
    plt.plot(x_val, y_val, color="black", linewidth=1.5, linestyle="-")
    plt.xlim(0, 4)
    ymin = min(y_val)
    ymax = max(y_val)
    if ymax > y_val[0] and ymax == y_val[-1]:
        maxindex = round(len(x_val) / 2.5)
        ymax = max(y_val[0:maxindex])
    if ymin == y_val[-1] and ymin < -1:
        ymin = -1
    plt.ylim(ymin - 0.1, ymax + 0.1)
    plt.xlabel(r'Radius $u=r/a_{0}$')
    plt.ylabel(r'Eigenfunction $\psi(u)$')
    plt.title(title, ha='center', fontsize=14)
    plt.legend()
    plt.grid()
    plt.savefig('files/' + str(file_name) + '.eps',
                format='eps', dpi=1000)  # Optional line: saves the plot as .eps file


def ploteigenfcs(y1, y2, y3, x, title, file_name, plots_no): # TODO: Plot arbitrary no. of functions
    '''
    Plots three different functions onto one figure. Currently the color used are
    black, blue, and red
    '''
    plt.figure(plots_no)
    plt.plot(x, y1, color="black", linewidth=1.5,
             linestyle="-", label="Ground State")
    plt.plot(x, y2, color="blue", linewidth=1.5,
             linestyle="-", label="First Excited State")
    plt.plot(x, y3, color="red", linewidth=1.5,
             linestyle="-", label="Second Excited State")
    ymax = 0
    ymin = 0
    for y_val in [y1, y2, y3]:
        ymax = max(max(y_val), ymax)
        ymin = min(min(y_val), ymin)
    plt.ylim(ymin, ymax)
    plt.xlim(0, 250)
    plt.xlabel(r'Radius $r$ [nm]')
    plt.ylabel(r'Eigenfunction $\psi(u)$')
    plt.title(title, ha='center', fontsize=14)
    plt.legend()
    plt.grid()
    plt.savefig('files/' + str(file_name) + '.eps', format='eps', dpi=1000)
    return plots_no + 1


def numerical_slv():
    '''
    Executes the numerical solver script and generates plots
    '''
    plots_no = 0

    eigenfc = shooting(1.0, 0, ENERGY, 3000)
    normalized = normalization(eigenfc)
    eigenfc2 = shooting(0, 1.0, ENERGY2, 3000)
    normalized2 = normalization(eigenfc2)
    eigenfc3 = shooting(1.0, 0, ENERGY3, 3000)
    normalized3 = normalization(eigenfc3)

    plots_no = ploteigenfcs(normalized, normalized2, normalized3, eigenfc[2],
                            r'\textbf{Eigenfunctions vs Radial Distance}', 'tad', plots_no)
    plt.show()


if __name__ == "__main__":
    numerical_slv()
