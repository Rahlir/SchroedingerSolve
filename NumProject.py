import matplotlib.pyplot as plt
from matplotlib import rc
import math

rc('text', usetex=True)

beta = 64
energy = 0.0980245145  # Ground state energy
energy2 = 0.38272399  # First excited state energy
energy3 = 0.807899  # Second excited state energy
unbounden = 2
pts = 30000
deltau = 4/pts
plotsno = 0  # Global variable to keep track of number of figures plotted in one session


'''
A function used to return a potential. When u < 1.0, then potential is 0. When u = 1.0, then potential is 1/2. And
when u > 1.0, then potential is 1
'''
def potential(position):
    if position < 1/2:
        return 0
    elif position == 1/2:
        return 1/2
    else:
        return 1

'''
First function to be called when generating a new eigenfunction. Arguments are
used to specify initial conditions, energy epsilon, and number of points 
(spacing)
'''
def shooting(initwave, initderiv, energy, noofpoints):
    # Initializing 3d array to store wave, derivatives, and position
    wave = [[initwave], [initderiv], [0]]
    # Calculate delta x based on number of points requested
    delta = 2.5/noofpoints
    return generate(wave, energy, delta)


'''
Function that is actually performing all the calculations for the new 
eigenfunction
'''
def generate(wave, energy, delta):
    while wave[2][-1] <= 2.5:
        # Equation 2 in code form
        newwave = wave[0][-1] + delta * wave[1][-1]
        # Equation 1 in code form
        newder = wave[1][-1] - delta * beta * \
                    (energy - potential(wave[2][-1])) * wave[0][-1]
        wave[0].append(newwave)
        wave[1].append(newder)
        # Calculating new u and adding it to the position array
        wave[2].append(wave[2][-1]+delta)
    return wave


'''
Function that returns the area under the graph of a mathematical 
function passed as an argument. The area is calculated from x=0
to x=xmax
'''
def areaunder(y, x, xmax):
    area = 0
    for i in range(0, len(y)-1):
        if x[i] < xmax:
            # The middle value between to y data points is used
            area += (y[i] + (y[i+1]-y[i])/2) * (x[i+1]-x[i])
    area = area*2
    print("Area "+str(area))
    return area


def probability(wave):
    return [x**2 for x in wave]


def normalization(wave):
    # Calculating the probability distribution of the wave
    probabilitydstr = probability(wave[0])
    constant = 1/math.sqrt(areaunder(probabilitydstr, wave[2], 4))
    # List comprehension is used to divide all data points by the
    # constant = 1/K from Equation 5
    return [x*constant for x in wave[0]]


def ploteigenfc(y, x, title, id):
    global plotsno
    plotsno += 1
    plt.figure(plotsno)
    plt.plot(x, y, color="black", linewidth=1.5, linestyle="-")
    plt.xlim(0, 4)
    ymin = min(y)
    ymax = max(y)
    if ymax > y[0] and ymax == y[-1]:
        maxindex = round(len(x)/2.5)
        ymax = max(y[0:maxindex])
    if ymin == y[-1] and ymin < -1:
        ymin = -1
    plt.ylim(ymin-0.1, ymax+0.1)
    plt.xlabel(r'Radius $u=r/a_{0}$')
    plt.ylabel(r'Eigenfunction $\psi(u)$')
    plt.title(title, ha='center', fontsize=14)
    plt.legend()
    plt.grid()
    plt.savefig('files/'+str(id) + '.eps', format='eps', dpi=1000)  # Optional line: saves the plot as .eps file


def ploteigenfcs(y1, y2, y3, x, title, id):
    plt.figure(1)
    plt.plot(x, y1, color="black", linewidth=1.5, linestyle="-", label="Ground State")
    plt.plot(x, y2, color="blue", linewidth=1.5, linestyle="-", label="First Excited State")
    plt.plot(x, y3, color="red", linewidth=1.5, linestyle="-", label="Second Excited State")
    ymax = 0
    ymin = 0
    for y in [y1, y2, y3]:
        ymax = max(max(y), ymax)
        ymin = min(min(y), ymin)
    plt.ylim(ymin, ymax)
    plt.xlim(0, 250)
    plt.xlabel(r'Radius $r$ [nm]')
    plt.ylabel(r'Eigenfunction $\psi(u)$')
    plt.title(title, ha='center', fontsize=14)
    plt.legend()
    plt.grid()
    plt.savefig('files/'+str(id) + '.eps', format='eps', dpi=1000)

eigenfc = shooting(1.0, 0, energy, 3000)
normalized = normalization(eigenfc)
areaunder(probability(normalized), eigenfc[2], 1.8)
eigenfc[2] = [x * 100 for x in eigenfc[2]]

eigenfc2 = shooting(0, 1.0, energy2, 3000)
normalized2 = normalization(eigenfc2)
areaunder(probability(normalized2), eigenfc2[2], 1.8)

eigenfc3 = shooting(1.0, 0, energy3, 3000)
normalized3 = normalization(eigenfc3)
areaunder(probability(normalized3), eigenfc3[2], 2.5)

ploteigenfcs(normalized, normalized2, normalized3, eigenfc[2], r'\textbf{Eigenfunctions vs Radial Distance}',
             'tad')

plt.show()