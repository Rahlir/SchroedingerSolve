import matplotlib.pyplot as plt
from matplotlib import rc
import math
# The code here is almost identical to the other file 'NumProject.py'. Only some minor changes related to the new form
# of TISE can be noticed

rc('text', usetex=True)

# Eigenvalues for energies
energy1 = 1.0186929792
energy2 = 2.33800744
energy13 = 10.04012447596
plotsno = 0  # Global variable to keep track of number of figures plotted in one session


'''
First function to be called when generating a new eigenfunction. Arguments are
used to specify initial conditions, energy epsilon, and number of points 
(spacing)
'''
def shooting(initwave, initderiv, energy, noofpoints):
    # Initializing 3d array to store wave, derivatives, and position
    wave = [[initwave], [initderiv], [0]]
    # Calculate delta x based on number of points requested
    delta = 20/noofpoints
    return generate(wave, energy, delta)


'''
Function that is actually performing all the calculations for the new 
eigenfunction
'''
def generate(wave, energy, delta):
    while wave[2][-1] <= 20:
        # Equation 2 in code form
        newwave = wave[0][-1] + delta * wave[1][-1]
        # Equation 1 in code form
        newder = wave[1][-1] + delta * ((wave[2][-1]) - energy) * wave[0][-1]
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
    constant = 1/math.sqrt(areaunder(probabilitydstr, wave[2], 20))
    # List comprehension is used to divide all data points by the
    # constant = 1/K from Equation 5
    return [x*constant for x in wave[0]]


def ploteigenfc(y, ynorm, x, title, id):
    global plotsno
    plotsno += 1
    plt.figure(plotsno)
    plt.plot(x, y, color="black", linewidth=1.5, linestyle="-", label="Not Normalized")
    plt.plot(x, ynorm, color="blue", linewidth=1.5, linestyle="-", label="Normalized")
    plt.xlim(0, 20)
    ymin = min(y)
    ymax = max(y)
    if ymax > y[0] and ymax == y[-1]:
        maxindex = round(len(x)/2.5)
        ymax = max(y[0:maxindex])
    if ymin == y[-1] and ymin < -1:
        ymin = -0.1
    plt.ylim(ymin-0.1, ymax+0.1)
    plt.xlabel(r'Displacement $u$')
    plt.ylabel(r'Eigenfunction $\psi(u)$')
    plt.title(title, ha='center', fontsize=14)
    plt.legend()
    plt.grid()
    plt.savefig(str(id) + '.eps', format='eps', dpi=1000)  # Optional line: saves the plot as .eps file


eigenfc = shooting(0, 1, energy13, 200000)
normalized = normalization(eigenfc)
areaunder(probability(normalized), eigenfc[2], 20)

ploteigenfc(eigenfc[0], normalized, eigenfc[2], r'\textbf{Constant Force Potential: Thirteenth Excited State}' '\n' 
                                                '(Energy $\epsilon = ' + str(energy13) + '$)', "constantforce"+str(energy13))
plt.show()