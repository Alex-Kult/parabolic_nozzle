from math import sin, cos, tan, atan, sqrt, radians, degrees
import nozzle_lib as noz
import numpy as np
import matplotlib.pyplot as plt
##########################

# Input Variables
rad_t = 8.405 #Throat radius [mm]
gamma = 1.16 #Ratio of specific heats
mach_e = 2.6541 #Exit Mach number
ratio_area = 5.1359 #Ratio of the exit area to the throat area

frac_con = 1 #Fraction of length of a 15 deg conical nozzle (Can go down to 70% without significant efficiency losses to save weight)

##########################
# Initial Expansion Angle
mu = noz.prandtl_meyer(gamma, mach_e)
ang_exp = 0.5*mu

# Initial circular curve
ratio_curve = 0.382 #Based on Georgia Tech
rad_c = ratio_curve*rad_t

# Nozzle Length
length_full = rad_t/tan(radians(15)) * (sqrt(ratio_area) - 1 + 1.5*(1/cos(radians(15)) - 1))
length = frac_con*length_full

# End Points of Parabola
coord_x_ce = rad_c*sin(ang_exp)
coord_y_ce = rad_t + rad_c*(1 - cos(ang_exp))

coord_x_e = length
coord_y_e = sqrt(ratio_area)*rad_t

# Parabola Coefficients
tan_f_1 = tan(ang_exp)**2 + 1
# Numerator
numerator = (
    -((coord_x_ce * coord_y_ce**2) / (tan_f_1)) 
    + coord_x_ce * coord_y_ce**2 
    + (2 * coord_x_ce * coord_y_ce * coord_y_e) / (tan_f_1) 
    - 2 * coord_x_ce * coord_y_ce * coord_y_e 
    - (coord_x_ce * coord_y_e**2) / (tan_f_1) 
    + coord_x_ce * coord_y_e**2 
    + (coord_y_ce**3 * sqrt(1 - 1 / (tan_f_1))) / sqrt(tan_f_1) 
    + (coord_y_ce**2 * coord_x_e) / (tan_f_1) 
    - coord_y_ce**2 * coord_x_e 
    - (3 * coord_y_ce**2 * coord_y_e * sqrt(1 - 1 / (tan_f_1))) / sqrt(tan_f_1) 
    - (2 * coord_y_ce * coord_x_e * coord_y_e) / (tan_f_1) 
    + 2 * coord_y_ce * coord_x_e * coord_y_e 
    + (3 * coord_y_ce * coord_y_e**2 * sqrt(1 - 1 / (tan_f_1))) / sqrt(tan_f_1) 
    + (coord_x_e * coord_y_e**2) / (tan_f_1) 
    - coord_x_e * coord_y_e**2 
    - (coord_y_e**3 * sqrt(1 - 1 / (tan_f_1))) / sqrt(tan_f_1))
# Denominator
denominator = (
    (coord_x_ce**2) / (tan_f_1) 
    - coord_x_ce**2 
    - (2 * coord_x_ce * coord_x_e) / (tan_f_1) 
    + 2 * coord_x_ce * coord_x_e 
    + (coord_y_ce**2) / (tan_f_1) 
    - (2 * coord_y_ce * coord_y_e) / (tan_f_1) 
    + (coord_x_e**2) / (tan_f_1) 
    - coord_x_e**2 
    + (coord_y_e**2) / (tan_f_1))
# Calculate Coefficients
a = numerator / denominator
b = (-a*coord_x_ce + a*coord_x_e + (coord_y_ce-coord_y_e)**2)**2 / (2*(coord_y_e - coord_y_ce))**2 - a*coord_x_e
c = coord_y_ce - sqrt(coord_x_ce*a + b)

# Exit Angle
ang_e = atan(a/(2*sqrt(a*coord_x_e + b)))

# Equations
def circle(x):
    eq_c = -np.sqrt(rad_c**2 - x**2) + rad_t + rad_c
    return eq_c

def parabola(x):
    eq_p = np.sqrt(a*x + b) + c
    return eq_p

# Graph Nozzle
x_c = np.linspace(0, coord_x_ce, 1000)
y_c = circle(x_c)

x_p = np.linspace(coord_x_ce, coord_x_e, 1000)
y_p = parabola(x_p)

# Printing Geometric Parameters
print(f"Initial Circular Curve Parameters:\n\tRadius: {rad_c} [mm]\n\tArc Length: {rad_c*ang_exp} [mm]\n\tCenter: (0, {rad_t + rad_c}) [mm]")
print(f"Parabolic Curve Parameters:\n\tEquation in the form: y = (Ax + B)^1/2 + C\n\t\tA = {a}\n\t\tB = {b}\n\t\tC = {c}")
print(f"Initial Expansion Angle: {degrees(ang_exp)} [deg]")
print(f"Exit Expansion Angle: {degrees(ang_e)} [deg]")
print(f"Nozzle Length: {coord_x_e} [mm]")

# Ploting Nozzle Geometry
plt.plot(x_c, y_c)
plt.plot(x_p, y_p)
plt.ylim(bottom=0)
plt.title("Nozzle Wall Curve")
plt.xlabel("Distance from Throat [mm]")
plt.ylabel("Distance from Centerline [mm]")
plt.show()