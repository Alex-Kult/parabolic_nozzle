from math import *

def prandtl_meyer(gamma, mach_e):
    nu = sqrt((gamma + 1)/(gamma - 1))*atan(sqrt((gamma - 1)/(gamma + 1)*(mach_e**2-1))) - atan(sqrt(mach_e**2 - 1))
    return nu