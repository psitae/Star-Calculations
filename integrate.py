# The following code is a test program to try to find the error in the magnitude calculations

from astropy import analytic_functions as fxn
from astropy import constants as cst, units as u
import numpy as np

np.errstate(over='ignore')
wv1 = 100e-9
wv2 = 100e-6
temp = int(input())

# simpson's integration technique n=1000
curve = lambda l,t: 2*cst.h*cst.c**2/( u.m * l )**5  *  1/(np.exp( cst.h*cst.c/(l*cst.k_B*t*u.K*u.m) ) - 1 )
total = curve(wv1, temp) + curve(wv2, temp)
h = (wv2 - wv1) / 1000 
for i in range(1, 1000, 2):
    total += 4 * curve(wv1 + i * h , temp)
for i in range(2, 999, 2):
    total += 2 * curve(wv1 + i * h , temp)
total = total * h * u.m / 3

# rectangular integration
total2 = 0
for i in np.linspace(wv1, wv2, 1000):
    total2 += curve(i, temp)*h*u.m

# this is to test the validity of the curve function
first = 2e-8
one = curve(1e-8, temp)
two = curve(2e-8, temp)
while one < two:
    first += 1e-8
    one, two = two, curve(first, temp)
 

print("Simpson's Integration: ", total)
print("Rectagle Integration: ", total2)
print("Stefan-Boltzmann law: ", cst.sigma_sb*(temp*u.K)**4)
print("Peak at ", first)
print("Peak should be at ", cst.b_wien/(temp*u.K))
print("curve value 502nm", curve (502e-9, temp)/1e12, "Trillion")
print("Exact: ", fxn.blackbody_lambda(5020, temp).to(u.W/(u.m**3*u.sr)))
print("curve value 150nm", curve (150e-6, temp)/1e12, "Trillion")
print("Exact: ", fxn.blackbody_lambda(150*u.um, temp).to(u.W/(u.m**3*u.sr)))
