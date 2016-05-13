import numpy as np
from astropy import constants as cst, units as u


def mags(band , temp = 6000, distance = 10, radius = 1):
    # calculates the Johnson-Cousins magnitudes
    # reference values are from at www.astrophysicsspectator.com/topics/observation/MagnitudesAndColors.html
    print("\nEntering Mags Function\n" + "Calculating for", band, temp, sep=" ")    
    zero_flux_erg = [3.98e-5, 6.95e-5, 3.63e-5, 2.254e-5, 1.196e-5]
    website_units = u.erg/(u.s*u.um*u.cm**2)
    flux_units = u.W/(u.m**3)
    zero_flux = website_units.to(flux_units, zero_flux_erg) * flux_units * u.m
    #              center    range
    info = {"UV": (3.50e-7, 7.00e-8, zero_flux[0]), 
            "B":  (4.38e-7, 9.85e-8, zero_flux[1]),
            "V":  (5.465e-7, 8.7e-7, zero_flux[2]),
            "R":  (6.47e-7, 1.515e-7, zero_flux[3]),
            "I":  (7.865e-7, 1.09e-7, zero_flux[4])}
    print ("Ranges: ", info[band], sep=" ")
    wv1 = info[band][0] - info[band][1]/2
    wv2 = info[band][0] + info[band][1]/2

    # simpsons integration technique n=1000
    curve = lambda l,t: 2*cst.h*cst.c**2/( u.m * l )**5  *  1/(np.exp( cst.h*cst.c/(l*cst.k_B*t*u.K*u.m) ) - 1 )
    total = curve(wv1, temp) + curve(wv2, temp)
    h = (wv2 - wv1) / 1000 
    for i in range(1, 1000, 2):
        total += 4 * curve(wv1 + i * h , temp)
    for i in range(2, 999, 2):
        total += 2 * curve(wv1 + i * h , temp)
    total = total * h * u.m / 3
    total *= np.pi # this accounts for directional emmission (sr^-1) )
    print ("Total is ", total)

    # now that total emmission is found, flux recieved is calculated
    # distance is in parsecs and radius is in solar radii
    radius *= u.Rsun.to(u.pc)
    recieved = (radius/distance)**2 * total
    print("Ratio: ", (radius/distance)**2)
    print("Recieved Flux: ", recieved )
    print("Reference Flux: ", info[band][2] * info[band][1])
    mag = -2.5 * np.log10( recieved/(info[band][2]*info[band][1]) )
    print( "Magnitude: ", mag)
    #abs_mag = mag - 5 * np.log10(10/distance)
    return mag

print ("Vega: ", mags("B", 9602, 7.68, 2.362) - mags("V", 9602, 7.68, 2.362))
print ("Betelgeuse: ", mags("B", 3300, 197, 1100) - mags("V", 3300, 197, 1100))
print ("Sun: ", mags("B", 5772, 4.85e-6) - mags("V", 5772, 4.85e-6))
