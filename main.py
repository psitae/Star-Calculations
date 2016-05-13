#!/usr/bin/env python      

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from  matplotlib.figure import Figure
import numpy as np
from astropy import constants as cst, units as u
from tkinter import *


class App():
    
    def make_graph(self, temp):
        
        temp = int(temp)
        mu = np.linspace(1e-7, 2e-6, 1e3, endpoint=True)
        flux = 2*cst.h*cst.c**2/( u.m * mu )**5  *  1/(np.exp( cst.h*cst.c/(mu*cst.k_B*temp*u.K*u.m) ) - 1 )
        
        # chart
        xax = np.linspace(0, 2, 1e3)
        self.ax.plot(xax, flux)
        self.canvas.show()
        
        # call mags function for 4 bands
        print("Calculating ...")
        global bands 
        for i in range(4):
            self.band_mags[i].set( str(self.mags( bands[i], temp, self.distance.get(), self.radius.get() )) )

    def clear_graph(self):

        for i in range(4):
            self.band_mags[i].set("Cleared") 

        self.f.delaxes(self.ax)
        self.ax = self.f.add_subplot(111, xlabel="Wavelength (um)", ylabel="Flux ( W/(m^3*sr) )")

        self.canvas.show()

    def mags(self, band, temp, distance = 10, radius = 1):

        # calculates the Johnson-Cousins magnitudes
        # reference values are from at www.astrophysicsspectator.com/topics/observation/MagnitudesAndColors.html

        #print("\nEntering Mags Function\n" + "Calculating for", band, temp, sep=" ")    
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
        #print ("Ranges: ", info[band], sep=" ")
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
        #print ("Total is ", total)

        # now that total emmission is found, flux recieved is calculated
        # distance is in parsecs and radius is in solar radii
        radius *= u.Rsun.to(u.pc)
        recieved = (radius/distance)**2 * total
        #print("Ratio: ", (radius/distance)**2)
        #print("Recieved Flux: ", recieved )
        #print("Reference Flux: ", info[band][2] * info[band][1])
        mag = -2.5 * np.log10( recieved/(info[band][2]*info[band][1]) )
        mag = "%s" % float("%.3g" % mag)
        #print( "Magnitude: ", mag)
        return mag

    def __init__(self, master):
        
        # title
        title = Label(master, height=4, text="Blackbody Graph and Magnitudes")
        
        # frames
        interface = Frame(master)
        self.chart = Frame(master)
        bottom = Frame(master)
        title.grid(columnspan=2)
        self.chart.grid(row=1, column=0)
        interface.grid(row=1, column=1)
        bottom.grid(columnspan=2)
 
        # blank graph
        self.f = Figure(figsize=(5, 4), dpi = 100)
        self.ax = self.f.add_subplot(111, xlabel="Wavelength (um)", ylabel="Flux ( W/(m^3*sr) )")
        self.canvas = FigureCanvasTkAgg(self.f, self.chart)
        self.canvas.show()
        self.canvas._tkcanvas.grid()

 
        # buttons
        self.gbutton = Button(bottom, text="Graph", command= lambda:self.make_graph( e.get() )  )
        self.clearer = Button(bottom, text="Clear", command=self.clear_graph)
        self.quitter = Button(bottom, text="Quit", command=interface.quit)
        self.quitter.grid(column=0)

        self.clearer.grid(row=0, column=1)
        self.gbutton.grid(row=0, column=2)

        # input
        etitle = Label(interface, text="Temperature (K)")
        etitle.grid()
        e = Entry(interface)
        e.grid()
        e.delete(0, END)
        e.insert(0, "5600")
        
        # scales
        self.distance = Scale(interface, label="Distance from Earth (pc)", from_=1, to=200, orient=HORIZONTAL, length=300) 
        self.distance.set(10)
        self.distance.grid()
        self.radius = Scale(interface, label="Radius (Rsol)", from_=.1, to=100, orient=HORIZONTAL, length = 300, resolution=.5)
        self.radius.set(1)
        self.radius.grid()
       
        # check
        # self.box = Checkbutton(frame1, text="Expand", variable=v, onvalue=1)
        # self.box.grid(side=RIGHT)

        # magnitudes
        mag_frame = Frame(interface)
        mag_frame.grid()
        magtitle = Label(mag_frame, text="Magnitudes")
        magtitle.grid(columnspan=20)

        global bands
        self.band_mags = [StringVar() for n in range(4)]
        self.mag_label = [Label(mag_frame, text=bands[n] + ": ").grid(row = 20 + n, column = 0, columnspan=2) for n in range(4)]
        self.mag_vals = [Label(mag_frame, textvariable=self.band_mags[n]).grid(row = 20 + n, column = 2) for n in range(4)]


bands = ["UV", "B", "V", "R"]
root = Tk()
App(root)

root.mainloop()
            
