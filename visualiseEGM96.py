# -*- coding: utf-8 -*-
"""
Created on Sat Feb 28 19:31:24 2015

Plot the gravity potential of the EGM96 gravity field model at a spherical surface.
Ignore tides for simplicity.

EGM96 coefficients available at:
ftp://cddis.gsfc.nasa.gov/pub/egm96/general_info/readme.egm96

Reference for the mathematics:
http://www.ngs.noaa.gov/PUBS_LIB/EGM96_GEOID_PAPER/egm96_geoid_paper.html
    
@author: Alek
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot, matplotlib.cm, matplotlib.ticker, matplotlib.font_manager, scipy.special, numpy, math, csv

"""
===============================================================================
    PROBLEM SETUP.
===============================================================================
"""
GM = 3986004.415E8 # Earth's gravity-mass constant from EGM96, m**3/s**2.
a = 6378136.3 # Equatorial scale factor from EGM96, m.
N_POINTS = 50 # Number of lattitudes and longitudes used to plot the geoid.
latitudes = numpy.linspace(0, 2*numpy.pi, N_POINTS) # Geocentric latitudes and longitudes where the geoid will be visualised.
longitudes = numpy.linspace(0, numpy.pi, N_POINTS)
radius = 6378136.3 # Radius at which the equipotential will be computed, m.
MAX_DEGREE = 4 # Maximum degree of the geopotential to visualise.

" Read the coefficients. "
degrees = []; orders = []; CcoeffsTemp = []; ScoeffsTemp = [];
with open("EGM96coefficients", "r") as egm96file:
    reader = csv.reader(egm96file, delimiter=" ")
    for row in reader:
        degrees.append( row[1] ) # There will be some "  " in row, the delimiter isn't always " ", sometimes it's "  "...
        orders.append( row[2] )
        CcoeffsTemp.append( row[3] )
        ScoeffsTemp.append( row[4] )

degrees=map(int, degrees); orders=map(int, orders); CcoeffsTemp=map(float, CcoeffsTemp); ScoeffsTemp=map(float, ScoeffsTemp); # Change to numbers from str.

" Parse C and S coefficients to an easily usable format. "
# Store a list of coefficients corresponding to the given degree of len( no. orders corresponding to this degree ).
Ccoeffs = {0:[1],1:[0]}; Scoeffs ={0:[0],1:[0]} # Initial coefficients for spherical Earth.
for i in range(len(degrees)): # Initialise emoty lists.
    Ccoeffs[degrees[i]] = []
    Scoeffs[degrees[i]] = []

for i in range(len(degrees)): # Store the coefficients.
    Ccoeffs[degrees[i]].append( CcoeffsTemp[i] )
    Scoeffs[degrees[i]].append( ScoeffsTemp[i] )
    
" Compute the gravitational potential at the desired locations. "
lats,longs = numpy.meshgrid(latitudes, longitudes)
gravitationalPotentials = numpy.zeros( lats.shape ) # Gravitational potentials computed with the given gravity model.

for i in range(lats.shape[0]):
    for j in range(lats.shape[1]):
        for degree in range(0, MAX_DEGREE+1): # Go through all the desired orders and compute the geoid corrections to the sphere.
            temp = 0. # Contribution to the potential from the current degree and all corresponding orders.
            legendreCoeffs = scipy.special.legendre(degree) # Legendre polynomial coefficients corresponding to the current degree.
            for order in range(degree): # Go through all the orders corresponding to the currently evaluated degree.
                temp += legendreCoeffs[order] * numpy.cos(lats[i,j]) * (Ccoeffs[degree][order]*numpy.cos( order*longs[i,j] ) + Scoeffs[degree][order]*numpy.sin( order*longs[i,j] ))
            
            gravitationalPotentials[i,j] += math.pow(a/radius, degree) * temp # Add the contribution from the current degree.
        
gravitationalPotentials *= -GM/radius # Final correction.

"""
===============================================================================
    PLOT FORMATTING.
===============================================================================
"""
ticksFontSize = 15
labelsFontSize = 30
titleFontSize = 34

matplotlib.rc('xtick', labelsize=ticksFontSize) 
matplotlib.rc('ytick', labelsize=ticksFontSize)

"""
===============================================================================
    FIGURE THAT SHOWS THE SPHERICAL HARMONICS.
===============================================================================
"""
fig = matplotlib.pyplot.figure(figsize=(12,8))
ax = Axes3D(fig)
ax.set_aspect("equal")
ax.view_init(elev=45., azim=45.)
ax.set_xlim([-1.5*radius, 1.5*radius])
ax.set_ylim([-1.5*radius, 1.5*radius])
ax.set_zlim([-1.5*radius, 1.5*radius])

gravitationalPotentialsPlot = gravitationalPotentials/gravitationalPotentials.max() # Normalise to [0 1]

" Plot a sphere. "
Xs = radius * numpy.outer(numpy.cos(latitudes), numpy.sin(longitudes))
Ys = radius * numpy.outer(numpy.sin(latitudes), numpy.sin(longitudes))
Zs = radius * numpy.outer(numpy.ones(latitudes.size), numpy.cos(longitudes))
equipotential = ax.plot_surface(Xs, Ys, Zs, facecolors=matplotlib.cm.jet(gravitationalPotentialsPlot), rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)

fig.show()