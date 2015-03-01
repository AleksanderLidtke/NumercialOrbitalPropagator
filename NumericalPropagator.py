# -*- coding: utf-8 -*-
"""
Created on Sat Feb 28 22:08:40 2015

@author: alek
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot, matplotlib.cm, matplotlib.ticker, matplotlib.font_manager, scipy.integrate, numpy, csv
import nrlmsise_00_header, nrlmsise_00

def readEGM96Coefficients():
    """ Read the EGM96 gravitational field model coefficients from EGM96coefficients
    file and parse them to be used with computeGravitationalPotential functions.
    
    Returns
    ----------
    2-tuple of the C and S coefficnients of EGM96 model. They are stored in dictionaries
        of list. The keys are degrees of the potential expansion and the values
        of the list entries are the coefficients corresponding to the orders for
        the expansion to a given degree.
    
    Reference
    ----------
    EGM96 coefficients have been downloaded from:
        ftp://cddis.gsfc.nasa.gov/pub/egm96/general_info/readme.egm96
    """
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
    Ccoeffs = {0:[1],1:[0]}; Scoeffs ={0:[0],1:[0]}; # Initial coefficients for spherical Earth.
    for i in range(len(degrees)): # Initialise emoty lists.
        Ccoeffs[degrees[i]] = []
        Scoeffs[degrees[i]] = []
    
    for i in range(len(degrees)): # Store the coefficients.
        Ccoeffs[degrees[i]].append( CcoeffsTemp[i] )
        Scoeffs[degrees[i]].append( ScoeffsTemp[i] )
        

"""
===============================================================================
    PROBLEM SETUP.
===============================================================================
"""
" Gravity model settings. "
GM = 3986004.415E8 # Earth's gravity constant from EGM96, m**3/s**2.
EarthRadius = 6378136.3 # Earth's equatorial radius from EGM96, m.

" Atmospheric density model settings. "
NRLMSISEflags = nrlmsise_00_header.nrlmsise_flags()
NRLMSISEaph = nrlmsise_00_header.ap_array() #TODO this should contain the following:
# * Array containing the following magnetic values:
# *   0 : daily AP
# *   1 : 3 hr AP index for current time
# *   2 : 3 hr AP index for 3 hrs before current time
# *   3 : 3 hr AP index for 6 hrs before current time
# *   4 : 3 hr AP index for 9 hrs before current time
# *   5 : Average of eight 3 hr AP indicies from 12 to 33 hrs 
# *           prior to current time
# *   6 : Average of eight 3 hr AP indicies from 36 to 57 hrs 
# *           prior to current time
F10_7A = 70 # 81-day average F10.7.
F10_7 = 180 # Daily F10.7 for the previous day.
MagneticIndex = 40 # Daily magnetic index.

" Initial conditions of the satellite. "
satelliteMass = 1000. # kg
Cd = 2.2 # Drag coefficient, dimensionless.
dragArea = 5.0 # Area exposed to atmospheric drag, m2.

def calculateDragAcceleration(stateVec, satMass):
    #    " Prepare the atmospheric density model inputs. "
#    #TODO - calculate the altitude, latitude, longitude
#    altitude_km = numpy.linalg.norm(stateVector[:3])/1000.0
#    NRLMSISEinput = nrlmsise_00_header.nrlmsise_input(year=0, doy=0, sec=0.0, alt=altitude_km, g_lat=0.0, g_long=0.0, lst=0.0, f107A=F10_7A, f107=F10_7, ap=MagneticIndex, ap_a=NRLMSISEaph)
#    nrlmsise_00_header.lstCalc( NRLMSISEinput ) # Calculate the local solar time.
#    
#    " Use the calculated atmospheric density to compute the drag force. "
#    NRLMSISEoutpt = nrlmsise_00_header.nrlmsise_output(); nrlmsise_00.gtd7(NRLMSISEinput, NRLMSISEflags, NRLMSISEoutpt);
#    atmosphericDensity = NRLMSISEoutpt.d[5]/1000.0 # Change from gm/cm3 to kg/m3
    dragForce = numpy.zeros(3) #-0.5*atmosphericDensity*dragArea*Cd* numpy.power(stateVector[3:],2)
    return dragForce/satMass

def calculateGravityAcceleration(stateVec):
    """ Calculate the acceleration due to gractiy acting on the satellite at
    a given state (3 positions and 3 velocities).
    Arguments
    ----------
    numpy.ndarray of shape (1,6) with three Cartesian positions and three velocities
        in an inertial reference frame in metres and metres per second, respectively.
    Returns
    ----------
    numpy.ndarray of shape (1,3) with three Cartesian components of the acceleration
        in m/s2 given in an inertial reference frame.
    """
    #TODO need to change to Earth-fixed frame to compute the gravity acceleration, then back to inertial
    gravityAcceleration = -GM/stateVec[:3].dot(stateVec[:3]) * stateVec[:3]/numpy.linalg.norm(stateVec[:3])
    return gravityAcceleration
    
def computeRateOfChangeOfState(stateVector, epoch, satMass):
    """ Compute the rate of change of the state vector.
    Arguments
    ----------
    stateVector - numpy.ndarray of shape (1,6) with three Cartesian positions
        and three velocities given in an inertial frame of reference.
    epoch - float corresponding to the epoch at which the rate of change is to
        be computed.
    satMass - float corresponding to the satellite mass.
    Returns
    ----------
    numpy.ndarray of shape (1,6) with the rates of change of position and velocity
        in the same inertial frame as the one in which stateVector was given.
    """
    gravityAcceleration = calculateGravityAcceleration(stateVector) # A vector of the gravity force from EGM96 model.
    dragAcceleration = calculateDragAcceleration(stateVector, satMass) #  A vector of the drag computed with NRLMSISE-00.
    
    stateDerivatives = numpy.zeros(6);
    stateDerivatives[:3] = stateVector[3:]; # Velocity is the rate of change of position.
    stateDerivatives[3:] = dragAcceleration+gravityAcceleration # Compute the acceleration i.e. the rate of change of velocity.
    return stateDerivatives

def calculateCircularPeriod(stateVec):
    """ Calculate the orbital period of a circular, Keplerian orbit passing through
    the state vector (3 positions and velocities).
    Arguments
    ----------
    numpy.ndarray of shape (1,3) with three Cartesian positions and velocities,
        in mtres and m/s, respectively.
    Returns
    ----------
    Orbital period of a circular orbit corresponding to the supplied state vector
        in seconds.
    """
    return 2*numpy.pi*numpy.sqrt(numpy.power(numpy.linalg.norm(stateVec[:3]),3)/GM)
    
"""
===============================================================================
    PROPAGATE THE ORBIT NUMERICALLY.
===============================================================================
"""
state_0 = numpy.array([EarthRadius+500.0e3,0.,0.,0.,0.,0.]) # Initial state vector with Cartesian positions and velocities in m and m/s.
state_0[5] = numpy.sqrt( (GM+satelliteMass)/numpy.linalg.norm(state_0[:3]) ) # Assume we're starting from a circular orbit.
initialOrbitalPeriod = calculateCircularPeriod(state_0) # Orbital period of the initial circular orbit.

epochsOfInterest = numpy.arange(0, 10*initialOrbitalPeriod, 10) # Times at which the solution is to be computed. Same unit as the time unit of the accelerations and velocities (seconds).
propagatedStates = scipy.integrate.odeint(computeRateOfChangeOfState, state_0, epochsOfInterest, args=(satelliteMass,)) # State vectors at the epochs of interest.

altitudes = [ (numpy.linalg.norm(x[:3])-EarthRadius) for x in propagatedStates] # Altitudes above spherical Earth...
specificEnergies = [ numpy.linalg.norm(x[3:])*numpy.linalg.norm(x[3:]) - GM*satelliteMass/numpy.linalg.norm(x[:3]) for x in propagatedStates] # ...and corresponding specific orbital energies.

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
    FIGURE THAT SHOWS THE EARTH AND SATELLITE TRAJECTORY.
===============================================================================
"""
fig = matplotlib.pyplot.figure(figsize=(12,8))
ax = Axes3D(fig)
ax.set_aspect("equal")
ax.view_init(elev=45., azim=45.)
ax.set_xlim([-1.5*EarthRadius, 1.5*EarthRadius]); ax.set_ylim([-1.5*EarthRadius, 1.5*EarthRadius]); ax.set_zlim([-1.5*EarthRadius, 1.5*EarthRadius]);

" Plot a sphere that represents the Earth. "
N_POINTS = 50 # Number of lattitudes and longitudes used to plot the geoid.
latitudes = numpy.linspace(0, numpy.pi, N_POINTS) # Geocentric latitudes and longitudes where the geoid will be visualised.
longitudes = numpy.linspace(0, 2*numpy.pi, N_POINTS)
Xs = EarthRadius * numpy.outer(numpy.cos(latitudes), numpy.sin(longitudes))
Ys = EarthRadius * numpy.outer(numpy.sin(latitudes), numpy.sin(longitudes))
Zs = EarthRadius * numpy.outer(numpy.ones(latitudes.size), numpy.cos(longitudes))
earthSurface = ax.plot_surface(Xs, Ys, Zs, rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)

" Plot the trajectory. "
ax.plot(propagatedStates[:,0],propagatedStates[:,1],propagatedStates[:,2], c='r', lw=4)

fig.show()

"""
===============================================================================
    FIGURE THAT SHOWS THE ALTITUDE AND ORBITAL ENERGY EVOLUTION.
===============================================================================
"""
fig2 = matplotlib.pyplot.figure(figsize=(12,8)); ax2=fig2.gca(); ax2_2 = ax2.twinx();
ax2.set_xlabel(r"$Time\ elapsed\ (s)$", fontsize=labelsFontSize)
ax2.set_ylabel(r"$Altitude\ above\ spherical\ Earth\ (m)$", fontsize=labelsFontSize)
ax2_2.set_ylabel(r"$Specific\ orbital\ energy\ (m^2 s^{-2})$", fontsize=labelsFontSize)
ax2.grid(True, which='both')

altPlot=ax2.plot(epochsOfInterest, altitudes, c='r', lw=4)
enPlot=ax2_2.plot(epochsOfInterest, specificEnergies, c='b', lw=4)

fig2.show()