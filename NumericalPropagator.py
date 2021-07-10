# -*- coding: utf-8 -*-
"""
Created on Sat Feb 28 22:08:40 2015

@author: alek
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot, matplotlib.cm, matplotlib.ticker, matplotlib.font_manager, scipy.integrate, numpy, csv, math
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
    
    # Change to numbers from str.
    degrees = [int(x) for x in degrees]
    orders = [int(x) for x in orders]
    CcoeffsTemp = [float(x) for x in CcoeffsTemp]
    ScoeffsTemp = [float(x) for x in ScoeffsTemp]
    
    " Parse C and S coefficients to an easily usable format. "
    # Store a list of coefficients corresponding to the given degree of len( no. orders corresponding to this degree ).
    Ccoeffs = {0:[1],1:[0,0]}; Scoeffs ={0:[0],1:[0,0]}; # Initial coefficients for spherical Earth. C_10, C_11, and S_11 are 0 if the origin is at the geocentre.
    for i in range(len(degrees)): # Initialise emoty lists.
        Ccoeffs[degrees[i]] = []
        Scoeffs[degrees[i]] = []
    
    for i in range(len(degrees)): # Store the coefficients.
        Ccoeffs[degrees[i]].append( CcoeffsTemp[i] )
        Scoeffs[degrees[i]].append( ScoeffsTemp[i] )
        
    return Ccoeffs, Scoeffs
        
def RungeKutta4(X,t,dt,rateOfChangeFunction):
    """ Use fourth-order Runge-Kutta numericla integration method to propagate
    a system of differential equations from state X, expressed at time t, by a
    time increment dt. Evolution of the sytstem is given by the rateOfChangeFunction
    that gives the rates of change of all the components of state X.
    Arguments
    ----------
    X - numpy.ndarray of shape (1,6) with three Cartesian positions and three velocities
        in an inertial reference frame in metres and metres per second, respectively.
    t - float, epoch at which state X is defined
    dt - float, epoch increment to which the state X is to be propagated
    rateOfChangeFunction - function that returns a numpy.ndarray of shape (1,3)
        with three Cartesian components of the acceleration in m/s2 given in an
        inertial reference frame. Its arguments are the state X and epoch t.
    Returns
    ----------
    numpy.ndarray of shape (1,6) with three Cartesian positions and three velocities
        in an inertial reference frame in metres and metres per second, respectively,
        propagated to time t+dt.
    """
    dxdt = rateOfChangeFunction(X,t)
    k0 = dt*dxdt
    k1 = dt*rateOfChangeFunction(X+k0/2.,t+dt/2.)
    k2 = dt*rateOfChangeFunction(X+k1/2.,t+dt/2.)
    k3 = dt*rateOfChangeFunction(X+k2,t+dt)
    return X + (k0+2.*k1+2.*k2+k3)/6.
    
"""
===============================================================================
    PROBLEM SETUP.
===============================================================================
"""
" Gravity model settings. "
GM = 3986004.415E8 # Earth's gravity constant from EGM96, m**3/s**2.
EarthRadius = 6378136.3 # Earth's equatorial radius from EGM96, m.
MAX_DEGREE = 2 # Maximum degree of the geopotential harmocic expansion to use.
Ccoeffs, Scoeffs = readEGM96Coefficients() # Get the gravitational potential exampnsion coefficients.

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

"""
===============================================================================
    PROPAGATION FUNCTIONS.
===============================================================================
"""
def calculateGeocentricLatLon(stateVec, epoch):
    """ Calculate the geocentric co-latitude (measured from the noth pole not
    equator), longitude and radius corresponding to the state vector given in
    inertial frame at a certain epoch.
    Arguments
    ----------
    stateVec - numpy.ndarray of shape (1,6) with three Cartesian positions and
        three velocities in an inertial reference frame in metres and metres
        per second, respectively.
    epoch - float with the UT epoch corresponding to the stateVect.
    Returns
    ----------
    3-tuple of floats with geocentric latitude, longitude and radius in radians
        and distance units of stateVec.
    References
    ----------
    Conversions taken from:
    http://agamenon.tsc.uah.es/Asignaturas/it/rd/apuntes/RxControl_Manual.pdf
    """
    #TODO need to change to Earth-fixed frame to compute the gravity acceleration, then back to inertial
    r = numpy.linalg.norm(stateVec[:3])
    colat = math.pi/2.0 - stateVec[2]/r
    lon = math.atan( stateVec[1]/stateVec[0] )
    return colat,lon,r

def calculateDragAcceleration(stateVec, ecpoh, satMass):
    """ Calculate the acceleration due to atmospheric draf acting on the
    satellite at a given state (3 positions and 3 velocities) and epoch.
    Use NRLMSISE2000 atmospheric model with globally defined solar activity
    proxies:
        F10_7A - 81-day average F10.7.
        F10_7 - daily F10.7 for the previous day.
        MagneticIndex - daily magnetic index AP.
        NRLMSISEaph - nrlmsise_00_header.ap_array with magnetic values.#
    
    Arguments
    ----------
    numpy.ndarray of shape (1,6) with three Cartesian positions and three
        velocities in an inertial reference frame in metres and metres per
            second, respectively.
    epoch - float corresponding to the epoch at which the rate of change is
        to be computed.
    Returns
    ----------
    numpy.ndarray of shape (1,3) with three Cartesian components of the
        acceleration in m/s2 given in an inertial reference frame.
    """
    #    " Prepare the atmospheric density model inputs. "
#    #TODO - calculate the altitude, latitude, longitude
    altitude_km = numpy.linalg.norm(stateVec[:3])/1000.0 #TODO this isn't altitude in km, but radius in km. Is this OK?
    NRLMSISEinput = nrlmsise_00_header.nrlmsise_input(year=0, doy=0, sec=0.0, alt=altitude_km, g_lat=0.0, g_long=0.0, lst=0.0, f107A=F10_7A, f107=F10_7, ap=MagneticIndex, ap_a=NRLMSISEaph)
    nrlmsise_00_header.lstCalc( NRLMSISEinput ) # Calculate the local solar time.
    
    " Use the calculated atmospheric density to compute the drag force. "
    NRLMSISEoutpt = nrlmsise_00_header.nrlmsise_output(); nrlmsise_00.gtd7(NRLMSISEinput, NRLMSISEflags, NRLMSISEoutpt);
    atmosphericDensity = NRLMSISEoutpt.d[5]/1000.0 # Change from gm/cm3 to kg/m3
    dragForce = -0.5*atmosphericDensity*dragArea*Cd* numpy.power(stateVec[3:],2) # Drag foce in Newtons.
    return dragForce/satMass

def calculateGravityAcceleration(stateVec, epoch):
    """ Calculate the acceleration due to gravtiy acting on the satellite at
    a given state (3 positions and 3 velocities). Ignore satellite's mass,
    i.e. use a restricted two-body problem.
    Arguments
    ----------
    numpy.ndarray of shape (1,6) with three Cartesian positions and three
        velocities in an inertial reference frame in metres and metres per
            second, respectively.
    epoch - float corresponding to the epoch at which the rate of change is
        to be computed.
    Returns
    ----------
    numpy.ndarray of shape (1,3) with three Cartesian components of the
        acceleration in m/s2 given in an inertial reference frame.
    """
    " Compute geocentric latitude and longitude. "
    colatitude,longitude,geocentricRadius = calculateGeocentricLatLon(stateVec, epoch)
    
    " Find the gravitational potential at the desired point. "
    gravitationalPotential = 0.0 # Potential of the gravitational field at the stateVec location.
    for degree in range(0, MAX_DEGREE+1): # Go through all the desired orders and compute the geoid corrections to the sphere.
        temp = 0. # Contribution to the potential from the current degree and all corresponding orders.
        legendreCoeffs = scipy.special.legendre(degree) # Legendre polynomial coefficients corresponding to the current degree.
        for order in range(degree+1): # Go through all the orders corresponding to the currently evaluated degree.
            if (abs(colatitude-math.pi/2. <= 1E-16)) or (abs(colatitude-3*math.pi/2. <= 1E-16)): # We're at the equator, cos(colatitude) will be zero and things will break.
                temp += legendreCoeffs[order] *         1.0          * (Ccoeffs[degree][order]*math.cos( order*longitude ) + Scoeffs[degree][order]*math.sin( order*longitude ))
            else:
                temp += legendreCoeffs[order] * math.cos(colatitude) * (Ccoeffs[degree][order]*math.cos( order*longitude ) + Scoeffs[degree][order]*math.sin( order*longitude ))
                
        gravitationalPotential += math.pow(EarthRadius/numpy.linalg.norm(stateVec[:3]), degree) * temp # Add the contribution from the current degree.
        
    gravitationalPotential *= GM/EarthRadius # Final correction.

    " Compute the acceleration due to the gravity potential at the given point. "
    gravityAcceleration = gravitationalPotential/numpy.linalg.norm(stateVec[:3]) * (-1.*stateVec[:3]/numpy.linalg.norm(stateVec[:3])) # First divide by the radius to get the acceleration value, then get the direction (towards centre of the Earth).

    return gravityAcceleration
    
def computeRateOfChangeOfState(stateVector, epoch):
    """ Compute the rate of change of the state vector.
    Arguments
    ----------
    stateVector - numpy.ndarray of shape (1,6) with three Cartesian positions
        and three velocities given in an inertial frame of reference.
    epoch - float corresponding to the epoch at which the rate of change is
        to be computed.
    Returns
    ----------
    numpy.ndarray of shape (1,6) with the rates of change of position and velocity
        in the same inertial frame as the one in which stateVector was given.
    """
    gravityAcceleration = calculateGravityAcceleration(stateVector, epoch) # A vector of the gravity force from EGM96 model.
    dragAcceleration = calculateDragAcceleration(stateVector, epoch, satelliteMass) #  A vector of the drag computed with NRLMSISE-00.
    
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
    return 2*math.pi*numpy.sqrt(math.pow(numpy.linalg.norm(stateVec[:3]),3)/GM)
    
"""
===============================================================================
    PROPAGATE THE ORBIT NUMERICALLY.
===============================================================================
"""
" Initial conditions. "
state_0 = numpy.array([EarthRadius+500.0e3,0.,0.,0.,0.,0.]) # Initial state vector with Cartesian positions and velocities in m and m/s.
state_0[5] = numpy.sqrt( (GM+satelliteMass)/numpy.linalg.norm(state_0[:3]) ) # Assume we're starting from a circular orbit.
initialOrbitalPeriod = calculateCircularPeriod(state_0) # Orbital period of the initial circular orbit.

Ccoeffs = {0:[1],1:[0,0],2:[-0.484165371736E-03,0,0],3:[0,0,0,0]};
Scoeffs = {0:[0],1:[0,0],2:[0,0,0],3:[0,0,0,0]}
colatitude,longitude,geocentricRadius = calculateGeocentricLatLon(state_0, 0)

gravitationalPotential = 0.0 # Potential of the gravitational field at the stateVec location.
for degree in range(0, MAX_DEGREE+1): # Go through all the desired orders and compute the geoid corrections to the sphere.
        temp = 0. # Contribution to the potential from the current degree and all corresponding orders.
        legendreCoeffs = scipy.special.legendre(degree) # Legendre polynomial coefficients corresponding to the current degree.
        for order in range(degree+1): # Go through all the orders corresponding to the currently evaluated degree.
            if colatitude-math.pi/2. <= 1E-16: # We're at the equator, cos(colatitude) will be zero and things will break.
                temp += legendreCoeffs[order] *         1.0          * (Ccoeffs[degree][order]*math.cos( order*longitude ) + Scoeffs[degree][order]*math.sin( order*longitude ))
            else:
                temp += legendreCoeffs[order] * math.cos(colatitude) * (Ccoeffs[degree][order]*math.cos( order*longitude ) + Scoeffs[degree][order]*math.sin( order*longitude ))
        
        gravitationalPotential += math.pow(EarthRadius/geocentricRadius, degree) * temp # Add the contribution from the current degree.
    
gravitationalPotential *= GM/geocentricRadius # Final correction.
gravityAcceleration = gravitationalPotential/geocentricRadius * (-1.*state_0[:3]/numpy.linalg.norm(state_0[:3])) # First divide by the radius to get the acceleration value, then get the direction (towards centre of the Earth).

" Propagation time settings. "
INTEGRATION_TIME_STEP_S = 10.0 # Time step at which the trajectory will be propagated.
epochsOfInterest = numpy.arange(0, 10*initialOrbitalPeriod, INTEGRATION_TIME_STEP_S) # Times at which the solution is to be computed. Same unit as the time unit of the accelerations and velocities (seconds).
propagatedStates = numpy.zeros( (epochsOfInterest.shape[0],6) ) # State vectors at the  epochs of interest.

" Actual numerical propagation main loop. "
propagatedStates[0,:] = state_0 # Apply the initial condition.
for i in range(1,epochsOfInterest.shape[0]): # Propagate the state to all the desired epochs statring from state_0.
    propagatedStates[i,:] = RungeKutta4(propagatedStates[i-1], epochsOfInterest[i-1],
                                        INTEGRATION_TIME_STEP_S, computeRateOfChangeOfState)
    #TODO check if altitude ins't too low. For some reason the object will start gaining energy 
    
" Compute quantities derived from the propagated state vectors. "
altitudes = [ (numpy.linalg.norm(x[:3])-EarthRadius) for x in propagatedStates] # Altitudes above spherical Earth...
specificEnergies = [ numpy.linalg.norm(x[3:])*numpy.linalg.norm(x[3:]) -
                    GM*satelliteMass/numpy.linalg.norm(x[:3]) for x in propagatedStates] # ...and corresponding specific orbital energies.

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
ax.set_aspect('auto') #TODO change 3D axes aspect ratio to equal, which isn't supported now. Current workaround is set scale_xyz below.
ax.view_init(elev=45., azim=45.)
figRange = 1.5*EarthRadius
ax.set_xlim([-figRange, figRange])
ax.set_ylim([-figRange, figRange])
ax.set_zlim([-figRange, figRange])
ax.auto_scale_xyz([-figRange, figRange], [-figRange, figRange], [-figRange, figRange])

" Plot a sphere that represents the Earth. "
N_POINTS = 20 # Number of lattitudes and longitudes used to plot the geoid.
latitudes = numpy.linspace(0, math.pi, N_POINTS) # Geocentric latitudes and longitudes where the geoid will be visualised.
longitudes = numpy.linspace(0, 2*math.pi, N_POINTS)
Xs = EarthRadius * numpy.outer(numpy.cos(latitudes), numpy.sin(longitudes))
Ys = EarthRadius * numpy.outer(numpy.sin(latitudes), numpy.sin(longitudes))
Zs = EarthRadius * numpy.outer(numpy.ones(latitudes.size), numpy.cos(longitudes))
earthSurface = ax.plot_surface(Xs, Ys, Zs, rstride=1, cstride=1, linewidth=0,
                               antialiased=False, shade=False, alpha=0.5)

" Plot the trajectory. "
ax.plot(propagatedStates[:,0],propagatedStates[:,1],propagatedStates[:,2], c='r', lw=4)

fig.show()

"""
===============================================================================
    FIGURE THAT SHOWS THE ALTITUDE AND ORBITAL ENERGY EVOLUTION.
===============================================================================
"""
fig2 = matplotlib.pyplot.figure(figsize=(12,8))
ax2=fig2.gca()
ax2_2 = ax2.twinx();
ax2.set_xlabel(r"$Time\ elapsed\ (s)$", fontsize=labelsFontSize)
ax2.set_ylabel(r"$Altitude\ above\ spherical\ Earth\ (m)$", fontsize=labelsFontSize)
ax2_2.set_ylabel(r"$Specific\ orbital\ energy\ (m^2 s^{-2})$", fontsize=labelsFontSize)
ax2.grid(True, which='both')

ax2.plot(epochsOfInterest, altitudes, c='r', lw=4, ls='-')
ax2_2.plot(epochsOfInterest, specificEnergies, c='b', lw=4, ls='-')

fig2.show()

"""
===============================================================================
    FIGURE SHOWING EVOLUTION OF THE POSITION COMPONENTS OVER TIME.
===============================================================================
"""
fig3, axarr = matplotlib.pyplot.subplots(3, sharex=True, figsize=(12,8))
axarr[0].grid(linewidth=2); axarr[1].grid(linewidth=2); axarr[2].grid(linewidth=2);
axarr[0].tick_params(axis='both',reset=False,which='both',length=5,width=1.5)
axarr[1].tick_params(axis='both',reset=False,which='both',length=5,width=1.5)
axarr[2].tick_params(axis='both',reset=False,which='both',length=5,width=1.5)

axarr[2].set_xlabel(r'$Time\ elapsed\ (s)$',fontsize=labelsFontSize)
axarr[0].set_ylabel(r'$X\ (m)$',fontsize=labelsFontSize)
axarr[1].set_ylabel(r'$Y\ (m)$',fontsize=labelsFontSize)
axarr[2].set_ylabel(r'$Z\ (m)$',fontsize=labelsFontSize)

axarr[0].plot(epochsOfInterest, propagatedStates[:,0], c='r', lw=4, ls='-')
axarr[1].plot(epochsOfInterest, propagatedStates[:,1], c='r', lw=4, ls='-')
axarr[2].plot(epochsOfInterest, propagatedStates[:,2], c='r', lw=4, ls='-')

fig3.show()