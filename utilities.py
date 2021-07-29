#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 14:27:17 2021

@author: alek
"""
import numpy, math, matplotlib.pyplot
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

GM = 3986004.415E8 # Earth's gravity constant from EGM96, m**3/s**2.

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

def inertialToVNC(position, velocity):
    """ Convert position and corresponding velocity in an inertial frame of reference (TEME, ECI or similar) to 
    Velocity-Normal-Co-normal frame.\n
    Parameters
    ----------
    @param position, velcity - dimensional position and velocity in the inertial frame. Must have the same units.\n
    Returns
    ----------
    @return - 3x3 rotation matrix from inertial to VNC frames - DOES NOT INCLUDE TRANSLATION BETWEEN ORIGINS.\n"""
    velocityUnit_inertial = velocity/math.sqrt(velocity.dot(velocity)) # Unit velocity vector expressed in inertial frame.
    normalUnit_inertial = numpy.cross (position, velocityUnit_inertial)
    normalUnit_inertial = normalUnit_inertial/math.sqrt(normalUnit_inertial.dot(normalUnit_inertial)) # Normal unit vector expressed in inertial frame.
    coNormalUnit_inertial = numpy.cross( normalUnit_inertial, velocityUnit_inertial)
    normalUnit_inertial = normalUnit_inertial/math.sqrt(normalUnit_inertial.dot(normalUnit_inertial)); # Co-normal (along radial) unit vector expressed in inertial frame.
    
    xAxisUnitVectorInertial=numpy.array([1.,0.,0.]) # Unit vecotrs of the inertial reference frame.
    yAxisUnitVectorInertial=numpy.array([0.,1.,0.])
    zAxisUnitVectorInertial=numpy.array([0.,0.,1.])

    row0 = numpy.array([velocityUnit_inertial.dot(xAxisUnitVectorInertial), velocityUnit_inertial.dot(yAxisUnitVectorInertial), velocityUnit_inertial.dot(zAxisUnitVectorInertial)])
    row1 = numpy.array([normalUnit_inertial.dot(xAxisUnitVectorInertial), normalUnit_inertial.dot(yAxisUnitVectorInertial), normalUnit_inertial.dot(zAxisUnitVectorInertial)])
    row2 = numpy.array([coNormalUnit_inertial.dot(xAxisUnitVectorInertial), coNormalUnit_inertial.dot(yAxisUnitVectorInertial), coNormalUnit_inertial.dot(zAxisUnitVectorInertial)])
    
    rotationMatrix = numpy.vstack( (row0, row1, row2) ) # Form a full 3x3 transformation matrix.
    return rotationMatrix

def inertialToRTC(position, velocity):
    """ Convert position and corresponding velocity in an inertial frame of reference (TEME, ECI or similar) to 
    Radial - Transverse - In-track (also called Radial - Along-track - Cross-ctrack) frame.\n
    Parameters
    ----------
    @param position, velcity - dimensional position and velocity in the inertial frame. Must have the same units.\n
    Returns
    ----------
    @return - 3x3 rotation matrix from inertial to RTS frames - DOES NOT INCLUDE TRANSLATION BETWEEN ORIGINS.\n"""
    radialUnit_inertial = position/math.sqrt(position.dot(position)) # Unit radial vector expressed in inertial frame.
    crossUnit_inertial = numpy.cross(position, velocity)
    crossUnit_inertial = crossUnit_inertial/math.sqrt(crossUnit_inertial.dot(crossUnit_inertial)) # Cross-track unit vector expressed in the inertial reference frame.
    transverseUnit_inertial = numpy.cross(radialUnit_inertial, crossUnit_inertial)
    transverseUnit_inertial = transverseUnit_inertial/math.sqrt(transverseUnit_inertial.dot(transverseUnit_inertial)) # Transverse unit vector expressed in the inertial reference frame.

    xAxisUnitVectorInertial=numpy.array([1.,0.,0.]) # Unit vecotrs of the inertial reference frame.
    yAxisUnitVectorInertial=numpy.array([0.,1.,0.])
    zAxisUnitVectorInertial=numpy.array([0.,0.,1.])

    row0 = numpy.array([radialUnit_inertial.dot(xAxisUnitVectorInertial), radialUnit_inertial.dot(yAxisUnitVectorInertial), radialUnit_inertial.dot(zAxisUnitVectorInertial)])
    row1 = numpy.array([transverseUnit_inertial.dot(xAxisUnitVectorInertial), transverseUnit_inertial.dot(yAxisUnitVectorInertial), transverseUnit_inertial.dot(zAxisUnitVectorInertial)])
    row2 = numpy.array([crossUnit_inertial.dot(xAxisUnitVectorInertial), crossUnit_inertial.dot(yAxisUnitVectorInertial), crossUnit_inertial.dot(zAxisUnitVectorInertial)])
    
    rotationMatrix = numpy.vstack( (row0, row1, row2) ) # Form a full 3x3 transformation matrix.
    del row0, row1, row2
    return rotationMatrix

def getCovarianceEllipsePoints(covarianceMatrix, noPoints=100):
    """ Generate x and y points that define an uncertainty ellipse, given a 2x2
    covarianceMatrix.\n
    Parameters
    ----------
    @param covarianceMatrix - 2x2 covariance matrix.\n
    @param maxRadius - maximum radius allowed for the ellipse.\n
    @param noPoints - number of points which will be generated for the ellipse.\n
    Returns
    ----------
    @return (x,y,z) - 3-tuple of points lying on the ellispe, all of length noPoints, in numpy.array format. z has default value of 0 througout"""
    angles=numpy.linspace(0., 2*math.pi, noPoints) # Angles around a circle.
    
    eigenValuesTEMP, eigenVectors = numpy.linalg.eig (covarianceMatrix ) # Compute eigen-stuff
    eigenValues=numpy.array([ [eigenValuesTEMP[0],0.],[0.,eigenValuesTEMP[1]]])
    xy = numpy.hstack( (numpy.cos(angles.reshape(-1,1)),numpy.sin(angles.reshape(-1,1))) ).dot( numpy.sqrt(eigenValues) ).dot(eigenVectors.T)
    x = xy[:,0]
    y = xy[:,1]
    z = numpy.zeros(len(x))
    return (x,y,z)
    
def plotCovarianceEllipsoid(covarianceMatrix, centreLocation=numpy.array([0.,0.,0.]), axes3d=None, noPoints=100, colour='b'):
    """ Plots a covariance ellipsoid in the given figure. Include the 1-sigma and 3-sigma probability levels
    with varying degree of translucence.\n
    Parameters
    ----------
    @param covarianceMatrix - 3x3 numpy.array which descrbes the covaraince ellipsoid.\n
    @param centreLocation - location of the centre point of the ellipsoid. Assumed to be the origin.\n
    @param axes3d - axes of the figure where the ellipsoid is to be plotted. If None a new figure will be created.\n
    @param noPoints - number of points to be found on the ellispe in each dimension (noPoints^3 in total). Keep it low to conserve memory but high enough to get the desired shape accuracy.\n
    @param colour - coulour designator supported by matplotlib.Axes.plot_wireframe that will be used to plot the ellipsoid surface."""
    Cxy=covarianceMatrix[0:2,0:2]
    Cyz=covarianceMatrix[1:3,1:3]
    Czx=numpy.array([[covarianceMatrix[2,2],covarianceMatrix[2,0]],[covarianceMatrix[0,2],covarianceMatrix[0,0]]])
    
    if axes3d is None: # Create a figure if neccesarry and ge the axes.
        fig3dCovarianceEllipsoid=matplotlib.pyplot.figure()
        axes3d=fig3dCovarianceEllipsoid.gca(projection='3d')
        
    U, s, rotation = numpy.linalg.svd(covarianceMatrix)
    radii = numpy.sqrt(s)
    radii3 = 3*numpy.sqrt(s)
    u = numpy.linspace(0.0, 2.0 * numpy.pi, noPoints) # Find the coordinates which are evenly spaced around the whole ellipsoid.
    v = numpy.linspace(0.0, numpy.pi, noPoints)
    x = radii[0] * numpy.outer(numpy.cos(u), numpy.sin(v)) # And compute the locations of the points at those coordinates for 1-sigma.
    y = radii[1] * numpy.outer(numpy.sin(u), numpy.sin(v))
    z = radii[2] * numpy.outer(numpy.ones_like(u), numpy.cos(v))
    x3 = radii3[0] * numpy.outer(numpy.cos(u), numpy.sin(v)) # And 3-sigma.
    y3 = radii3[1] * numpy.outer(numpy.sin(u), numpy.sin(v))
    z3 = radii3[2] * numpy.outer(numpy.ones_like(u), numpy.cos(v))
    for i in range(len(x)): # Also apply a translation and rotation to those points.
        for j in range(len(x)):
            [x[i,j],y[i,j],z[i,j]] = numpy.dot([x[i,j],y[i,j],z[i,j]], rotation) + centreLocation
            [x3[i,j],y3[i,j],z3[i,j]] = numpy.dot([x3[i,j],y3[i,j],z3[i,j]], rotation) + centreLocation
    
    xX,yY,zeros=getCovarianceEllipsePoints(Cxy) # Compute and plot the ellipsoid in the x-y plane.
    axes3d.plot(xX+centreLocation[0],yY+centreLocation[1],zeros+centreLocation[2],colour,label='Cxy')
    
    xY,yZ,zeros=getCovarianceEllipsePoints(Cyz) # Now the y-z plane...
    axes3d.plot(zeros+centreLocation[0],xY+centreLocation[1],yZ+centreLocation[2],colour,label='Cyz')
    
    xZ,yX,zeros=getCovarianceEllipsePoints(Czx) # ...and, finally, the x-z plane.
    axes3d.plot(yX+centreLocation[0],zeros+centreLocation[1],xZ+centreLocation[2],colour,label='Czx')

    axes3d.plot_wireframe(x, y, z,  rstride=4, cstride=4, color=colour, alpha=0.5) # Plot the 1-sigma ellipsoid surface.
    axes3d.plot_wireframe(x3, y3, z3,  rstride=4, cstride=4, color=colour, alpha=0.1) # Also include the 3-sigma one.

def coe2rv(sma,ecc,incl,raan,argp,ta):
    """ Convert classical orbital elements into Cartesian ones in the same frame
    of reference.
    Parameters
    ----------
    * sma - semi-major axis in m
    * ecc - eccentricity
    * incl - inclination in degrees
    * raan - right ascension of ascending node in degrees
    * argp - argument of perigee in degrees
    * ta - true anomaly in degrees
    
    Returns
    ----------
    numpy.ndarray of shape (6,1) with three Cartesian positions and three velocities
        in km and km/sec, respectively.
    """
    slr = sma*(1.-ecc*ecc) # Semilatus rectum in km
    
    # Check that RAAN and true anomaly are in 0 to 360 deg format.
    if raan < 0: raan = raan+360
    if ta < 0: ta = ta+360
        
    incl = numpy.deg2rad(incl) # Need these in radians.
    raan = numpy.deg2rad(raan)
    argp = numpy.deg2rad(argp)
    ta = numpy.deg2rad(ta)
        
    # Postion and velocity vectors in the perifocal frame.
    posPF = numpy.zeros(3); velPF = numpy.zeros(3);
    posPF = (slr/(1.+ecc*numpy.cos(ta))) * numpy.array([numpy.cos(ta),numpy.sin(ta),0])
    velPF = numpy.sqrt(GM/slr) * numpy.array([-numpy.sin(ta),ecc+numpy.cos(ta),0])
    
    # Rotate perifocal to inertial frame.
    # First rotation about the angular momentum vector (angle: argument of perigee)
    ROT3_w = numpy.array([[ numpy.cos(-argp), numpy.sin(-argp), 0.],
                          [-numpy.sin(-argp), numpy.cos(-argp), 0.],
                          [               0.,                0, 1.]])
    
    # Second rotation about the line of nodes vector (angle: inclination).
    ROT1_i = numpy.array([[1.,                0.,               0.],
                          [0.,  numpy.cos(-incl), numpy.sin(-incl)],
                          [0., -numpy.sin(-incl), numpy.cos(-incl)]])
    
    # Thirds rotation about Z ECI frame axis (angle: right ascension)
    ROT3_raan = numpy.array([[ numpy.cos(-raan), numpy.sin(-raan), 0.],
                             [-numpy.sin(-raan), numpy.cos(-raan), 0.],
                             [               0.,                0, 1.]])
               
    # Combined rotation matrix
    pf2ijk = ROT3_raan.dot(ROT1_i).dot(ROT3_w);

    return numpy.hstack( (pf2ijk.dot(posPF), pf2ijk.dot(velPF)) )