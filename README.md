# NumercialOrbitalPropagator
A numerical orbital propagator written in Python. It can propagate a satellite orbit forward in time by
computing the force acting on the satellite, finding the corresponding acceleration, taking a step in time,
and repeating the process.

It uses the NRLMSISE-00 atmospheric model available at:
https://github.com/DeepHorizons/Python-NRLMSISE-00
and EGM96 geopotential expansion available at:
http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html
