Introduction
============

- This is a Python project for the simulation the sky-map in ultra-long 
  wavelength.
  The output can be stored either in memory by means of a numpy array, 
  or in a HDF5 format file.
- This software can simply run as a single process on a single compute node,
  but for higher performance, it can also use the Message Passing Interface
  (MPI) to run data processing tasks distributed and parallelly on multiple
  computie nodes \ supercomputers.
- It can be used in upcoming DSL project as input sky map to simulating 
  the observed sky map.
- It can fulfill the processing tasks of fitting spectral index of galactic synchrontron,       interpolate sky map from  408MHz to any lower frequency, and to produce the    skymap in any frequency considering the absorption effect, especially in very-long wavelength. 

