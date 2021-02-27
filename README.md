# propagatormethods
Matlab scripts for analytic calculations and kinetic Monte Carlo simulations of diffusing particles on networks.

----------------

This code is associated with Zubenelgenubi C. Scott, Aidan I. Brown, Saurabh S. Mogre, Laura M. Westrate, Elena F. Koslover. *Diffusive search and trajectories on spatial networks: a propagator approach*. Manuscript submitted for publication at EPJE.

Please cite the article above if using the code.

Any questions should be directed to Elena Koslover (ekoslover@physics.ucsd.edu).

--------------


You need to provide a network structure file containing the following:
- lines of format "NODE ID x y" or "NODE ID x y z" for 2D or 3D networks, respectively
- lines of format "EDGE ID N1 N2 len" where ID is an integer for the edge ID, N1 and N2 are IDs of nodes connected by the edge. "len" is the edge length (optional, can be omitted). If no edge length is supplied for a given edge, the straight-line distance between the two points is calculated.

An example network structure file is provided in example.net

Run example.m to calculate the mean first passage time and related quantities for several different scenarios, all using the example network:

1. Analytic MFPT in the case of two target nodes. This gives the MFPT to hit the first of the targets for a diffusing particle. The resulting MFPT values for each starting node are given in the array MFPTs.
2. Analytic MFPT in the case of 3 partially absorbing edges. This example gives the MFPT for a diffusing particle to be absorbed by the finite reaction rate edges. The resulting MFPT values for each starting edge are given in the array MFPTs.
3. Simulated MFPT with two target nodes (same as 1, so should be comparable.) This method uses the kinetic Monte Carlo simulations.
4. Simulated diffusing particle trajectories on the network with no targets. Showcases the ability to save particle positions at pre-specified save times and uses edge propagation when necessary.
5. Simulated particle pair reaction times. This gives the time to first encounter for two particles diffusing on the network.


