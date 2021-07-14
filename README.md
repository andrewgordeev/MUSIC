**Notes for gluon test version**:
-----------------------------

Compile the code as with the public version of MUSIC, then run it with input file 'music_input_gluontest' to run the test. This performs a single timestep of evolution on an initial cylindrically symmetric Gaussian temperature profile (max temperature = 600 MeV, width = 5 fm). 

The script 'tests/pce/plot-velocity.py' plots the velocity distribution after one timestep and compares it to the result given by Mathematica, with the plot saved in 'tests/pce/plots.'

The script 'tests/pce/plot-eos.py' plots any quantities from the equation of state, saving plots in the same folder.
