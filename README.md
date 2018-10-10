# BeamFilmBestFit
The perl script BF_FindBestParameters.pl is used to compare recombination data simulated using MADpatterns (https://mwhite4.github.io/MADpatterns/) to experimental recombination datasets.

Specifically, it will compare the CoC curves, CO frequencies and Event Distributions from the output files of "analyze_events_on_linear_objects" run on recombination data simulated using "crossover_simulation" with the CoC curve, CO frequencies and Event Distribution of the experimental data.

Simulated data are ranked on the fit of the CoC curves, CO frequencies and Event Distributions to those of the experimental data. The best fit simulation is that with the lowest rank sum.

The script requires the user to provide as arguments
 1. The directory containing the simulated and experimental analysis files 
 2. The file name of the experimental analysis file
 3. Whether the experimental data are from bivalents (0) or gametes (1) e.g. cytological data or backcross data)

