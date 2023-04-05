

lin.csv: parameter space with events where objects don't collide (if black-holes)

tempgen.py: Generate temps, find amplitudes, get the strain value at a given time

Gen_duration.ipynb: Find duration of templates
              used tempgen.py to generate templates. Finds the time interval which covers the strain that is 0.95 times the amplitude.
              (result in /Important_plots/Duration_PS.png)
              
Template_slice_solve.ipynb: Fix the issue where chaning the length of times gives different SNR and Horizon Distances
              Found that the problem was that the injection time changed slightly when the template is cropped. So injecting template at different times 
              gives different SNR.
              
              
Horizon.ipynb: Found horizon distances (result in /Important_plots/Horizon_PS.png)
              
