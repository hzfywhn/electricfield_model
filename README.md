Spatial modeling of electric fields using Lattice Kriging

Before running the code, you should have:
1. SuperDARN gridded velocity in text grid2 file, contact Virginia Tech SuperDARN group for the required input files
2. Weimer potential output (text file), note that the Weimer model can be found at https://zenodo.org/record/2530324

Then run codes in the following sequence:
1. txt2nc_ion_drift, to convert SuperDARN text files to netcdf (20150317.north.grd2 -> ion_drift_north.nc)
2. convert_ion_drift_to_E, to convert ion drifts to electric fields and calculate mlt (ion_drift_north.nc -> E_north.nc)
3. txt2nc_potential, to convert Weimer text files to netcdf (potentials.txt -> weimer.nc)
4. extend_potential, to extend Weimer potential to lower latitudes (weimer.nc -> weimer_ext.nc)
5. prepare_input, to combine SuperDARN electric fields and extended Weimer potentials and produce formated input for modeling (north_in.nc, north_out_grid.nc)
6. electricfield_model, to predict potentials from existing electric field observations (north_out_5d.nc)

The final output north_out_5d.nc is the one containing predicted potentials at high latitudes