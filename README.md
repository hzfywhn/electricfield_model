Spatial modeling of electric fields using Lattice Kriging

Before running the code, you should have:
1. SuperDARN gridded velocity in text grid2 file. See 20150317.north.grd2 for a sample. Contact Virginia Tech SuperDARN group for the required input files.
2. IMF required for the Weimer model can be obtained from SPDF OMNIWeb (https://omniweb.gsfc.nasa.gov/).
3. Weimer potential output (text file) after running the Weimer model. Weimer model can be found at https://zenodo.org/record/2530324. imfs.txt and potentialgrid.txt are the sample input for gridpotentials.sav, and potentials.txt is the output.

Then run codes in the following sequence:
1. txt2nc_ion_drift, to convert SuperDARN text files to netcdf (20150317.north.grd2 -> ion_drift_north.nc)
2. calculate_mlt_ion_drift, to include magnetic local time information (ion_drift_mlt_north.nc)
3. (optional) plot_ion_drift, to draw vector plots of ion drifts and compare to SuperDARN map velocity tool
4. convert_ion_drift_to_E, to convert ion drifts to electric fields and calculate mlt (ion_drift_mlt_north.nc -> E_north.nc)
5. txt2nc_potential, to convert Weimer text files to netcdf (potentials.txt -> weimer.nc)
6. extend_potential, to extend Weimer potential to lower latitudes (weimer.nc -> weimer_ext.nc)
7. prepare_input, to combine SuperDARN electric fields and extended Weimer potentials and produce formated input for modeling (north_in.nc, north_out_grid.nc)
8. electricfield_model, to predict potentials from existing electric field observations (north_out_5d.nc)

The final output north_out_5d.nc is the one containing predicted potentials at high latitudes