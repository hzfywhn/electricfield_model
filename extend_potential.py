from netCDF4 import Dataset
from numpy import argmin, isclose, exp


# Weimer potential stops at 50 MLAT
# extend the potential pattern to lower latitudes
# input: Weimer potential
# output: extended Weimer potential to a lower latitude


datain = Dataset(filename='weimer.nc')
time = datain['time']
timeunits = time.units
time = time[:].filled()
mlat = datain['mlat'][:].filled()
mlt = datain['mlt'][:].filled()
pot = datain['potential'][:].filled()
datain.close()

ntime = len(time)
nmlat = len(mlat)
nmlt = len(mlt)

# exponential factor for the decay function
# current value is chosen to make the potential decrease fast
# and let the value at 30 MLAT be small (but non-zero)
expfac = 5

for itime in range(ntime):
    for imlt in range(nmlt):
# find the first non-zero element at each MLT from 90 MLAT
        first = argmin(isclose(pot[itime, :, imlt], 0))

# exponential decay outward (to lower latitudes)
        pot[itime, 0: first, imlt] = pot[itime, first, imlt] / exp((mlat[first] - mlat[0: first]) / expfac)

dataout = Dataset(filename='weimer_ext.nc', mode='w')
dataout.createDimension(dimname='time', size=ntime)
dataout.createDimension(dimname='mlat', size=nmlat)
dataout.createDimension(dimname='mlt', size=nmlt)
time_out = dataout.createVariable(varname='time', datatype='i4', dimensions='time')
mlat_out = dataout.createVariable(varname='mlat', datatype='i4', dimensions='mlat')
mlt_out = dataout.createVariable(varname='mlt', datatype='f4', dimensions='mlt')
pot_out = dataout.createVariable(varname='potential', datatype='f4', dimensions=('time', 'mlat', 'mlt'))
time_out.units = timeunits
time_out[:] = time
mlat_out[:] = mlat
mlt_out[:] = mlt
pot_out[:] = pot
dataout.close()