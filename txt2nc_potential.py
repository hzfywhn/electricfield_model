from numpy import loadtxt, array, ndarray
from datetime import datetime
from fortranformat import FortranRecordReader
from netCDF4 import Dataset


# input: Weimer potential (text) from gridpotentials.sav
# output: Weimer potential (netcdf)


Year = 2015
Month = 3
Day = 17
start_time = datetime(year=Year, month=Month, day=Day)

# the same format as used for the Weimer model
imfs = loadtxt(fname='imfs.txt')
year = imfs[:, 4].astype(int)
month = imfs[:, 5].astype(int)
day = imfs[:, 6].astype(int)
hour = imfs[:, 7].astype(int)
ntime = len(year)

# MLAT/MLT grids for Weimer outputs
grid = open(file='potentialgrid.txt')
grid.readline()
mlt = array(grid.readline().split()).astype(float)
mlat = array(grid.readline().split()).astype(float)
grid.close()
nmlt = len(mlt)
nmlat = len(mlat)
np = nmlt * nmlat

# Weimer output has a fixed amount of values per line (13)
nfull = 13
nsub = np % nfull
full_line = FortranRecordReader('({:d}F8.3)'.format(nfull))
sub_line = FortranRecordReader('({:d}F8.3)'.format(nsub))

# read the Weimer potential (from gridpotentials.sav) and put into the correct 2d array position
potential = ndarray(shape=(ntime, nmlat, nmlt))
pot = ndarray(shape=np)
datain = open(file='potentials.txt')
itime = 0
i = 0
for line in datain:
    if i+nfull < np:
# current line is part of the model output
        pot[i: i+nfull] = full_line.read(line)
        i += nfull
    else:
# current line marks the end of the model output, reformat to mlat/mlt grids
        pot[i: i+nsub] = sub_line.read(line)
        potential[itime, :, :] = pot.reshape((nmlat, nmlt))
        itime += 1
        i = 0
datain.close()

dataout = Dataset(filename='weimer.nc', mode='w')
dataout.createDimension(dimname='time', size=ntime)
dataout.createDimension(dimname='mlat', size=nmlat)
dataout.createDimension(dimname='mlt', size=nmlt)
time_out = dataout.createVariable(varname='time', datatype='i4', dimensions='time')
mlat_out = dataout.createVariable(varname='mlat', datatype='f4', dimensions='mlat')
mlt_out = dataout.createVariable(varname='mlt', datatype='f4', dimensions='mlt')
potential_out = dataout.createVariable(varname='potential', datatype='f4', dimensions=('time', 'mlat', 'mlt'))
time_out.units = 'minutes since {:4d}-{:02d}-{:02d} 00:00:00'.format(Year, Month, Day)
for itime in range(ntime):
    time_out[itime] = int((datetime(year=year[itime], month=month[itime], day=day[itime], hour=hour[itime]) - start_time).total_seconds() / 60)
mlat_out[:] = mlat
mlt_out[:] = mlt
potential_out[:] = potential * 1e3
dataout.close()