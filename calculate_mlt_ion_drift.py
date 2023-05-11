from netCDF4 import Dataset, num2date
from datetime import datetime, timedelta
from numpy import ndarray
from aacgmv2 import convert_mlt


# calculate magnetic local time for given ut and magnetic locations
# input: SuperDARN ion drifts
# output: SuperDARN ion drifts with magnetic local time


hemi = 'north'

datain = Dataset(filename='ion_drift_{:s}.nc'.format(hemi))
time = datain['time']
timeunits = time.units
time = time[:].filled()
nrec = datain['nrec'][:].filled()
gmlat = datain['gmlat'][:].filled()
gmlong = datain['gmlong'][:].filled()
kvect = datain['kvect'][:].filled()
vlos = datain['vlos'][:].filled()
vlos_sd = datain['vlos_sd'][:].filled()
datain.close()

ntime = len(time)
maxnrec = nrec.max()

start_date = datetime.strptime(num2date(times=0, units=timeunits).strftime('%Y-%m-%d'), '%Y-%m-%d')

# aacgmv2 mlt
gmlt = ndarray(shape=(ntime, maxnrec))
for itime in range(ntime):
    n = nrec[itime]
    gmlt[itime, 0: n] = convert_mlt(arr=gmlong[itime, 0: n], dtime=start_date+timedelta(minutes=time[itime].item()))

dataout = Dataset(filename='ion_drift_mlt_{:s}.nc'.format(hemi), mode='w')
dataout.createDimension(dimname='time', size=ntime)
dataout.createDimension(dimname='nrec', size=maxnrec)
time_out = dataout.createVariable(varname='time', datatype='i4', dimensions='time')
nrec_out = dataout.createVariable(varname='nrec', datatype='i4', dimensions='time')
gmlat_out = dataout.createVariable(varname='gmlat', datatype='f4', dimensions=('time', 'nrec'))
gmlong_out = dataout.createVariable(varname='gmlong', datatype='f4', dimensions=('time', 'nrec'))
gmlt_out = dataout.createVariable(varname='gmlt', datatype='f4', dimensions=('time', 'nrec'))
kvect_out = dataout.createVariable(varname='kvect', datatype='f4', dimensions=('time', 'nrec'))
vlos_out = dataout.createVariable(varname='vlos', datatype='f4', dimensions=('time', 'nrec'))
vlos_sd_out = dataout.createVariable(varname='vlos_sd', datatype='f4', dimensions=('time', 'nrec'))
time_out.units = timeunits
time_out[:] = time
nrec_out[:] = nrec
gmlat_out[:] = gmlat
gmlong_out[:] = gmlong
gmlt_out[:] = gmlt
kvect_out[:] = kvect
vlos_out[:] = vlos
vlos_sd_out[:] = vlos_sd
dataout.close()