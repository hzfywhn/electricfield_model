from netCDF4 import Dataset, num2date
from datetime import datetime, timedelta
from numpy import ndarray
from aacgmv2 import convert_mlt, convert_latlon
from pyIGRF import igrf_value


# calculate magnetic local time for given ut and magnetic locations, also calculate electric fields using E=-vxB
# input: SuperDARN ion drifts
# output: SuperDARN eletric fields


hemi = 'south'

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

gmlt = ndarray(shape=(ntime, maxnrec))
Elos = ndarray(shape=(ntime, maxnrec))
Elos_sd = ndarray(shape=(ntime, maxnrec))
for itime in range(ntime):
    n = nrec[itime]
    dt = start_date + timedelta(minutes=time[itime].item())

# aacgmv2 mlt
    gmlt[itime, 0: n] = convert_mlt(arr=gmlong[itime, 0: n], dtime=dt)

# find magnetic fields at corresponding locations using IGRF
    for i in range(n):
        glat, glon, galt = convert_latlon(in_lat=gmlat[itime, i], in_lon=gmlt[itime, i], height=100, dtime=dt, method_code='A2G')
        _, _, _, _, _, B, _ = igrf_value(lat=glat, lon=glon, alt=galt, year=dt.year+dt.timetuple().tm_yday/365)
        B = abs(B) * 1e-9
        Elos[itime, i] = vlos[itime, i] * B
        Elos_sd[itime, i] = vlos_sd[itime, i] * B

# in northern hemisphere, magnetic fields are pointing downward
# so electric field line-of-sight directions are on the right
if hemi == 'north':
    Ekvect = kvect + 90

# in southern hemisphere, line-of-sight directions are opposite
if hemi == 'south':
    Ekvect = kvect - 90

dataout = Dataset(filename='E_{:s}.nc'.format(hemi), mode='w')
dataout.createDimension(dimname='time', size=ntime)
dataout.createDimension(dimname='nrec', size=maxnrec)
time_out = dataout.createVariable(varname='time', datatype='i', dimensions='time')
nrec_out = dataout.createVariable(varname='nrec', datatype='i', dimensions='time')
mlat_out = dataout.createVariable(varname='gmlat', datatype='f', dimensions=('time', 'nrec'))
mlt_out = dataout.createVariable(varname='gmlt', datatype='f', dimensions=('time', 'nrec'))
Ekvect_out = dataout.createVariable(varname='Ekvect', datatype='f', dimensions=('time', 'nrec'))
Elos_out = dataout.createVariable(varname='Elos', datatype='f', dimensions=('time', 'nrec'))
Elos_sd_out = dataout.createVariable(varname='Elos_sd', datatype='f', dimensions=('time', 'nrec'))
time_out.units = timeunits
time_out[:] = time
nrec_out[:] = nrec
mlat_out[:] = gmlat
mlt_out[:] = gmlt
Ekvect_out[:] = Ekvect
Elos_out[:] = Elos
Elos_sd_out[:] = Elos_sd
dataout.close()