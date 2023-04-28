from datetime import datetime, timedelta
from numpy import ndarray, array
from netCDF4 import Dataset


# input: SuperDARN gridded ion drifts (grid2 text)
# output: SuperDARN ion drifts (netcdf)


start_date = datetime(year=2015, month=3, day=17)
stop_date = datetime(year=2015, month=3, day=18)
dt = timedelta(days=1)
hemi = 'north'

ntime = 718 * (stop_date - start_date).days

# MAXRECS dimension keeps the record of line-of-sight ion drifts at each 2-minute interval
MAXRECS = 2000

time = ndarray(shape=ntime, dtype=int)
nrec = ndarray(shape=ntime, dtype=int)
gmlat = ndarray(shape=(ntime, MAXRECS))
gmlong = ndarray(shape=(ntime, MAXRECS))
kvect = ndarray(shape=(ntime, MAXRECS))
vlos = ndarray(shape=(ntime, MAXRECS))
vlos_sd = ndarray(shape=(ntime, MAXRECS))

cur_date = start_date
while cur_date < stop_date:
    datain = open(file='{:s}.{:s}.grd'.format(cur_date.strftime('%Y%m%d'), hemi))

    itime = 0
    for hour in range(24):
        for minute in range(1, 60, 2):
# no records at 00:00:00-00:02:00 and 23:58:00-24:00:00, skip
            if (hour == 0 and minute == 1) or (hour == 23 and minute == 59):
                continue

            t = hour*60 + minute
            rec = array(datain.readline().split()).astype(int)
            assert (rec[0] == cur_date.year and rec[1] == cur_date.month and rec[2] == cur_date.day and rec[5] == 0 and
                rec[6] == cur_date.year and rec[7] == cur_date.month and rec[8] == cur_date.day and rec[11] == 0 and
                rec[3]*60 + rec[4] == t - 1 and rec[9]*60 + rec[10] == t + 1)

            cur_time = datetime(year=cur_date.year, month=cur_date.month, day=cur_date.day, hour=hour, minute=minute)
            time[itime] = int((cur_time - start_date).total_seconds() / 60)

# number of stations can be found in the preamble
            datain.readline()
            nstation = int(datain.readline().split()[0])
            for _ in range(3):
                datain.readline()
            for _ in range(nstation):
                datain.readline()

# number of records can be found in the data section
            nrec[itime] = int(datain.readline().split()[0])
            for _ in range(3):
                datain.readline()
            for irec in range(nrec[itime]):
                rec = datain.readline().split()
                gmlong[itime, irec] = float(rec[0])
                gmlat[itime, irec] = float(rec[1])
                kvect[itime, irec] = float(rec[2])
                vlos[itime, irec] = float(rec[6])
                vlos_sd[itime, irec] = float(rec[7])

            itime += 1
    datain.close()
    cur_date += dt

maxnrec = nrec.max()

dataout = Dataset(filename='ion_drift_{:s}.nc'.format(hemi), mode='w')
dataout.createDimension(dimname='time', size=ntime)
dataout.createDimension(dimname='nrec', size=maxnrec)
time_out = dataout.createVariable(varname='time', datatype='i4', dimensions='time')
nrec_out = dataout.createVariable(varname='nrec', datatype='i4', dimensions='time')
gmlat_out = dataout.createVariable(varname='gmlat', datatype='f4', dimensions=('time', 'nrec'))
gmlong_out = dataout.createVariable(varname='gmlong', datatype='f4', dimensions=('time', 'nrec'))
kvect_out = dataout.createVariable(varname='kvect', datatype='f4', dimensions=('time', 'nrec'))
vlos_out = dataout.createVariable(varname='vlos', datatype='f4', dimensions=('time', 'nrec'))
vlos_sd_out = dataout.createVariable(varname='vlos_sd', datatype='f4', dimensions=('time', 'nrec'))
time_out.units = 'minutes since {:s} 00:00:00'.format(start_date.strftime('%Y-%m-%d'))
time_out[:] = time
nrec_out[:] = nrec
gmlat_out[:] = gmlat[:, 0: maxnrec]
gmlong_out[:] = gmlong[:, 0: maxnrec]
kvect_out[:] = kvect[:, 0: maxnrec]
vlos_out[:] = vlos[:, 0: maxnrec]
vlos_sd_out[:] = vlos_sd[:, 0: maxnrec]
dataout.close()