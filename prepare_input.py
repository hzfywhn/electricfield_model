from numpy import pi, deg2rad, linspace, meshgrid, cos, sin, ndarray, argmin, abs, interp
from netCDF4 import Dataset
from scipy.interpolate import griddata


# prepare input files for modeling
# input: SuperDARN line-of-sight electric fields (from calculate_mlt) and Weimer potential (from extend_potential)
# output: formated input file and output grids


hemi = 'north'
filename_bg = 'weimer_ext.nc'
filename_superdarn = 'E_{:s}.nc'.format(hemi)

# output grid
xylim_out = pi*2/9
dxy_out = deg2rad(0.5)

# internal grid size and resolution used for data preparation
xylim = pi/3
dxy = deg2rad(0.2)

# create a high-resolution internal grid for interpolation
nxy = int(xylim*2/dxy) + 1
xy = linspace(start=-xylim, stop=xylim, num=nxy)
x2d, y2d = meshgrid(xy, xy)
xi = (x2d.flatten(), y2d.flatten())

# Weimer model has been extended to 30 MLAT using an exponential decay function
datain = Dataset(filename=filename_bg)
time_bg = datain['time'][:].filled()
mlat = datain['mlat'][:].filled()
mlt = datain['mlt'][:].filled()
potential = datain['potential'][:].filled()
datain.close()

# used for coordinate transform, from magnetic coordinates to model coordinates
t2d, r2d = meshgrid(mlt*pi/12, pi/2-deg2rad(mlat))

# model coordinate system is a normalized flattened spherical coordinate centered at the north pole
x2d = r2d * cos(t2d)
y2d = r2d * sin(t2d)
points = (x2d.flatten(), y2d.flatten())

# reformat Weimer potential at a very fine grid
ntime_bg = len(time_bg)
pot = ndarray(shape=(ntime_bg, nxy, nxy))
for itime in range(ntime_bg):
    pot[itime, :, :] = griddata(points=points, values=potential[itime, :, :].flatten(), xi=xi).reshape((nxy, nxy))

# calculate Weimer electric fields on the fine grid
a = 6.371e6
Ex = ndarray(shape=(ntime_bg, nxy, nxy))
Ex[:, :, 1: nxy-1] = -(pot[:, :, 2: nxy] - pot[:, :, 0: nxy-2]) / (2*dxy*a)
Ey = ndarray(shape=(ntime_bg, nxy, nxy))
Ey[:, 1: nxy-1, :] = -(pot[:, 2: nxy, :] - pot[:, 0: nxy-2, :]) / (2*dxy*a)

# note that SuperDARN measurements have been converted to electric fields at this step
datain = Dataset(filename=filename_superdarn)
time = datain['time']
timeunits = time.units
time = time[:].filled()
nrec = datain['nrec'][:].filled()
gmlat = datain['gmlat'][:].filled()
gmlt = datain['gmlt'][:].filled()
Ekvect = datain['Ekvect'][:].filled()
Elos = datain['Elos'][:].filled()
Elos_sd = datain['Elos_sd'][:].filled()
datain.close()

# the southern hemisphere is modeled by looking downward from the north pole (across the globe)
if hemi == 'south':
    gmlat = -gmlat
    Ekvect -= 180

r = pi/2 - deg2rad(gmlat)
t = gmlt * pi/12
x = r * cos(t)
y = r * sin(t)

# the azimuth is counted as the angle between x axis in model coordinate
# SuperDARN measurements count azimuth from local magnetic north
# the transform is valid in the northern hemisphere
# but as gmlat and Ekvect have been switched to the northern hemisphere, this transform is valid for both hemispheres
azim = t - deg2rad(Ekvect) + pi

# add fake observations near the lower latitude boundary to prevent the model from explosion
# 24 (mlt) x 2 (mlat) x 2 (direction) fake observations are typically enough to do
t2d, r2d = meshgrid(linspace(start=0, stop=24, endpoint=False, num=24)*pi/12, pi/2-deg2rad([35, 40]))
x2d_lb = (r2d * cos(t2d)).flatten()
y2d_lb = (r2d * sin(t2d)).flatten()
cnt_lb = len(x2d_lb)
nrec_lb = cnt_lb * 2

x_lb = ndarray(shape=nrec_lb)
x_lb[0: cnt_lb] = x2d_lb
x_lb[cnt_lb: nrec_lb] = x2d_lb
y_lb = ndarray(shape=nrec_lb)
y_lb = ndarray(shape=nrec_lb)
y_lb[0: cnt_lb] = y2d_lb
y_lb[cnt_lb: nrec_lb] = y2d_lb

# the line-of-sight directions are chosen to be toward right (Ex) and upward (Ey)
azim_lb = ndarray(shape=nrec_lb)
azim_lb[0: cnt_lb] = 0
azim_lb[cnt_lb: nrec_lb] = pi/2

# find the corresponding line-of-sight electric fields at fake observation locations
ntime = len(time)
Elos_lb = ndarray(shape=(ntime, nrec_lb))
for i in range(cnt_lb):
    ix = argmin(abs(xy - x2d_lb[i]))
    iy = argmin(abs(xy - y2d_lb[i]))
    Elos_lb[:, i] = interp(x=time, xp=time_bg, fp=Ex[:, iy, ix])
    Elos_lb[:, i+cnt_lb] = interp(x=time, xp=time_bg, fp=Ey[:, iy, ix])

Elos_sd_lb = 0.1

# combine real observations with fake observations to form a complete observation set
maxnrec_all = nrec.max() + nrec_lb
nrec_all = ndarray(shape=ntime, dtype=int)
x_all = ndarray(shape=(ntime, maxnrec_all))
y_all = ndarray(shape=(ntime, maxnrec_all))
azim_all = ndarray(shape=(ntime, maxnrec_all))
Elos_all = ndarray(shape=(ntime, maxnrec_all))
Elos_sd_all = ndarray(shape=(ntime, maxnrec_all))
for itime in range(ntime):
    n = nrec[itime]

    x_all[itime, 0: n] = x[itime, 0: n]
    y_all[itime, 0: n] = y[itime, 0: n]
    azim_all[itime, 0: n] = azim[itime, 0: n]
    Elos_all[itime, 0: n] = Elos[itime, 0: n]
    Elos_sd_all[itime, 0: n] = Elos_sd[itime, 0: n]

    x_all[itime, n: n+nrec_lb] = x_lb
    y_all[itime, n: n+nrec_lb] = y_lb
    azim_all[itime, n: n+nrec_lb] = azim_lb
    Elos_all[itime, n: n+nrec_lb] = Elos_lb[itime, :]
    Elos_sd_all[itime, n: n+nrec_lb] = Elos_sd_lb

    nrec_all[itime] = n + nrec_lb

# find the prior estimate of electric fields at corresponding locations
Elos_prior = ndarray(shape=(ntime, maxnrec_all))
for itime in range(ntime):
    for i in range(nrec_all[itime]):
        ix = argmin(abs(xy - x_all[itime, i]))
        iy = argmin(abs(xy - y_all[itime, i]))
        ex = interp(x=time[itime], xp=time_bg, fp=Ex[:, iy, ix])
        ey = interp(x=time[itime], xp=time_bg, fp=Ey[:, iy, ix])
        Elos_prior[itime, i] = ex*sin(azim_all[itime, i]) + ey*cos(azim_all[itime, i])

# create output grids at a desired resolution
nxy_out = int(xylim_out*2/dxy_out) + 1
xy_out = linspace(start=-xylim_out, stop=xylim_out, num=nxy_out)

# note that the prior is potential (not electric field)
prior = ndarray(shape=(ntime, nxy_out, nxy_out))
for itime in range(ntime):
    for i in range(nxy_out):
        ix = argmin(abs(xy - xy_out[i]))
        for j in range(nxy_out):
            iy = argmin(abs(xy - xy_out[j]))
            prior[itime, j, i] = interp(x=time[itime], xp=time_bg, fp=pot[:, iy, ix])

dataout = Dataset(filename=hemi+'_in.nc', mode='w')
dataout.createDimension(dimname='time', size=ntime)
dataout.createDimension(dimname='nrec', size=maxnrec_all)
time_out = dataout.createVariable(varname='time', datatype='i4', dimensions='time')
nrec_out = dataout.createVariable(varname='nrec', datatype='i4', dimensions='time')
x_out = dataout.createVariable(varname='x', datatype='f4', dimensions=('time', 'nrec'))
y_out = dataout.createVariable(varname='y', datatype='f4', dimensions=('time', 'nrec'))
azim_out = dataout.createVariable(varname='azim', datatype='f4', dimensions=('time', 'nrec'))
Elos_out = dataout.createVariable(varname='Elos', datatype='f4', dimensions=('time', 'nrec'))
Elos_sd_out = dataout.createVariable(varname='Elos_sd', datatype='f4', dimensions=('time', 'nrec'))
Elos_prior_out = dataout.createVariable(varname='Elos_prior', datatype='f4', dimensions=('time', 'nrec'))
time_out.units = timeunits
time_out[:] = time
nrec_out[:] = nrec_all
x_out[:] = x_all
y_out[:] = y_all
azim_out[:] = azim_all
Elos_out[:] = Elos_all
Elos_sd_out[:] = Elos_sd_all
Elos_prior_out[:] = Elos_prior
dataout.close()

dataout = Dataset(filename=hemi+'_out_grid.nc', mode='w')
dataout.createDimension(dimname='time', size=ntime)
dataout.createDimension(dimname='x', size=nxy_out)
dataout.createDimension(dimname='y', size=nxy_out)
time_out = dataout.createVariable(varname='time', datatype='i4', dimensions='time')
x_out = dataout.createVariable(varname='x', datatype='f4', dimensions='x')
y_out = dataout.createVariable(varname='y', datatype='f4', dimensions='y')
prior_out = dataout.createVariable(varname='prior', datatype='f4', dimensions=('time', 'y', 'x'))
time_out.units = timeunits
time_out[:] = time
x_out[:] = xy_out
y_out[:] = xy_out
prior_out[:] = prior
dataout.close()