from netCDF4 import Dataset
from numpy import ndarray, interp, deg2rad, pi, cos, sin, meshgrid, linspace
from scipy.interpolate import griddata


hemi = 'north'

data = Dataset(filename='../superdarn/ion_drift_mlt_{:s}.nc'.format(hemi))
time = data['time']
timeunits = time.units
time = time[:].filled()
nrec = data['nrec'][:].filled()
gmlat = data['mlat'][:].filled()
gmlt = data['mlt'][:].filled()
Ekvect = data['Ekvect'][:].filled()
Elos = data['E_median'][:].filled()
Elos_sd = data['E_sd'][:].filled()
data.close()

if hemi == 'south':
    gmlat = -gmlat
    Ekvect -= 180

bg = Dataset(filename='../weimer/weimer_ext.nc')
time_bg = bg['time'][:].filled()
mlat_bg = bg['mlat'][:].filled()
mlt_bg = bg['mlt'][:].filled()
potential_bg = bg['potential'][:].filled()
bg.close()

ntime = len(time)
nmlat_bg = len(mlat_bg)
nmlt_bg = len(mlt_bg)

pot = ndarray(shape=(ntime, nmlat_bg, nmlt_bg))
for imlat in range(nmlat_bg):
    for imlt in range(nmlt_bg):
        pot[:, imlat, imlt] = interp(x=time, xp=time_bg, fp=potential_bg[:, imlat, imlt])

a = 6.371e6
theta = deg2rad(mlat_bg)
phi = mlt_bg * pi/12

Ex = ndarray(shape=(ntime, nmlat_bg, nmlt_bg))
Ex[:, :, 0] = pot[:, :, 1] - pot[:, :, nmlt_bg-1]
Ex[:, :, 1: nmlt_bg-1] = pot[:, :, 2: nmlt_bg] - pot[:, :, 0: nmlt_bg-2]
Ex[:, :, nmlt_bg-1] = pot[:, :, 0] - pot[:, :, nmlt_bg-2]
for imlat in range(nmlat_bg):
    Ex[:, imlat, :] /= cos(theta[imlat])
Ex /= -a * (phi[2] - phi[0])

Ey = -(pot[:, 2: nmlat_bg, :] - pot[:, 0: nmlat_bg-2, :]) / (a * (theta[2] - theta[0]))

t, r = meshgrid(phi, pi/2-theta)
x = r * cos(t)
y = r * sin(t)
coor = (x.flatten(), y.flatten())

t, r = meshgrid(phi, pi/2-theta[1: nmlat_bg-1])
x = r * cos(t)
y = r * sin(t)
coor_Ey = (x.flatten(), y.flatten())

mlt_lb, mlat_lb = meshgrid(linspace(start=0, stop=24, endpoint=False, num=24), [35, 40])
mlt_lb = mlt_lb.flatten()
mlat_lb = mlat_lb.flatten()
cnt_lb = len(mlt_lb)
r = pi/2 - deg2rad(mlat_lb)
t = mlt_lb * pi/12
xi = (r*cos(t), r*sin(t))

nrec_lb = cnt_lb * 2
gmlat_lb = ndarray(shape=(ntime, nrec_lb))
gmlt_lb = ndarray(shape=(ntime, nrec_lb))
Ekvect_lb = ndarray(shape=(ntime, nrec_lb))
Elos_lb = ndarray(shape=(ntime, nrec_lb))
Elos_sd_lb = ndarray(shape=(ntime, nrec_lb))
for itime in range(ntime):
    Ex_lb = griddata(points=coor, values=Ex[itime, :, :].flatten(), xi=xi)
    Ey_lb = griddata(points=coor_Ey, values=Ey[itime, :, :].flatten(), xi=xi)

    gmlat_lb[itime, 0: cnt_lb] = mlat_lb
    gmlat_lb[itime, cnt_lb: nrec_lb] = mlat_lb
    gmlt_lb[itime, 0: cnt_lb] = mlt_lb
    gmlt_lb[itime, cnt_lb: nrec_lb] = mlt_lb
    Ekvect_lb[itime, 0: cnt_lb] = 90
    Ekvect_lb[itime, cnt_lb: nrec_lb] = 0
    Elos_lb[itime, 0: cnt_lb] = Ex_lb
    Elos_lb[itime, cnt_lb: nrec_lb] = Ey_lb
    Elos_sd_lb[itime, :] = 0.1

maxnrec_all = nrec.max() + nrec_lb
nrec_all = ndarray(shape=ntime, dtype=int)
gmlat_all = ndarray(shape=(ntime, maxnrec_all))
gmlt_all = ndarray(shape=(ntime, maxnrec_all))
Ekvect_all = ndarray(shape=(ntime, maxnrec_all))
Elos_all = ndarray(shape=(ntime, maxnrec_all))
Elos_sd_all = ndarray(shape=(ntime, maxnrec_all))
for itime in range(ntime):
    n = nrec[itime]

    gmlat_all[itime, 0: n] = gmlat[itime, 0: n]
    gmlt_all[itime, 0: n] = gmlt[itime, 0: n]
    Ekvect_all[itime, 0: n] = Ekvect[itime, 0: n]
    Elos_all[itime, 0: n] = Elos[itime, 0: n]
    Elos_sd_all[itime, 0: n] = Elos_sd[itime, 0: n]

    gmlat_all[itime, n: n+nrec_lb] = gmlat_lb[itime, :]
    gmlt_all[itime, n: n+nrec_lb] = gmlt_lb[itime, :]
    Ekvect_all[itime, n: n+nrec_lb] = Ekvect_lb[itime, :]
    Elos_all[itime, n: n+nrec_lb] = Elos_lb[itime, :]
    Elos_sd_all[itime, n: n+nrec_lb] = Elos_sd_lb[itime, :]

    nrec_all[itime] = n + nrec_lb

ex = ndarray(shape=(ntime, maxnrec_all))
ey = ndarray(shape=(ntime, maxnrec_all))
for itime in range(ntime):
    n = nrec_all[itime]
    r = pi/2 - deg2rad(gmlat_all[itime, 0: n])
    t = gmlt_all[itime, 0: n] * pi/12
    xi = (r*cos(t), r*sin(t))

    ex[itime, 0: n] = griddata(points=coor, values=Ex[itime, :, :].flatten(), xi=xi)
    ey[itime, 0: n] = griddata(points=coor_Ey, values=Ey[itime, :, :].flatten(), xi=xi)

r = pi/2 - deg2rad(gmlat_all)
t = gmlt_all * pi/12
azim = deg2rad(Ekvect_all)

data_in = Dataset(filename=hemi+'_in.nc', mode='w')
data_in.createDimension(dimname='time', size=ntime)
data_in.createDimension(dimname='nrec', size=maxnrec_all)
time_out = data_in.createVariable(varname='time', datatype='i', dimensions='time')
nrec_out = data_in.createVariable(varname='nrec', datatype='i', dimensions='time')
x_out = data_in.createVariable(varname='x', datatype='f', dimensions=('time', 'nrec'))
y_out = data_in.createVariable(varname='y', datatype='f', dimensions=('time', 'nrec'))
azim_out = data_in.createVariable(varname='azim', datatype='f', dimensions=('time', 'nrec'))
Elos_out = data_in.createVariable(varname='Elos', datatype='f', dimensions=('time', 'nrec'))
Elos_sd_out = data_in.createVariable(varname='Elos_sd', datatype='f', dimensions=('time', 'nrec'))
Elos_prior_out = data_in.createVariable(varname='Elos_prior', datatype='f', dimensions=('time', 'nrec'))
time_out.units = timeunits
time_out[:] = time
nrec_out[:] = nrec_all
x_out[:] = r * cos(t)
y_out[:] = r * sin(t)
azim_out[:] = t - azim + pi
Elos_out[:] = Elos_all
Elos_sd_out[:] = Elos_sd_all
Elos_prior_out[:] = ex*sin(azim) + ey*cos(azim)
data_in.close()

nxy = 161
xy = linspace(start=-pi*2/9, stop=pi*2/9, num=nxy)
x, y = meshgrid(xy, xy)
xi = (x.flatten(), y.flatten())
prior = ndarray(shape=(ntime, nxy, nxy))
for itime in range(ntime):
    prior[itime, :, :] = griddata(points=coor, values=pot[itime, :, :].flatten(), xi=xi).reshape((nxy, nxy))

data_out = Dataset(filename=hemi+'_out_grid.nc', mode='w')
data_out.createDimension(dimname='time', size=ntime)
data_out.createDimension(dimname='x', size=nxy)
data_out.createDimension(dimname='y', size=nxy)
time_out = data_out.createVariable(varname='time', datatype='i', dimensions='time')
x_out = data_out.createVariable(varname='x', datatype='f', dimensions='x')
y_out = data_out.createVariable(varname='y', datatype='f', dimensions='y')
prior_out = data_out.createVariable(varname='prior', datatype='f', dimensions=('time', 'y', 'x'))
time_out.units = timeunits
time_out[:] = time
x_out[:] = xy
y_out[:] = xy
prior_out[:] = prior
data_out.close()
