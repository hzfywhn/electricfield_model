from netCDF4 import Dataset, num2date
from datetime import datetime, timedelta
from numpy import pi, deg2rad, arange, cos, sin
from matplotlib import use
from matplotlib.pyplot import figure


# draw quiver plots of ion drifts in magnetic latitude/local time coordinates
# for easy comparison to http://vt.superdarn.org/tiki-index.php?page=Map+Velocity+Tool


hemi = 'north'

data = Dataset(filename='ion_drift_mlt_{:s}.nc'.format(hemi))
time = data['time']
timeunits = time.units
time = time[:].filled()
nrec = data['nrec'][:].filled()
gmlat = data['gmlat'][:].filled()
gmlt = data['gmlt'][:].filled()
kvect = data['kvect'][:].filled()
vlos = data['vlos'][:].filled()
data.close()

ntime = len(time)
start_date = datetime.strptime(num2date(times=0, units=timeunits).strftime('%Y-%m-%d'), '%Y-%m-%d')

# pi/2 is used to rotate the map to set zero local time at bottom
t = gmlt*pi/12 - pi/2

# note the calculation of azimuth is different for two hemispheres (not just a minus sign)
if hemi == 'north':
    azim = pi - deg2rad(kvect) + t
    r = pi/2 - deg2rad(gmlat)
    yticklabels = [90, 70, 50]
    yticks = pi/2 - deg2rad(yticklabels)
if hemi == 'south':
    azim = deg2rad(kvect) + t
    r = pi/2 + deg2rad(gmlat)
    yticklabels = [-90, -70, -50]
    yticks = pi/2 + deg2rad(yticklabels)

axislim = (-yticks[-1], yticks[-1])
xticks = arange(start=0, stop=pi*2, step=pi/2)
xticklabels = ['{:02d}'.format(ihr) for ihr in arange(start=0, stop=24, step=6)]

use('Agg')
fig = figure()
for itime in range(ntime):
    n = nrec[itime]
    dt = start_date + timedelta(minutes=time[itime].item())

    x = r[itime, 0: n] * cos(t[itime, 0: n])
    y = r[itime, 0: n] * sin(t[itime, 0: n])
    u = vlos[itime, 0: n] * cos(azim[itime, 0: n])
    v = vlos[itime, 0: n] * sin(azim[itime, 0: n])

# initialize an empty canvas
    axes = fig.add_subplot(polar=True)
    axes.set_theta_zero_location(loc='S')
    axes.set_xticks(ticks=xticks)
    axes.set_xticklabels(labels=xticklabels)
    axes.set_ylim(bottom=yticks[0], top=yticks[-1])
    axes.set_yticks(ticks=yticks)
    axes.set_yticklabels(labels=yticklabels)
    axes.set_title(label='{:s} {:s}'.format(hemi, dt.strftime("%Y/%m/%d %H:%M")))

# add arrows to the plot
    fig.add_axes(rect=axes.get_position(), aspect='equal', frame_on=False, xlim=axislim, xticks=[], ylim=axislim, yticks=[]).quiver(x, y, u, v, scale=1e4, width=2e-3)

    fig.savefig(fname='ion_drift_{:s}/{:s}.png'.format(hemi, dt.strftime("%Y%m%d%H%M")), bbox_inches='tight')
    fig.clear()