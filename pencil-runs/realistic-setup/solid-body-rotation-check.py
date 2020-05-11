import pencil
import numpy
from matplotlib import pyplot
from matplotlib import cm
import sys
import os

ivar = -1
pvar = "pvar.dat"
if(len(sys.argv) > 1):
    ivar = int(sys.argv[1])
    pvar = "PVAR" + sys.argv[1]

data        = pencil.read_var(ivar=ivar,trimall=True,quiet=True)
xgrid   = data.x ; nx = len(xgrid) ; dxgrid  = data.dx ; x0 = xgrid[0] ; x1 = xgrid[-1]
ygrid   = data.y ; ny = len(ygrid) ; dygrid  = data.dy ; y0 = ygrid[0] ; y1 = ygrid[-1]
zgrid   = data.z ; nz = len(zgrid) ; dzgrid  = data.dz ; z0 = zgrid[0] ; z1 = zgrid[-1]

z_ref = 0.0
if(len(sys.argv) > 2):
    z_ref = float(sys.argv[2])
iz_ref = numpy.where(numpy.abs(z_ref - zgrid) == numpy.abs(z_ref - zgrid).min())[0][0]
ux = numpy.transpose(data.ux[iz_ref,:,:])
uy = numpy.transpose(data.uy[iz_ref,:,:])
uz = numpy.transpose(data.uz[iz_ref,:,:])

print(iz_ref)

omega  = 1.05
ux_rot = numpy.zeros((nx,ny))
uy_rot = numpy.zeros((nx,ny))
for ix in range(nx):
    for iy in range(ny):
        x_l  = xgrid[ix]                ; y_l    = ygrid[iy]
        ux_l = ux[ix,iy]                ; uy_l   = uy[ix,iy]
        r_l  = (x_l**2 + y_l**2)**(0.5) ; phi_l  = numpy.arctan2(y_l,x_l)
        u_l  = omega/r_l
        ux_rot[ix,iy] = u_l*numpy.sin(phi_l) ; uy_rot[ix,iy] = u_l*numpy.cos(phi_l)

fig, ((ax1,ax2,ax3)) = pyplot.subplots(1,3)
fig.set_size_inches(19.2,10.8)
fig.set_dpi(100)

contour1 = ax1.pcolormesh(xgrid,ygrid,numpy.abs(ux_rot - ux))
ax1.set_xlim([x0,x1])
ax1.set_ylim([y0,y1])
ax1.set_aspect("equal")
fig.colorbar(contour1,ax=ax1)
contour2 = ax2.pcolormesh(xgrid,ygrid,numpy.abs(uy_rot - uy))
ax2.set_xlim([x0,x1])
ax2.set_ylim([y0,y1])
ax2.set_aspect("equal")
fig.colorbar(contour2,ax=ax2)
contour3 = ax3.pcolormesh(xgrid,ygrid,numpy.abs(uz))
ax3.set_xlim([x0,x1])
ax3.set_ylim([y0,y1])
ax3.set_aspect("equal")
fig.colorbar(contour3,ax=ax3)

ax1.set_xlabel(r"x",fontsize=24)
ax1.set_ylabel(r"y",fontsize=24)
ax2.set_xlabel(r"x",fontsize=24)
ax2.set_ylabel(r"y",fontsize=24)
ax3.set_xlabel(r"x",fontsize=24)
ax3.set_ylabel(r"y",fontsize=24)

for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(19)
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(19)
for tick in ax2.xaxis.get_major_ticks():
    tick.label.set_fontsize(19)
for tick in ax2.yaxis.get_major_ticks():
    tick.label.set_fontsize(19)
for tick in ax3.xaxis.get_major_ticks():
    tick.label.set_fontsize(19)
for tick in ax3.yaxis.get_major_ticks():
    tick.label.set_fontsize(19)

pyplot.show()
