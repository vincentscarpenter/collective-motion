import pencil
import numpy
from matplotlib import pyplot
from matplotlib import patches
import sys

parameters = pencil.read_param()
x0 = parameters.xyz0[0]
y0 = parameters.xyz0[1]
z0 = parameters.xyz0[2]
x1 = parameters.xyz1[0]
y1 = parameters.xyz1[1]
z1 = parameters.xyz1[2]

xbox_length =      x1 - x0 
xbox_center = 0.5*(x1 + x0)
ybox_length =      y1 - y0 
ybox_center = 0.5*(y1 + y0)
zbox_length =      z1 - z0 
zbox_center = 0.5*(z1 + z0)
var          = "pvar.dat"
ivar         = -1
if len(sys.argv) == 2:
    ivar = int(sys.argv[1])
    var  = "PVAR" + str(ivar)
elif len(sys.argv) == 5:
    ivar = int(sys.argv[1])
    var  = "PVAR" + str(ivar)
    xbox_length = float(sys.argv[2])
    ybox_length = float(sys.argv[3])
    zbox_length = float(sys.argv[4])
    x_out       = (xbox_center - 0.5*xbox_length) < x0 or (xbox_center + 0.5*xbox_length) > x1
    y_out       = (ybox_center - 0.5*ybox_length) < y0 or (ybox_center + 0.5*ybox_length) > y1
    z_out       = (zbox_center - 0.5*zbox_length) < z0 or (zbox_center + 0.5*zbox_length) > z1
    if(x_out or y_out or z_out):
        print("Box extends beyond the domain. Please pick a box that fits within the domain.")
        sys.exit("Exiting...")
elif len(sys.argv) == 8:
    ivar = int(sys.argv[1])
    var  = "PVAR" + str(ivar)
    if(sys.argv[2][0] == 'm'):
        xbox_center = -float(sys.argv[2][1:])
    else:
        xbox_center = float(sys.argv[2])
    xbox_length = float(sys.argv[3])
    if(sys.argv[4][0] == 'm'):
        ybox_center = -float(sys.argv[4][1:])
    else:
        ybox_center = float(sys.argv[4])
    ybox_length = float(sys.argv[5])
    if(sys.argv[6][0] == 'm'):
        zbox_center = -float(sys.argv[6][1:])
    else:
        zbox_center = float(sys.argv[6])
    zbox_length = float(sys.argv[7])
    x_out       = (xbox_center - 0.5*xbox_length) < x0 or (xbox_center + 0.5*xbox_length) > x1
    y_out       = (ybox_center - 0.5*ybox_length) < y0 or (ybox_center + 0.5*ybox_length) > y1
    z_out       = (zbox_center - 0.5*zbox_length) < z0 or (zbox_center + 0.5*zbox_length) > z1
    if(x_out or y_out or z_out):
        print("Box extends beyond the domain. Please pick a box that fits within the domain.")
        sys.exit("Exiting...")

pdata = pencil.read_pvar(varfile=var)
data  = pencil.read_var(ivar=ivar,trimall=True)

xp     = pdata.xp
yp     = pdata.yp
zp     = pdata.zp
vpx    = pdata.vpx
vpy    = pdata.vpy
vpz    = pdata.vpz
xgrid  = data.x
ygrid  = data.y
zgrid  = data.z
dxgrid = data.dx
dygrid = data.dy
dzgrid = data.dz
vz     = numpy.transpose(data.uz)

xbox   = numpy.where(numpy.abs(xp - xbox_center) <= xbox_length/2.0)
ybox   = numpy.where(numpy.abs(yp - ybox_center) <= ybox_length/2.0)
zbox   = numpy.where(numpy.abs(zp - zbox_center) <= zbox_length/2.0)
xybox  = numpy.intersect1d(xbox,ybox)
xyzbox = numpy.intersect1d(xybox,zbox)
xpar_ref = xp[xyzbox]
ypar_ref = yp[xyzbox]
zpar_ref = zp[xyzbox]
vpz_ref  = vpz[xyzbox]
npar_box = len(xyzbox)

xbox_gas = numpy.where(numpy.abs(xgrid - xbox_center) <= xbox_length/2.0)
ybox_gas = numpy.where(numpy.abs(ygrid - ybox_center) <= ybox_length/2.0)
zbox_gas = numpy.where(numpy.abs(zgrid - zbox_center) <= zbox_length/2.0)

gas_reference_velocity = numpy.zeros(npar_box)

#gas_vz = numpy.zeros(npar_box)
#particles_excluded = False
#npar_excluded      = 0
#for ipar in range(npar_box):
#    xpar     = xpar_ref[ipar]
#    ypar     = ypar_ref[ipar]
#    zpar     = zpar_ref[ipar]
#    ixp      = numpy.where(numpy.abs(xpar-xgrid) == numpy.min(numpy.abs(xpar-xgrid)))[0][0]
#    iyp      = numpy.where(numpy.abs(ypar-ygrid) == numpy.min(numpy.abs(ypar-ygrid)))[0][0]
#    izp      = numpy.where(numpy.abs(zpar-zgrid) == numpy.min(numpy.abs(zpar-zgrid)))[0][0]
#    if(ixp == x0 or ixp == x1 or iyp == y0 or iyp == y1 or izp == z0 or izp == z1):
#        particles_excluded = True
#        numpy.delete(gas_vz,ipar)
#        numpy.delete(vpz_ref,ipar)
#        numpy.delete(xpar_ref,ipar)
#        numpy.delete(ypar_ref,ipar)
#        numpy.delete(zpar_ref,ipar)
#        print("WARNING: a particle in the box is too close to the boundary and has been excluded.")
#        npar_excluded+=1
#        continue
#
# From the Fortran shit
#
#    dxp0     = (xpar-xgrid[ixp])/dxgrid
#    dyp0     = (ypar-ygrid[iyp])/dygrid
#    dzp0     = (zpar-zgrid[izp])/dzgrid
#    fac_x_m1 = 0.5*(0.5-dxp0)**2
#    fac_x_00 = 0.75-dxp0**2
#    fac_x_p1 = 0.5*(0.5+dxp0)**2
#    fac_y_m1 = 0.5*(0.5-dyp0)**2
#    fac_y_00 = 0.75-dyp0**2
#    fac_y_p1 = 0.5*(0.5+dyp0)**2
#    fac_z_m1 = 0.5*(0.5-dzp0)**2
#    fac_z_00 = 0.75-dzp0**2
#    fac_z_p1 = 0.5*(0.5+dzp0)**2
#    gas_vz[ipar] = fac_x_00*fac_y_00*fac_z_00*vz[ixp,iyp,izp] + fac_x_00*fac_y_00*(vz[ixp,iyp,izp+1]*fac_z_p1 + vz[ixp,iyp,izp-1]*fac_z_m1) + fac_x_00*fac_z_00*(vz[ixp,iyp+1,izp]*fac_y_p1 + vz[ixp,iyp-1,izp]*fac_y_m1) + fac_y_00*fac_z_00*(vz[ixp+1,iyp,izp]*fac_x_p1 + vz[ixp-1,iyp,izp]*fac_x_m1) + fac_x_p1*fac_y_p1*(vz[ixp+1,iyp+1,izp+1]*fac_z_p1 + vz[ixp+1,iyp+1,izp-1]*fac_z_m1) + fac_x_p1*fac_y_m1*(vz[ixp+1,iyp-1,izp+1]*fac_z_p1 + vz[ixp+1,iyp-1,izp-1]*fac_z_m1) + fac_x_m1*fac_y_p1*(vz[ixp-1,iyp+1,izp+1]*fac_z_p1 + vz[ixp-1,iyp+1,izp-1]*fac_z_m1) + fac_x_m1*fac_y_m1*(vz[ixp-1,iyp-1,izp+1]*fac_z_p1 + vz[ixp-1,iyp-1,izp-1]*fac_z_m1) + fac_x_00*fac_y_p1*(vz[ixp,iyp+1,izp+1]*fac_z_p1 + vz[ixp,iyp+1,izp-1]*fac_z_m1) + fac_x_00*fac_y_m1*(vz[ixp,iyp-1,izp+1]*fac_z_p1 + vz[ixp,iyp-1,izp-1]*fac_z_m1) + fac_y_00*fac_z_p1*(vz[ixp+1,iyp,izp+1]*fac_x_p1 + vz[ixp-1,iyp,izp+1]*fac_x_m1) + fac_y_00*fac_z_m1*(vz[ixp+1,iyp,izp-1]*fac_x_p1 + vz[ixp-1,iyp,izp-1]*fac_x_m1) + fac_z_00*fac_x_p1*(vz[ixp+1,iyp+1,izp]*fac_y_p1 + vz[ixp+1,iyp-1,izp]*fac_y_m1) + fac_z_00*fac_x_m1*(vz[ixp-1,iyp+1,izp]*fac_y_p1 + vz[ixp-1,iyp-1,izp]*fac_y_m1)
#
#print(str(npar_excluded) + " excluded.")
gas_vz_avg = numpy.mean(vz[xbox_gas[0][0]:xbox_gas[0][-1],ybox_gas[0][0]:ybox_gas[0][-1],zbox_gas[0][0]:zbox_gas[0][-1]])
print("Mean gas velocity in the box is:")
print(str(gas_vz_avg*2000.) + " mm/s")
vpz_avg = numpy.mean(vpz_ref - gas_reference_velocity)
print("Mean relative velocity of particles in the box is:")
print(str(vpz_avg*2000.) + " mm/s")

closeness = numpy.zeros(len(xyzbox))
for ipar_ref_loop in range(len(xyzbox)):
    xpar_ref_loop = xpar_ref[ipar_ref_loop]
    zpar_ref_loop = zpar_ref[ipar_ref_loop]
    for ipar_loop in range(len(xyzbox)):
        if(ipar_loop == ipar_ref_loop):
            continue
        else:
            xpar_loop = xpar_ref[ipar_loop]
            zpar_loop = zpar_ref[ipar_loop]
            closeness[ipar_ref_loop] += 0.001*((xpar_ref_loop - xpar_loop)**2 + (zpar_ref_loop - zpar_loop)**2)**(-0.5)
closeness_avg = numpy.mean(closeness)
print("Mean closeness in the box is:")
print(closeness_avg)

dtog = ((4./3.)*numpy.pi*(0.0000825)**3*len(xyzbox))*390./(xbox_length*ybox_length*zbox_length)
print("Dust to gas ratio in the box is:")
print(dtog)

fig, ((ax1,ax2)) = pyplot.subplots(1,2)

markersize  = 10.0

ax1.set_xlim([x0,x1])
ax1.set_ylim([z0,z1])
ax1.set_aspect("equal")

ax1.scatter(xp,zp,s=markersize)
ax1.add_patch(patches.Rectangle((xbox_center-0.5*xbox_length,zbox_center-0.5*zbox_length),xbox_length,zbox_length,edgecolor="black",facecolor="none",linewidth=2.0))

ax1.set_xlabel(r"$x$")
ax1.set_ylabel(r"$z$")
ax1.set_title("Average relative velocity for particles in this box.")

ax2.set_xlim([y0,y1])
ax2.set_ylim([z0,z1])

ax2.set_aspect("equal")

ax2.scatter(yp,zp,s=markersize)
ax2.add_patch(patches.Rectangle((ybox_center-0.5*ybox_length,zbox_center-0.5*zbox_length),ybox_length,zbox_length,edgecolor="black",facecolor="none",linewidth=2.0))

ax2.set_xlabel(r"$y$")
ax2.set_ylabel(r"$z$")
ax2.set_title("Average relative velocity for particles in this box.")

pyplot.show()
