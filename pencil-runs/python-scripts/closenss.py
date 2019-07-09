import pencil
import numpy
from matplotlib import pyplot
import sys
import os

base_path,script_name   = os.path.split(sys.argv[0])
scratch,simulation_name = os.path.split(base_path)
script_name             = script_name[:-3]

ivar = -1
pvar = "pvar.dat"
if len(sys.argv) > 1:
    ivar = int(sys.argv[1])
    pvar = "PVAR" + sys.argv[1]

parameters = pencil.read_param()
data       = pencil.read_var(ivar=ivar)
pdata      = pencil.read_pvar(varfile=pvar)

xp    = pdata.xp*1000.
yp    = pdata.yp*1000.
zp    = pdata.zp*1000.
ipars = pdata.ipars

closeness = numpy.zeros(len(ipars))
for i in range(len(ipars)):
    x_ref = xp[i]
    y_ref = yp[i]
    z_ref = zp[i]
    for j in range(len(ipars)):
        if(j == i):
            continue
        else:
            x = xp[j]
            y = yp[j]
            z = zp[j]
            closeness[i] = closeness[i] + ((x - x_ref)**2 + (y - y_ref)**2 + (z - z_ref)**2)**(-0.5)

print("Min closeness: " + str(closeness.min()))
print("Mean closeness: " + str(numpy.mean(closeness)))
print("Median closeness: " + str(numpy.median(closeness)))
print("Max closeness: " + str(closeness.max()))
