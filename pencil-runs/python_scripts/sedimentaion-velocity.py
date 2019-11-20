import pencil
import numpy
from matplotlib import pylab
import sys

params = pencil.read_param()
ts     = pencil.read_ts()

tau = 0.014
g   = 2.45

v_s      = -tau*g
v_z_max  = ts.vpzmax
v_z_min  = ts.vpzmin
#v_z_mean = ts.vpzm
time     = ts.t/tau

figure  = pylab.figure()
subplot = figure.add_subplot(111)
subplot.plot(time,v_z_max,color="red",label="Maximum z-velocity")
subplot.plot(time,v_z_min,color="green",label="Minimum z-velocity")
#subplot.plot(time,v_z_min,color="black",label="Mean z-velocity")
subplot.plot(time,numpy.repeat(v_s,len(time)),linestyle="--",color="grey",label="Isolated terminal velocity")

subplot.set_title("vpzmax vs t")
subplot.set_ylabel("vpz")
subplot.set_xlabel("t (friction times)")

subplot.legend()
pylab.show()
