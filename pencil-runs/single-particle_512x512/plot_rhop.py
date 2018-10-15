import pencil
import numpy
import matplotlib.pyplot as pyplot
import sys

var = -1

if len(sys.argv) == 2:
    var = int(sys.argv[1])

data = pencil.read_var(ivar=var,trimall=True)
#datap = pencil.read_pvar(ivar=var)
x2d,y2d = numpy.meshgrid(data.x,data.y)

fig, ((ax)) = pyplot.subplots(1,1)
ax.contourf(x2d,y2d,data.rhop,256)

ax.set_xlabel("x")
ax.set_ylabel("y")

pyplot.show()
