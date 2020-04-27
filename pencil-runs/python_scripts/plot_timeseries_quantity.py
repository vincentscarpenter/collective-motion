import pencil
import numpy
from matplotlib import pylab
import sys

log_y = False

ts = pencil.read_ts()
ts_dictionary = ts.__dict__

if(len(sys.argv) == 1):
    print("Need timeseries value to plot. Options are: ")
    print(ts_dictionary.keys())
    quantity_name = raw_input("Please make a selection, or type Q to quit: ")
    if(quantity_name == "Q"):
        sys.exit("Quitting...")
else:
    quantity_name = sys.argv[1]
try:
    quantity = ts_dictionary[quantity_name]
except KeyError:
    sys.exit(quantity_name + " is not the name of a quantity in the time series.")

t0 = ts.t[0]
t1 = ts.t[-1]
if(len(sys.argv) == 4):
    t0 = float(sys.argv[2])
    t1 = float(sys.argv[3])
    while(t0 < ts.t[0] or t1 > ts.t[-1]):
        print("Times are out of domain; simulation times go from " + str(ts.t[0]) + " to " + str(ts.t[-1]) + ". Please pick new limits, or enter Q to quit.")
        t0 = raw_input("Lower limit for time: ")
        if(t0 == "Q"):
            sys.exit("Quitting...")
        else:
            t0 = float(t0)
        t1 = raw_input("Upper limit for time: ")
        if(t1 == "Q"):
            sys.exit("Quitting...")
        else:
            t1 = float(t1)
times = numpy.where((ts.t > t0) & (ts.t < t1))

figure  = pylab.figure()
subplot = figure.add_subplot(111)
subplot.plot(ts.t[times],quantity[times])

subplot.set_title(quantity_name + " vs t")
subplot.set_ylabel(quantity_name)
subplot.set_xlabel("t (code units)")

pylab.show()
