import pencil
import numpy
from matplotlib import pylab
import sys

log_y = False

ts = pencil.read_ts()
ts_dictionary = ts.__dict__

start_params = pencil.read_param()
grav = start_params.gravx

#if(len(sys.argv) == 1):
#    print("Need timeseries value to plot. Options are: ")
#    print(ts_dictionary.keys())
#    quantity_name = raw_input("Please make a selection, or type Q to quit: ")
#    if(quantity_name == "Q"):
#        sys.exit("Quitting...")
#else:
#    quantity_name = sys.argv[1]

quantity_name = "vpxm"

try:
    quantity = ts_dictionary[quantity_name]
except KeyError:
    sys.exit(quantity_name + " is not the name of a quantity in the time series.")

figure  = pylab.figure()
subplot = figure.add_subplot(111)
time    = ts.t/0.014
subplot.plot(time,quantity,linestyle="none",marker=".",markeredgecolor="black",markerfacecolor="none",markersize=0.1)
subplot.plot(time,ts.t*2.45,linestyle="--",color="gray")
subplot.plot(time,numpy.repeat(2.45*0.014,len(time)),linestyle="-",color="green")

subplot.set_ylim([0.0,0.06])

if(len(sys.argv) == 3):
    log_y_str = sys.argv[2]
    if(log_y_str == "True" or log_y_str == "true" or log_y_str == "T"):
        log_y = True
    elif(not (log_y_str == "False" or log_y_str == "false" or log_y_str == "F")):
        print("Second argument does not match format (True, true, or T) or (False, false, or F).")
        print("Assuming log_y = False")
if(log_y):
    subplot.set_yscale("log")

subplot.set_title(quantity_name + " vs t")
subplot.set_ylabel(quantity_name)
subplot.set_xlabel("t (friction times)")

pylab.show()
