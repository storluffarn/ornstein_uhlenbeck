
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as colors
import time 

traw = np.loadtxt(open("time0.csv", "rb"), delimiter=",", skiprows=0)
xraw = np.loadtxt(open("xout0.csv", "rb"), delimiter=",", skiprows=0)
#lraw = np.loadtxt(open("lout.csv", "rb"), delimiter=",", skiprows=1)
#fraw = np.loadtxt(open("tomout.csv", "rb"), delimiter=",", skiprows=1)
#params = np.loadtxt(open("parameters.dat", "rb"), delimiter=",", skiprows=1)

nscale = 1.0e9  # WARNING this is *only* for plotting, danger!

#tdata = nscale * traw[:,0]
#supdata = nscale * traw[:,1]
#sdata = nscale * traw[:,1]
#xdata = nscale * xraw[:,0]
##ldata = nscale * lraw[:,0]
#fdata = nscale * fraw[:,2]

tdata = nscale * traw
xdata = nscale * xraw
#fdata = nscale * fraw

#fig0, ax = plt.subplots()
#ax.plot(tdata,supdata)
#ax.set(xlabel='t (ns)', ylabel='f (nm)')
#fig0.savefig("plots/plot_sup.png")

fig1, ax = plt.subplots()
#ax.plot(tdata,xdata[200:400],'.')
ax.plot(tdata,xdata,'.', markersize = 2)
ax.set(xlabel='t (ns)', ylabel='x (nm)')
fig1.savefig("plots/plot_x.png")

#
#fig2, ax = plt.subplots()
#ax.plot(tdata,ldata)
#ax.set(xlabel='t (ns)', ylabel='l (nm)')
#fig2.savefig("plots/plot_l.png")

#fig3, (ax1,ax2) = plt.subplots(nrows = 2, sharex = True)
#
#ax1.plot(tdata,ldata)
#ax1.set(ylabel='l (nm)')
#ax1.xaxis.set_major_locator(ticker.MultipleLocator(2.0))
#ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
#ax1.xaxis.grid(True, which='minor', linestyle='dotted')
#ax1.xaxis.grid(True, which='major', linestyle='dotted')
#ax1.yaxis.grid(True, which='major', linestyle='dotted')
#fig3.suptitle("k1 = {:.2e}, k2 = {:.2e}, t0 = {:.2e}, \n gamma_x = {:.2e}, gamma_l = {:.2e}".format(params[0],params[1], params[2] ,params[17],params[18]), fontsize=14)
#
#ax2.plot(tdata,xdata)
#ax2.set(xlabel='t (nm)', ylabel='x (nm)')
#ax2.xaxis.set_major_locator(ticker.MultipleLocator(20.0))
#ax2.xaxis.set_minor_locator(ticker.MultipleLocator(5.0))
#ax2.xaxis.grid(True, which='minor', linestyle='dotted')
#ax2.xaxis.grid(True, which='major', linestyle='dotted')
#ax2.yaxis.grid(True, which='major', linestyle='dotted')
#
#plt.subplots_adjust(hspace=0)
#
#timestamp = time.strftime("%Y%m%d-%H%M%S")
#
#fig3.savefig("plots/plot_lx{}".format(timestamp))

#fig4, ax = plt.subplots()
#ax.plot(tdata,fdata)
#ax.set(xlabel='t (ns)', ylabel='f (nN)')
#fig4.savefig("plots/plot_f.png")
