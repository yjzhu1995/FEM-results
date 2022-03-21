import numpy as np
# import pygmt
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from scipy.interpolate import RectBivariateSpline
from mpl_toolkits.axes_grid1 import make_axes_locatable



# main program 
m2cm = 100.0

file = "./displacements.dat"
nongrav = np.loadtxt(file, skiprows=0)

file = "/mnt/d/myprograms/benchmark/dipslip_w100/okada/displacement_Okada_wang.dat"
okada = np.loadtxt(file, skiprows=1)

fig = plt.figure()
plt.subplot(211)
plt.plot(nongrav[:,0]*111.2, nongrav[:,2]*m2cm, '-', color='green', label='layered')
plt.plot(okada[:,0]*111.2, okada[:,2]*m2cm, '--', color='black', label='okada')
plt.xlim(left=-200, right=200)
plt.xlabel('km')
plt.ylabel('E-W displacements (cm)')
plt.legend(loc='best')
plt.subplot(212)
plt.plot(nongrav[:,0]*111.2, nongrav[:,4]*m2cm, '-', color='green', label='layered')
plt.plot(okada[:,0]*111.2, okada[:,4]*m2cm, '--', color='black', label='okada')
plt.xlim(left=-200, right=200)
plt.xlabel('km')
plt.ylabel('vertcial displacements (cm)')
plt.legend(loc='best')
fig.show()
fig.savefig(fname='nongrav_compare.png',format='png',dpi=720)
