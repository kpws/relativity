import pylab as pl
import numpy as np
import plotKerrGeo
import kerr

M=1
a=0.9

h=kerr.Kerr(M,a)
plotKerrGeo.plotGeo(M,a)

for r in range(6,10):
    y=h.geodesic( [0,r,np.pi/2,0], [1,-1,0,10], np.linspace(0,100,1000))
    pl.plot( y[:,1] * np.cos(y[:,3]), y[:,1] * np.sin(y[:,3]) ,label='Falling from stationary')

pl.legend()
pl.show()
