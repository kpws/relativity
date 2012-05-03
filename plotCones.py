import pylab as pl
import numpy as np
import plotKerrGeo
import kerr

M=1
a=0.7

h=kerr.Kerr(M,a)

for r in range(3,6):
    y=h.geodesic( [0,r,0,0], [1,0,0,0], np.linspace(0,200,5000), free=0)
    print(y[0,:])
    print(y[10,:])
    print(h.invariant(y[0,:4],y[0,4:]))
    print(h.invariant(y[10,:4],y[10,4:]))
    pl.plot(y[:,7])
    pl.plot(y[:,3])
    #pl.plot( y[:,1] * np.cos(y[:,3]), y[:,1] * np.sin(y[:,3]),'-',label='Falling from stationary')
    #pl.plot( y[:,0],'-',label='Falling from stationary')

#plotKerrGeo.plotGeo(M,a)
pl.legend()
pl.show()
