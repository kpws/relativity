import pylab as pl
import numpy as np
import plotKerrGeo
import kerr

M=1
a=0.9*M

h=kerr.Kerr(M,a)

for r in range(-5,5):
    y=h.geodesic( [0,3,0,0], [1,0,0,0.2*r], np.linspace(0,25,10000), free=0)
    #print(y[0,:])
    #print(y[10,:])
    for i in range(-1,2):
        si=h.invariant(y[0,:4],y[0,4:],i)
        ei=h.invariant(y[-1,:4],y[-1,4:],i)
        print("Invariant "+str(i)+" varies"+
               " by "+str(100*abs((si-ei)/si))+" % from start to end.")
    #pl.plot(y[:,7])
    #pl.plot(y[:,3])
    pl.plot( y[:,1] * np.cos(y[:,3]), y[:,1] * np.sin(y[:,3]),'-',label='Falling from stationary')
    #pl.plot( y[:,0],'-',label='Falling from stationary')

plotKerrGeo.plotGeo(M,a)
pl.legend()
pl.show()
