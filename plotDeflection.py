import pylab as pl
import numpy as np
import plotKerrGeo
import kerr
import sympy as sp
import math

M=1
a=-0.9*M

e=1
h=kerr.Kerr(M,a)

for yi in range(-20,1):
    x=-15.0
    y=float(yi)*0.5
    r=np.sqrt(x**2+y**2)
    phi=np.pi-np.arctan(y/-x)
    p=[0,r,0,phi]
    u=sp.symbols('u')
    v=[sp.sympify(0),u*np.cos(phi),sp.sympify(0),-u*np.sin(phi)/r]
    sol=sp.solve(h.invariant(p,v,0)+e,u)
    if len(sol)>1: raise Exception('Ambigous speed')
    else: sol=sol[0]
    y=h.geodesic(p,[s.subs({u:sol}) for s in v],
            np.linspace(0,25,50000), free=0)
    for i in range(0,len(y),len(y[:,0])/5000):
        if abs(h.invariant(y[i,:4],y[i,4:])+1)>0.000003:
            y=y[:i-1,:]
            break
    lb=0;ub=len(y[:,0])
    while ub-lb>1:
        guess=(lb+ub)/2
        err=abs(h.invariant(y[guess,:4],y[guess,4:])+1)
        if err<0.000003 and not math.isnan(err):
            lb=guess
        else:
            ub=guess
    y=y[:lb-3,:]
    xx=np.cos(y[:,3])*y[:,1]
    yy=np.sin(y[:,3])*y[:,1]
    for i in range(1,len(y)-1):
        dx1=xx[i]-xx[i-1]
        dx2=xx[i+1]-xx[i]
        dy1=yy[i]-yy[i-1]
        dy2=yy[i+1]-yy[i]
        angle=np.arccos((dx1*dx2+dy1*dy2)/np.sqrt((dx1*dx1+dy1*dy1)*(dx2*dx2+dy2*dy2)))
        if y[i,1]<0.05*M or angle>0.2:
            y=y[:i+1,:]
            break
    if y[-1,1] < M*0.2:
        y[-1,1]=0
    for i in range(-1,2):
        si=h.invariant(y[0,:4],y[0,4:],i)
        ei=h.invariant(y[-2,:4],y[-2,4:],i)
        print("Invariant "+str(i)+" varies"+
               " by "+str(100*abs((si-ei)/si))+" % from start to end.")
    #pl.plot(y[:,7])
    #pl.plot(y[:,3])
    pl.plot( y[:,1] * np.cos(y[:,3]), y[:,1] * np.sin(y[:,3]),'k-')
    #pl.plot( y[:,0],'-',label='Falling from stationary')
#pl.plot([ h.invariant(y[i,:4],y[i,4:]) for i in range(len(y))])
   
plotKerrGeo.plotGeo(M,a, 'Trajectories for counterrotating particles with $e=1$')
pl.legend(loc=2)
pl.show()
