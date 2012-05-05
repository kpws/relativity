import pylab as pl
import numpy as np
import plotKerrGeo
import kerr
import sympy as sp

M=1
a=0.9*M

e=0.6
h=kerr.Kerr(M,a)

for yi in range(-50,1):
    x=-10.0
    y=float(yi)*.3
    r=np.sqrt(x**2+y**2)
    phi=np.pi-np.arctan(y/-x)
    p=[0,r,0,phi]
    u=sp.symbols('u')
    v=[sp.sympify(0),u*np.cos(phi),sp.sympify(0),-u*np.sin(phi)/r]
    sol=sp.solve(h.invariant(p,v,0)+e,u)
    if len(sol)>1: raise Exception('Ambigous speed')
    else: sol=sol[0]
    y=h.geodesic(p,[s.subs({u:sol}) for s in v],
            np.linspace(0,25,100000), free=0)
    #print(y[0,:])
    #print(y[10,:])
    lb=0;ub=len(y[:,0])
    for i in range(0,len(y),len(y[:,0])/5000):
        if abs(h.invariant(y[i,:4],y[i,4:])+1)>0.0003:
            ub=i
            break
    while ub-lb>1:
        guess=(lb+ub)/2
        if abs(h.invariant(y[guess,:4],y[guess,4:])+1)<0.0003:
            lb=guess
        else:
            ub=guess
    y=y[:ub,:]
    for i in range(len(y)):
        if y[i,1]<0.05*M:
            y=y[:i+1,:]
            break
    if y[-1,1] < M*0.2:
        y[-1,1]=0
    for i in range(-1,2):
        si=h.invariant(y[0,:4],y[0,4:],i)
        ei=h.invariant(y[-1,:4],y[-1,4:],i)
        print("Invariant "+str(i)+" varies"+
               " by "+str(100*abs((si-ei)/si))+" % from start to end.")
    #pl.plot(y[:,7])
    #pl.plot(y[:,3])
    pl.plot( y[:,1] * np.cos(y[:,3]), y[:,1] * np.sin(y[:,3]),'k-')
    #pl.plot( y[:,0],'-',label='Falling from stationary')

plotKerrGeo.plotGeo(M,a)
pl.legend()
pl.show()
