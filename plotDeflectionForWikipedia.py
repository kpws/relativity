import pylab as pl
import numpy as np
import plotKerrGeo
import kerr
import sympy as sp
import math

M=1.
e=1.
xStart=-10*M
yStart=-8*M
endTau=25

for cw in [False, True]:

    a=(1 if cw else -1)*0.9*M
    h=kerr.Kerr(M,a)
    pl.figure()
    for yi in range(-32,1,2):
        x=xStart
        y=float(yi)*0.5
        r=np.sqrt(x**2+y**2)
        phi=np.pi-np.arctan(y/-x)
        p=[0,r,0,phi]
        u=sp.symbols('u')
        v=[sp.sympify(0),u*np.cos(phi),sp.sympify(0),-u*np.sin(phi)/r]
        sol=sp.solve(h.invariant(p,v,0)+e,u)
        if len(sol)>1: raise Exception('Ambigous speed')
        else: sol=sol[0]
        tau,y=h.geodesic(p,[s.subs({u:sol}) for s in v], np.linspace(0,endTau,1e5), free=0)
        for i in range(-1,2):
            si=h.invariant(y[0,:4],y[0,4:],i)
            ei=h.invariant(y[-2,:4],y[-2,4:],i)
            print("Relative variation of invariant "+str(i)+
                    ": "+str(abs((si-ei)/si))+", from start to end.")
        
        plotCoords=zip( y[:,1] * np.cos(y[:,3]), y[:,1] * np.sin(y[:,3]))
        #remove all points outside plot because matplotlib doesn't crop the svn image correctly
        plotCoords=[p for p in plotCoords if xStart<=p[0]<=-xStart and yStart<=p[1]<=-yStart]
        pl.plot([p[0] for p in plotCoords],[ p[1] for p in plotCoords],'k-')
        pTau=[t for t in np.linspace(0,endTau, endTau+1) if t<=tau[-1]]
        plotCoords=zip(np.interp(pTau,tau,y[:,1]*np.cos(y[:,3])), np.interp(pTau,tau,y[:,1] * np.sin(y[:,3])))
        #remove all points outside plot because matplotlib doesn't crop the svn image correctly
        plotCoords=[p for p in plotCoords if xStart<=p[0]<=-xStart and yStart<=p[1]<=-yStart]
        pl.plot([p[0] for p in plotCoords],[ p[1] for p in plotCoords],'ko')
       
    plotKerrGeo.plotGeo(M,a)
    pl.legend(loc=2)
    pl.xlim(xStart,-xStart)
    pl.ylim(yStart,-yStart)
    pl.savefig(('cw' if cw else 'ccw')+'.svg', bbox_inches=0)
pl.show()
