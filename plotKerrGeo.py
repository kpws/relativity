import pylab as pl
import numpy as np
import sympy as sp

def plotGeo(M,a,text=''):
    rp=M+np.sqrt(M**2-a**2)
    rm=M-np.sqrt(M**2-a**2)
    re=M+np.sqrt(M**2-a**2*0)
    OmegaH=a/(2*M*rp)

    phi=np.linspace(-np.pi,np.pi,300)
    pl.plot(rp*np.cos(phi),rp*np.sin(phi),'b',label='Outer event horizon')
    pl.plot(rm*np.cos(phi),rm*np.sin(phi),'--r',label='Inner event horizon')
    pl.plot(re*np.cos(phi),re*np.sin(phi),'--g',label='Ergosphere boundary')

    pl.axis('equal')
    tickR=100
    tickMarks=['' if (i/2)*2!=i else '$'+sp.latex(sp.sympify(str(i)+'*M'))+'$' for i in range(-tickR,tickR+1)]
    pl.xticks( np.arange(-tickR,tickR+1), tickMarks )
    pl.yticks( np.arange(-tickR,tickR+1), tickMarks )

    axisR=re*1.2
    pl.axis([-axisR,axisR,-axisR,axisR])
    pl.xlabel('$r\cos(\phi)$')
    pl.ylabel('$r\sin(\phi)$')
    #pl.title('$\\theta=\pi/2$-plane of Kerr black-hole geometry, Boyer-Lindquist'+
    #'coordinates, $M=1$, $a='+str(a)+'$, '+text)
