import sympy as sp
import numpy as np
from scipy.integrate import odeint
from sympy.utilities.lambdify import lambdify
from itertools import product

class SpaceTime(object):
    def __init__(self):
        self.ginv=self.g.inv()
        d=self.d=range(len(self.x))
        
        self.makeCS()

        #Ricci curvature_a_b  through formula 2.33 in hartle
        self.R2=[[sum(self.CS[c][a][b].diff(self.x[c])
           -self.CS[c][a][c].diff(self.x[b])
           +sum(self.CS[c][a][b]*self.CS[e][c][e]
               -self.CS[e][a][c]*self.CS[c][b][e] for e in d) for c in d) for b in d] for a in d]

    #Christoffel symbol^a_b_g
    def makeCS(self):
        d=self.d
        self.CS=[[[(sum(self.ginv[e,a]*(self.g[e,b].diff(self.x[c])
                        +self.g[e,c].diff(self.x[b])
                        -self.g[b,c].diff(self.x[e])) for e in d)/2)
          for c in d] for b in d] for a in d]

    #only time-like geodesics implemente so far
    def geodesic(self,x,u,tau,null=False):
        d=self.d
        uu=float(sum(sum(u[i]*u[j]*self.g[i,j].subs(zip(self.x,x)) for j in d) for i
                in d))
        if not null:
            if uu>0:
                raise Exception('Velocity is space-like')
            for i in d:
                u[i]/=np.sqrt(-uu)
        n=len(x)
        y0=x+u
        CSf=[[[lambdify(self.x,self.CS[i][j][k]) for k in d]for j in d]for i in d]
        def yprim(y,dummy):
            return y[n:]+[ -sum( sum( CSf[i][j][k](*y[:n]) * y[n+k] for k in d) * y[n+j] for j in d) for i in d ]
        return odeint(yprim, y0, tau)








