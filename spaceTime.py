import sympy as sp
import numpy as np
from scipy.integrate import odeint
from sympy.utilities.lambdify import lambdify
from sympy.solvers import solve

class SpaceTime(object):
    def __init__(self):
        self.ginv=self.g.inv()
        d=range(len(self.x))
        
        self.makeCS()

        #Ricci curvature_a_b  through formula 2.33 in hartle
        self.R2=[[sum(self.CS[c][a][b].diff(self.x[c])
           -self.CS[c][a][c].diff(self.x[b])
           +sum(self.CS[c][a][b]*self.CS[e][c][e]
               -self.CS[e][a][c]*self.CS[c][b][e] for e in d) for c in d) for b in d] for a in d]

    def makeg(self,lineElement,diffs):
        d=range(len(self.x))
        self.g=sp.Matrix([[lineElement.diff(diffs[i],diffs[j])/2 for j in d]for i in d])

    #Christoffel symbol^a_b_g
    def makeCS(self):
        d=range(len(self.x))
        self.CS=[[[sum(self.ginv[e,a]*(self.g[e,b].diff(self.x[c])
                        +self.g[e,c].diff(self.x[b])
                        -self.g[b,c].diff(self.x[e])) for e in d)/2
          for c in d] for b in d] for a in d]

    #only time-like geodesics implemented so far
    def geodesic(self,x,u,tau,free=-1):
        d=range(len(self.x))
        if free!=-1:
            u[free]=sp.symbols('x')
            sol=solve(sum(sum(u[i]*u[j]*self.g[i,j].subs(zip(self.x,x)) for j in
                d) for i in d)+1, u[free])
            sol=[s for s in sol if s.is_real]
            if sol==[]:
                raise Exception("Couldn't find time-like solution")
            u[free]=float(max(sol))
        n=len(x)
        y0=x+u
        CSf=[[[lambdify(self.x,self.CS[i][j][k]) for k in d]for j in d]for i in d]
        def yprim(y,dummy):
            return [  y[n+i] if i < n else
                      -sum( sum( CSf[i-n][j][k](*y[:n]) * y[n+k] for k in d) * y[n+j] for j in d)
                           for i in range(2*n) ]
        return odeint(yprim, y0, tau)

    def invariant(self,x,u,k=-1):
        v = u if k==-1 else self.killing[k]
        return sum(sum(v[i]*u[j]*self.g[i,j].subs(zip(self.x,x)) for j in self.d)for i in self.d)
