import sympy as sp
import spaceTime
import pickle
from itertools import product

class Kerr(spaceTime.SpaceTime):
    def __init__(self,M0,a0):
        self.x=[self.t,self.r,self.theta,self.phi]=sp.symbols(['t','r','theta','phi'])
        d=self.d=range(len(self.x))
        M,a=sp.symbols(['M','a'])
        rs=2*M
        rho2=self.r**2 + a**2 * sp.cos(self.theta)**2
        delta=self.r**2 - rs*self.r + a**2
        self.g=sp.Matrix([[(delta - a**2 * sp.sin(self.theta)**2) / rho2, 0,
            0, -a * sp.sin(self.theta)**2 * delta / rho2 +
                              (self.r**2 +
                                  a**2)*a*sp.sin(self.theta)**2/rho2],
                          [0, -rho2 / delta, 0, 0],
                          [0, 0, -rho2, 0],
                          [-a * sp.sin(self.theta)**2 * delta / rho2 +
                              (self.r**2 +
                                  a**2)*a*sp.sin(self.theta)**2/rho2,
                              0, 0, a**2 * sp.sin(self.theta)**4 * delta /
                              rho2 - (self.r**2+a**2)**2 *
                              sp.sin(self.theta)**2/rho2]])
        try:
            with open('kerrCache') as f:
                cache=pickle.load(f)
                self.ginv=cache[0]
                self.CS=cache[1]
        except IOError:
            gdet=self.g.det().factor().trigsimp()
            cfm=self.g.cofactorMatrix()
            self.ginv=sp.Matrix([[
                (cfm[i,j].factor().trigsimp()/gdet).factor().trigsimp() for j in d] for i in d])
            self.makeCS()
            for i,j,k in product(d,d,d):
                self.CS[i][j][k]=self.CS[i][j][k].factor().trigsimp()
            with open('kerrCache','w') as f:
                pickle.dump([self.ginv,self.CS],f)
        subDic={M:M0,a:a0}
        self.g=self.g.subs(subDic)
        self.ginv=self.ginv.subs(subDic)
        for i,j,k in product(d,d,d):
            self.CS[i][j][k]=self.CS[i][j][k].subs(subDic)

