import sympy as sp
import spaceTime
import pickle
from itertools import product

class Kerr(spaceTime.SpaceTime):
    def __init__(self,M,a):
        #these are Kerrs original coords (advanced Eddington-Finkelstein for
        #rotation BG) but with the substitution chi=cos(theta)
        self.x=[u,r,chi,phi]=[self.u,self.r,self.chi,self.phi]=sp.symbols(['u','r','chi','phi'])
        d=self.d=range(len(self.x))
        M0,a0=sp.symbols(['M0','a0'])
        [du,dr,dchi,dphi]=sp.symbols(['du','dr','dchi','dphi'])
        self.makeg( -( 1 - 2*M0*r / (r**2 + a0**2*chi**2) ) * (du + a0*(1-chi**2)*dphi)**2 + 
                     2*(du + a0*(1-chi**2)*dphi ) * (dr + a0*(1-chi**2)*dphi) +
                     (r**2 + a0**2*chi**2) * (dchi**2/(1-chi**2) + (1-chi**2)*dphi**2) ,
                     [du,dr,dchi,dphi])
        try:
            with open('kerrCache') as f:
                cache=pickle.load(f)
                if self.g!=cache[0]:
                    raise IOError()
                self.ginv=cache[1]
                self.CS=cache[2]
        except IOError:
            gdet=self.g.det().factor()
            cfm=self.g.cofactorMatrix()
            self.ginv=sp.Matrix([[
                (cfm[i,j].factor().trigsimp()/gdet).factor() for j in d] for i in d])
            self.makeCS()
            for i,j,k in product(d,d,d):
                self.CS[i][j][k]=self.CS[i][j][k].factor()
            with open('kerrCache','w') as f:
                pickle.dump([self.g,self.ginv,self.CS],f)
        subDic={M0:M,a0:a}
        self.g=self.g.subs(subDic)
        self.ginv=self.ginv.subs(subDic)
        for i,j,k in product(d,d,d):
            self.CS[i][j][k]=self.CS[i][j][k].subs(subDic)
        self.killing=[[1,0,0,0],[0,0,0,1]]

