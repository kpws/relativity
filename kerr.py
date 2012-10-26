import sympy as sp
import spaceTime
import pickle
from itertools import product

class Kerr(spaceTime.SpaceTime):
    def __init__(self,M,a,coords="KerrChi"):
        M0,a0=sp.symbols(['M0','a0'])
        if coords=="KerrChi":
            #these are Kerrs original coords (advanced Eddington-Finkelstein but for
            #rotating BH) but with the substitution chi=cos(theta)
            self.x=[u,r,chi,phi]=[self.u,self.r,self.chi,self.phi]=sp.symbols(['u','r','chi','phi'])
            [du,dr,dchi,dphi]=sp.symbols(['du','dr','dchi','dphi'])
            self.makeg( -( 1 - 2*M0*r / (r**2 + a0**2*chi**2) ) * (du + a0*(1-chi**2)*dphi)**2 + 
                         2*(du + a0*(1-chi**2)*dphi ) * (dr + a0*(1-chi**2)*dphi) +
                         (r**2 + a0**2*chi**2) * (dchi**2/(1-chi**2) + (1-chi**2)*dphi**2) ,
                         [du,dr,dchi,dphi])
            self.killing=[[1,0,0,0],[0,0,0,1]]
        elif coords=='BLChi':
            self.x=[t,r,chi,phi]=[self.t,self.r,self.chi,self.phi]=sp.symbols(['t','r','chi','phi'])
            [dt,dr,dchi,dphi]=sp.symbols(['dt','dr','dchi','dphi'])
            rho2=r**2+a0**2*chi**2
            delta=r**2+a0**2-2*M0*r
            self.makeg( -( 1 - 2*M0*r / rho2)*dt**2
                        -(4*M0*a0*r*(1-chi**2))/rho2*dt*dphi
                        +rho2/delta*dr**2+rho2*dchi**2/(1-chi**2)
                        +(r**2+a0**2+(2*M0*r*a0**2*(1-chi**2))/rho2)*(1-chi**2)*dphi**2,
                         [dt,dr,dchi,dphi])
            self.killing=[[1,0,0,0],[0,0,0,1]]
        else:
            raise Exception("Kerr BH in "+coords+" is not implemented yet.")
        d=self.d=range(len(self.x))
        try:
            with open('kerrCache_'+coords) as f:
                cache=pickle.load(f)
                if self.g!=cache[0]:
                    raise IOError()
                self.ginv=cache[1]
                self.CS=cache[2]
                self.dCSdx=cache[3]
        except IOError:
            gdet=self.g.det().factor()
            cfm=self.g.cofactorMatrix()
            self.ginv=sp.Matrix([[
                (cfm[i,j].factor().trigsimp()/gdet).factor().trigsimp() for j in d] for i in d])
            self.makeCS()
            for i,j,k in product(d,d,d):
                self.CS[i][j][k]=self.CS[i][j][k].factor().trigsimp()
            self.makedCSdx()
            for i,j,k,l in product(d,d,d,d):
                self.dCSdx[i][j][k][l]=self.dCSdx[i][j][k][l].factor().trigsimp()
            with open('kerrCache_'+coords,'w') as f:
                pickle.dump([self.g,self.ginv,self.CS,self.dCSdx],f)
        subDic={M0:M,a0:a}
        self.g=self.g.subs(subDic)
        self.ginv=self.ginv.subs(subDic)
        for i,j,k in product(d,d,d):
            self.CS[i][j][k]=self.CS[i][j][k].subs(subDic)
        for i,j,k,l in product(d,d,d,d):
            self.dCSdx[i][j][k][l]=self.dCSdx[i][j][k][l].subs(subDic)
        self.makeNumericalFunctions()
