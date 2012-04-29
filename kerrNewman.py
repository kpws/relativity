import sympy as sp
import spaceTime

class KerrNewman(spaceTime.SpaceTime):
    def __init__(self,M,alpha,Q):
        self.x=[self.t,self.r,self.theta,self.phi]=sp.symbols(['t','r','theta','phi'])
        rs=2*M
        rq2=Q**2/(4*sp.pi)
        rho2=self.r**2 + alpha**2 * sp.cos(self.theta)**2
        delta=self.r**2 - rs*self.r + alpha**2 + rq2
        self.g=sp.Matrix([[(delta - alpha**2 * sp.sin(self.theta)**2) / rho2, 0,
            0, -alpha * sp.sin(self.theta)**2 * delta / rho2 +
                              (self.r**2 +
                                  alpha**2)*alpha*sp.sin(self.theta)**2/rho2],
                          [0, -rho2 / delta, 0, 0],
                          [0, 0, -rho2, 0],
                          [-alpha * sp.sin(self.theta)**2 * delta / rho2 +
                              (self.r**2 +
                                  alpha**2)*alpha*sp.sin(self.theta)**2/rho2,
                              0, 0, alpha**2 * sp.sin(self.theta)**4 * delta /
                              rho2 - (self.r**2+alpha**2)**2 *
                              sp.sin(self.theta)**2/rho2]])
        super(KerrNewman,self).__init__()
