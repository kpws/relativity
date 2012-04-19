import sympy as sp
import spaceTime

class KerrNewman(spaceTime.SpaceTime):
    def __init__(M,J,Q):
        self.x=[self.t,self.r,selft.theta,self.phi]=sp.symbols(['t','r','theta','phi'])
        rs=2*M
        rq2=Q**2/(4*sp.pi)
        alpha=J/M
        rho2=r**2 + alpha**2 * sp.cos(self.theta)**2
        self.g=sp.Matrix([[delta/rho2,0,0,0],
                          [0,1,0,0],
                          [0,0,1,0],
                          [0,0,0,1]])
        super(KerrNewman,self).__init__()
