import sympy as sp
import spaceTime

class Flat(spaceTime.SpaceTime):
    def __init__(self):
        self.x=sp.symbols(['t','x','y','z'])        
        self.g=sp.Matrix([[-1,0,0,0],
                     [0,1,0,0],
                     [0,0,1,0],
                     [0,0,0,1]])
        super(Flat,self).__init__()
