import sympy as sp
import spaceTime
import kerrNewman

class Kerr(kerrNewman.KerrNewman):
    def __init__(self,M,J):
        super(Kerr,self).__init__(M,J,0)
