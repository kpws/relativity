import sympy as sp
from kerr import *

[M,a]=sp.symbols(['M','a'])
h=Kerr(M,a)

sp.pprint(h.ginv[0,0])
