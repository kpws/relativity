import sympy as sp
from rel import *

c=[t,x,y,z]=sp.symbols(['t','x','y','z'])
d=range(len(c))

f,eps=sp.symbols(['f','eps'])
h=sp.Matrix([[0,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,0]])*f(x-t)*eps

sp.pprint()
g=eta+h
sp.pprint(g)

R=sp.simplify(R(c,CS(c,g)))
for a in d:
    for b in range(a+1):
        print('R('+str(c[a])+','+str(c[b])+')='+str(R[a][b].series(eps,n=2)))

