import sympy as sp


eta=sp.Matrix([[-1,0,0,0],
               [0,1,0,0],
               [0,0,1,0],
               [0,0,0,1]])

#Christoffel symbol^a_b_g
def CS(x,g):
    d=range(len(x))
    ginv=g.inv()
    return [[[sum(ginv[e,a]*(g[e,b].diff(x[c])
                    +g[e,c].diff(x[b])
                    -g[b,c].diff(x[e])) for e in d)/2
      for c in d] for b in d] for a in d]

#Ricci curvature_a_b  through formula 2.33 in book
def R(x,CS):
    d=range(len(x))
    return [[sum(CS[c][a][b].diff(x[c])
       -CS[c][a][c].diff(x[b])
       +sum(CS[c][a][b]*CS[e][c][e]
           -CS[e][a][c]*CS[c][b][e] for e in d) for c in d) for b in d] for a in d]

