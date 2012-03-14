import sympy as sp

#Coordinates
x=[t,r,theta,phi]=sp.symbols(['t','r','theta','phi'])
d=range(len(x))

#Schwarzschild metric_a_b
M=sp.symbols('M')
g=sp.Matrix([[-(1-2*M/r),0,0,0],
             [0,(1-2*M/r)**-1,0,0],
             [0,0,r**2,0],
             [0,0,0,r**2*sp.sin(theta)**2]])

#Christoffel symbol^alpha_beta_gamma
CS=[[[sum(g.inv()[e,a]*(g[e,b].diff(x[c])
                       +g[e,c].diff(x[b])
                       -g[b,c].diff(x[e])) for e in d)/2
      for c in d] for b in d] for a in d]

#Ricci curvature_a_b  through formula 2.33 in book
R=[[sum(CS[c][a][b].diff(x[c])
       -CS[c][a][c].diff(x[b])
       +sum(CS[c][a][b]*CS[e][c][e]
           -CS[e][a][c]*CS[c][b][e] for e in d) for c in d) for b in d] for a in d]

#Print the results
for a in d:
    for b in range(a+1):
        print('R('+str(x[a])+','+str(x[b])+')='+str(R[a][b].simplify()))
