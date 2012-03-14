#    Copyright 2012 Petter Saeterskog
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
