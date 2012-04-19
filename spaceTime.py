import sympy as sp

class SpaceTime(object):
    def __init__(self):
        self.ginv=self.g.inv()
        d=self.d=range(len(self.x))
        #Christoffel symbol^a_b_g
        self.CS=[[[sum(self.ginv[e,a]*(self.g[e,b].diff(self.x[c])
                        +self.g[e,c].diff(self.x[b])
                        -self.g[b,c].diff(self.x[e])) for e in d)/2
          for c in d] for b in d] for a in d]

        #Ricci curvature_a_b  through formula 2.33 in hartle
        self.R2=[[sum(self.CS[c][a][b].diff(self.x[c])
           -self.CS[c][a][c].diff(self.x[b])
           +sum(self.CS[c][a][b]*self.CS[e][c][e]
               -self.CS[e][a][c]*self.CS[c][b][e] for e in d) for c in d) for b in d] for a in d]

