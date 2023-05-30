import numpy as np
from numpy.linalg import *
from matplotlib import pyplot as plt
from scipy.spatial import ConvexHull

def gcde(a, b):
    # Extended euclidean algorithm
    if a == 0:
        return b, 0, 1
    else:
        gcd, x, y = gcde(b % a, a)
        return gcd, y - (b // a) * x, x

class Edge:
    def __init__(self,λ,c):
        self.λ = np.array(λ)
        self.c = c
        self.ends = []
        self.neighbors = []
    def __repr__(self):
        return f"({self.λ[0]} {self.λ[1]})·(x,y) {self.c:+} ≥ 0"

    def l(self,x):
        return np.matmul(self.λ , x) - self.c
    
    def find_default_probe(self):
        _,s,t = gcde(-self.λ[1], self.λ[0])
        p = np.array([-t,s])
        return np.sign(np.matmul(self.λ,p))*p
        

def intersect_edges(e1,e2):
    M = np.block([[e1.λ],[e2.λ]])
    if det(M)!=0:
        return np.matmul(inv(M),np.array([e1.c,e2.c]))
    return None

class Polygon:
    def __init__(self, points = None, edges = None):
        if points is not None:
            self.vertices = points
            self.update_edges()

        if edges is not None:
            self.edges = edges
            self.vertices = []
            self.update_vertices()

    def __repr__(self):
        s=""
        for e in self.edges:
            s = s + str(e) + "\n"
        return s

    def update_edges(self):
        print(self.vertices)
        hull = ConvexHull(self.vertices)
        self.vertices = hull.points[hull.vertices]
        self.edges = []
        for e in hull.equations:
            self.edges.append(Edge([round(v) for v in e[:2]],e[2]))
        for i in range(len(hull.simplices)):
            n0, n1 = tuple(hull.neighbors[i])
            self.edges[i].neighbors = [ self.edges[j] for j in hull.neighbors[i] ]
            self.edges[i].ends = [intersect_edges(self.edges[i], e2) for e2 in self.edges[i].neighbors ]

    def update_vertices(self):
        self.vertices = []
        for e0 in self.edges:
            for e1 in self.edges:
                a = intersect_edges(e0, e1)
                if a is not None:
                    self.vertices.append(a.copy())
        for e in self.edges:
            self.vertices = [ c for c in self.vertices if e.l(c) >= 0 ]
        print(self.vertices)
        self.update_edges()
    
    def draw(self):
        x,y = ([ c[i] for c in self.vertices + [self.vertices[0]] ] for i in range(2))
        fig,ax = plt.subplots()
        ax.grid()
        ax.plot(x,y,color='black')
        ax.set_aspect('equal')

