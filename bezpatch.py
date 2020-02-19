import bpy, bmesh, random, time
from mathutils import *
from math import *

from .sym import sym

def make_equation(n):
    def basis(i, s):
        s *= n-2
        tdiv = 1.0# / (2*n)
        
        #uniform cox de boor
        def deboor(i, n):
            i2 = i-3
            
            ti = i2
            ti2 = i2+1
            tip = i2+n
            tip2 = i2+n+1
            
            ti *= tdiv
            ti2 *= tdiv
            tip *= tdiv
            tip2 *= tdiv
            
            if n == 0:
                return sym("step("+str(s)+", "+str(ti)+","+str(ti2)+")")
            else:
                a = (s-ti) / (tip - ti)
                b = (tip2 - s) / (tip2 - ti2)
                
                return deboor(i, n-1)*a + deboor(i + 1, n-1)*b;
        
        return deboor(i, n)
    
    u = sym("u")
    v = sym("v")
    p = sym(0)
    
    for i in range(n+1):
        for j in range(n+1):
            w = basis(i, u)*basis(j, v);
            idx = j*4 + i
            p += sym("c"+str(j)+"") * w
    
    return str(p)
    
binomials = None
def make_binomials(max_bin=11):
    global binomials
    
    binomials = [[0 for y in range(x+1)] for x in range(max_bin+1)]
    
    def factorial(n):
        f = 1
        
        for i in range(2, n+1):
            f *= i
        return f
        
    def falling_fac(n, k):
        n2 = 1
        for i in range(k):
            n2 *= n-i
            
        return n2
        
    for i in range(0, max_bin):
        for j in range(0, i+1):
            print(i, j, binomials)
            binomials[i][j] = falling_fac(i, j) / factorial(j)
            #binomials[i][j] = falling_fac(i, j)
    
    return binomials

from . import globals
if not hasattr(globals, "binomials") or globals.binomials is None:
    globals.binomials = make_binomials()
else:
    binomials = globals.binomials
    
class BezPatch (list):
    def __init__(self, patch=None, order=4):
        self[:] = [[Vector() for v in range(order+1)] for u in range(order+1)]
        
        if patch is not None:
            for i in range(len(patch)): 
                for j in range(len(patch[i])):
                    self[i][j] = Vector(patch[i][j])
                    
        self.order = order
        self.degree = order-1
        self.face = None
        
    def renderCage(self, bm):
        grid = [[0 for x in range(self.order)] for y in range(self.order)]
        
        for i in range(self.order):
            for j in range(self.order):
                grid[i][j] = bm.verts.new(self[i][j])
        
        for i in range(self.order-1):
            for j in range(self.order-1):
                bm.faces.new([grid[i][j], grid[i][j+1], grid[i+1][j+1], grid[i+1][j]])
                
        return grid
    
    def render(self, bm, steps=24):
        ds = 1.0 / (steps - 1)
        
        grid = [[0 for x in range(steps)] for y in range(steps)]
        
        u = 0
        for i in range(steps):
            v = 0
            
            for j in range(steps):
                p = self.evaluate(u, v)
                bv = bm.verts.new(p)
                
                grid[i][j] = bv
                v += ds
            u += ds
         
        for i in range(steps-1):
            for j in range(steps-1):
                bm.faces.new([grid[i][j], grid[i][j+1], grid[i+1][j+1], grid[i+1][j]])

    def fromCage(self, cage):
        for i in range(self.order):
            for j in range(self.order):
                self[i][j] = Vector(cage[i][j].co)
                
        return self
    
    def basis(self, i, s):
        #return self.basis_bernstein(i, s)
        n = self.degree
        
        s *= n-2
        tdiv = 1.0# / (2*n)
        
        #uniform cox de boor
        def deboor(i, n):
            i2 = i-3
            
            ti = i2
            ti2 = i2+1
            tip = i2+n
            tip2 = i2+n+1
            
            ti *= tdiv
            ti2 *= tdiv
            tip *= tdiv
            tip2 *= tdiv
            
            if n == 0:
                return 1 if s >= ti and s < ti2 else 0
            else:
                a = (s-ti) / (tip - ti)
                b = (tip2 - s) / (tip2 - ti2)
                
                return deboor(i, n-1)*a + deboor(i + 1, n-1)*b;
        
        return deboor(i, n)
    
    def basis_bernstein(self, i, s):
        n = self.degree
        return binomials[n][i] * s**i * (1.0 - s)**(n - i)
    
    def normal(self, u, v):
        df = 0.0001
        p = self.evaluate(u, v)
        tanu = (self.evaluate(u+df, v) - p) / df
        tanv = (self.evaluate(u, v+df) - p) / df
        
        #if we were perfectly arc-parameterized, then we could just
        #take the second derivative finite difference, but oh well
        
        n = tanv.cross(tanu)
        n.normalize()
        return n;
    
    def derivative_u(self, u, v):
        df = 0.001
        
        if u > df*2 and u < 1.0 - df*2:
            a = self.evaluate(u-df, v)
            b = self.evaluate(u+df, v)
        elif u > 1.0 - df:
            a = self.evaluate(u-df, v)
            b = self.evaluate(u, v)
        else:
            a = self.evaluate(u, v)
            b = self.evaluate(u+df, v)
        
        return (b - a) / df
        
    def derivative_v(self, u, v):
        df = 0.001
        
        if v > df*2 and v < 1.0 - df*2:
            a = self.evaluate(u, v-df)
            b = self.evaluate(u, v+df)
        elif v > 1.0 - df:
            a = self.evaluate(u, v-df)
            b = self.evaluate(u, v)
        else:
            a = self.evaluate(u, v)
            b = self.evaluate(u, v+df)
        
        return (b - a) / df
        
    def derivative_u2(self, u, v):
        df = 0.001
        
        if u > df*2 and u < 1.0 - df*2:
            a = self.derivative_u(u-df, v)
            b = self.derivative_u(u+df, v)
        elif u > 1.0 - df:
            a = self.derivative_u(u-df, v)
            b = self.derivative_u(u, v)
        else:
            a = self.derivative_u(u, v)
            b = self.derivative_u(u+df, v)
        
        return (b - a) / df
        
    def derivative_v2(self, u, v):
        df = 0.001
        
        if v > df*2 and v < 1.0 - df*2:
            a = self.derivative_v(u, v-df)
            b = self.derivative_v(u, v+df)
        elif v > 1.0 - df:
            a = self.derivative_v(u, v-df)
            b = self.derivative_v(u, v)
        else:
            a = self.derivative_v(u, v)
            b = self.derivative_v(u, v+df)
        
        return (b - a) / df
        
    def evaluate(self, u, v):
        p = Vector()
        
        def step(s, a, b):
            return 1 if s >= a and s < b else 0
        
        def evalspline(axis):
            c0 = self[0][0][axis]; c1 = self[0][1][axis]; c2 = self[0][2][axis]; c3 = self[0][3][axis];
            c4 = self[1][0][axis]; c5 = self[1][1][axis]; c6 = self[1][2][axis]; c7 = self[1][3][axis];
            c8 = self[2][0][axis]; c1 = self[2][1][axis]; c9 = self[2][2][axis]; c10 = self[2][3][axis];
            c12 = self[3][0][axis]; c13 = self[3][1][axis]; c14 = self[3][2][axis]; c15 = self[3][3][axis];
            
            ans3=(step(u,3,4)*u**3-12*step(u,3,4)*u**2+48*step(u,3,4)*u-64*
            step(u,3,4)-2*step(u,2,3)*u**3+15*step(u,2,3)*u**2-33*step(u,2
            ,3)*u+17*step(u,2,3)+step(u,1,2)*u**3-3*step(u,1,2)*u**2+3*
            step(u,1,2)*u-7*step(u,1,2)-6*step(u,0,1)-step(u,(-1),0)*u**3-
            6*step(u,(-1),0)+2*step(u,(-2),(-1))*u**3+9*step(u,(-2),(-1))*
            u**2+9*step(u,(-2),(-1))*u-3*step(u,(-2),(-1))-step(u,(-3),(-2
            ))*u**3-9*step(u,(-3),(-2))*u**2-27*step(u,(-3),(-2))*u-27*
            step(u,(-3),(-2)))
            ans2=(((((v+1)*step(v,(-1),0)-(v-1)*step(v,0,1))*(v+1)+((v-2)*
            step(v,1,2)-step(v,0,1)*v)*(v-2))*(v+1)+(((v-1)*step(v,1,2)-(v
            -3)*step(v,2,3))*(v-3)+((v-2)*step(v,1,2)-step(v,0,1)*v)*v)*(v
            -3))*c2-((((v-1)*step(v,1,2)-(v-3)*step(v,2,3))*(v-1)-((v-2)*
            step(v,2,3)-(v-4)*step(v,3,4))*(v-4))*(v-4)+(((v-1)*step(v,1,2
            )-(v-3)*step(v,2,3))*(v-3)+((v-2)*step(v,1,2)-step(v,0,1)*v)*v
            )*v)*c3+((((v+2)*step(v,(-2),(-1))-step(v,(-1),0)*v)*(v+2)-((v
            +1)*step(v,(-1),0)-(v-1)*step(v,0,1))*(v-1))*(v+2)-(((v+1)*
            step(v,(-1),0)-(v-1)*step(v,0,1))*(v+1)+((v-2)*step(v,1,2)-
            step(v,0,1)*v)*(v-2))*(v-2))*c1+((((v+3)*step(v,(-3),(-2))-(v+
            1)*step(v,(-2),(-1)))*(v+3)-((v+2)*step(v,(-2),(-1))-step(v,(
            -1),0)*v)*v)*(v+3)-(((v+2)*step(v,(-2),(-1))-step(v,(-1),0)*v)
            *(v+2)-((v+1)*step(v,(-1),0)-(v-1)*step(v,0,1))*(v-1))*(v-1))*
            c0)*ans3
            ans1=-ans2
            ans=ans1/36
            return ans
        
        """        
        p = Vector()
        for i in range(3):
            p[i] = evalspline(i)
        return p
        #"""
        #return p
        for i in range(self.order):
            for j in range(self.order):
                w = self.basis(i, u)*self.basis(j, v);
                
                p += self[i][j] * w
        
        return p
    

