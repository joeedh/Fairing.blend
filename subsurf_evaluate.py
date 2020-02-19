import bpy, bmesh, random, time, struct
from math import *
from mathutils import *

import bpy, bmesh, time, random, os, os.path
from ctypes import *
from mathutils import *
from math import *
from .sym import sym
from .bezpatch import BezPatch
from .subsurf import subdivide2, subdivide
from . import globals

tn = 0
def getText(text="", co=[0,0,0], scale=0.1, mat="TextMaterial", only_check=False):
    global tn
    tn += 1
    
    scene = bpy.context.scene
    
    name = "_text" + str(tn)
    
    if name not in bpy.data.curves and only_check:
        return None
        
    if name not in bpy.data.curves:
        tx = bpy.data.curves.new(name, type="FONT")
    else:
        tx = bpy.data.curves[name]
    
    if name not in bpy.data.objects:
        ob = bpy.data.objects.new(name, object_data=tx)
        scene.collection.objects.link(ob)
    else:
        ob = bpy.data.objects[name]
        try: #try to add to scene, sanity check
            scene.collection.objects.link(ob)
        except:
            pass
    
    if mat in bpy.data.materials:
        mat = bpy.data.materials[mat]
        if len(tx.materials) == 0:
            tx.materials.append(mat)
        
    tx.body = text
    ob.location = Vector(co)
    ob.scale = Vector([scale, scale, scale])
    ob.update_tag()
    
    return ob

for i in range(32):
    ob = getText("", only_check=True)
    if ob != None:
        ob.location = Vector([100, 100, 100])
tn = 0
    
#"""
#https://svn.blender.org/svnroot/bf-blender/branches/pynodes/extern/qdune/primitives/CCSubdivision.cpp?p=37276
def makeBasis(matData):
    #return
    #make bspline evaluation basis
    buv = [[0 for x in range(16)] for y in range(16)]
    
    #// bspline basis (could use RiBSplineBasis, but want double prec by default)
    bsp =[[-1.0/6.0,     0.5,    -0.5, 1.0/6.0],
          [     0.5,    -1.0,     0.5,     0.0],
          [    -0.5,     0.0,     0.5,     0.0],
          [ 1.0/6.0, 4.0/6.0, 1.0/6.0,     0.0]];
    
    for i in range(16):
        d = i >> 2
        r = i & 3
        
        for v in range(4):
            for u in range(4):
                buv[i][v*4 + u] = bsp[u][d]*bsp[v][r]
    
    #max size needed for N==50
    tmp = [0 for x in range(1728*2)]

    eigen = matData[3:]
    #return
    for rn in range(len(eigen)-2):
        K = 2*(rn + 3) + 8
        
        for k in range(3):
            for k2 in range(K*16):
                tmp[k2] = 0
            
            idx = 0
            for i in range(K):
                for j in range(16):
                    sc = eigen[rn].phi[k][i + j*K]
                    
                    y4 = 0
                    while y4 < 16:
                        for x in range(4):
                            tmp[idx + y4 + x] += sc*buv[j][y4 + x]
                        y4 += 4
                
                idx += 16
            
            #now replace 'Phi' by tmp array
            for i in range(K*16):
                eigen[rn].phi[k][i] = tmp[i] 
            #memcpy(const_cast<double*>(&eigen[rn].x[k][0]), tmp, sizeof(double)*K*16);

               
def EvalSpline_du(C, i, u, v):
    c0 = C[i+0]; c1 = C[i+1]; c2 = C[i+2]; c3 = C[i+3]; c4 = C[i+4]; c5 = C[i+5];
    c6 = C[i+6]; c7 = C[i+7]; c8 = C[i+8]; c9 = C[i+9]; c10 = C[i+10]; c11 = C[i+11];
    c12 = C[i+12]; c13 = C[i+13]; c14 = C[i+14]; c15 = C[i+15]
    
    return ((((c0*v+c1)*v+c2)*v+c3)*u+((c4*v+c5)*v+c6)*v+c7)*u+((c8*v+ \
          c9)*v+c10)*v+c11+((((c0*v+c1)*v+c2)*v+c3)*u+((c4*v+c5)*v+c6)*v \
          +c7+(((c0*v+c1)*v+c2)*v+c3)*u)*u;
          
def EvalSpline_dv(C, i, u, v):
    c0 = C[i+0]; c1 = C[i+1]; c2 = C[i+2]; c3 = C[i+3]; c4 = C[i+4]; c5 = C[i+5];
    c6 = C[i+6]; c7 = C[i+7]; c8 = C[i+8]; c9 = C[i+9]; c10 = C[i+10]; c11 = C[i+11];
    c12 = C[i+12]; c13 = C[i+13]; c14 = C[i+14]; c15 = C[i+15]
    
    return (((3*c0*v**2+2*c1*v+c2)*u+(2*c4*v+c5)*v+(c4*v+c5)*v+c6)*u+ \
            c10+3*c8*v**2+2*c9*v)*u+3*c12*v**2+2*c13*v+c14

def EvalSpline_du2(C, i, u, v):
    c0 = C[i+0]; c1 = C[i+1]; c2 = C[i+2]; c3 = C[i+3]; c4 = C[i+4]; c5 = C[i+5];
    c6 = C[i+6]; c7 = C[i+7]; c8 = C[i+8]; c9 = C[i+9]; c10 = C[i+10]; c11 = C[i+11];
    c12 = C[i+12]; c13 = C[i+13]; c14 = C[i+14]; c15 = C[i+15]
    
    return 2*(3*c0*u*v**3+3*c1*u*v**2+3*c2*u*v+3*c3*u+c4*v**3+c5*v**2+c6*v+c7);
          
def EvalSpline_dv2(C, i, u, v):
    c0 = C[i+0]; c1 = C[i+1]; c2 = C[i+2]; c3 = C[i+3]; c4 = C[i+4]; c5 = C[i+5];
    c6 = C[i+6]; c7 = C[i+7]; c8 = C[i+8]; c9 = C[i+9]; c10 = C[i+10]; c11 = C[i+11];
    c12 = C[i+12]; c13 = C[i+13]; c14 = C[i+14]; c15 = C[i+15]
    
    return 2*((((3*c0*v+c1)*u+3*c4*v+c5)*u+3*c8*v+c9)*u+3*c12*v+c13);
     
def EvalSpline(C, i, u, v):
    u = max(u, 0.000001)
    v = max(v, 0.000001)
    
    """
    f := ((((((C0*v + C1)*v + C2)*v + C3) *u +
            (((C4*v + C5)*v + C6)*v + C7))*u +
            (((C8*v + C9)*v + C10)*v + C11))*u +
            (((C12*v + C13)*v + C14)*v + C15));
    df(f, u);
    df(f, v);
    
    """
    #"""
    return ((((((C[i+ 0]*v + C[i+ 1])*v + C[i+ 2])*v + C[i+ 3]) *u + \
            (((C[i+ 4]*v + C[i+ 5])*v + C[i+ 6])*v + C[i+ 7]))*u +\
            (((C[i+ 8]*v + C[i+ 9])*v + C[i+10])*v + C[i+11]))*u +\
            (((C[i+12]*v + C[i+13])*v + C[i+14])*v + C[i+15]));
    #"""
    sO = C[ i+0] + u*(C[ i+1] + u*(C[ i+2] + u*C[ i+3])) 
    s1 = C[ i+4] + u*(C[ i+5] + u*(C[ i+6] + u*C[ i+7])) 
    s2 = C[ i+8] + u*(C[ i+9] + u*(C[i+10] + u*C[i+11])) 
    s3 = C[i+12] + u*(C[i+13] + u*(C[i+14] + u*C[i+15]))
    return sO + v* (s1 + v* (s2 + v*s3)) 

def EvalSpline3(C, u, v):
    u = max(u, 0.000001)
    v = max(v, 0.000001)
    
    ret = Vector()
    for i in range(3):
        #"""
        ret[i] = ((((((C[  0][i]*v + C[  1][i])*v + C[  2][i])*v + C[  3][i]) *u + \
                (((C[  4][i]*v + C[  5][i])*v + C[  6][i])*v + C[  7][i]))*u +\
                (((C[  8][i]*v + C[  9][i])*v + C[ 10][i])*v + C[ 11][i]))*u +\
                (((C[ 12][i]*v + C[ 13][i])*v + C[ 14][i])*v + C[ 15][i]));
        #"""
        continue
        sO = C[  0][i] + u*(C[  1][i] + u*(C[  2][i] + u*C[  3][i])) 
        s1 = C[  4][i] + u*(C[  5][i] + u*(C[  6][i] + u*C[  7][i])) 
        s2 = C[  8][i] + u*(C[  9][i] + u*(C[ 10][i] + u*C[ 11][i])) 
        s3 = C[ 12][i] + u*(C[ 13][i] + u*(C[ 14][i] + u*C[ 15][i]))
        ret[i] = sO + v* (s1 + v* (s2 + v*s3)) 
    
    return ret

class Disk:
    def __init__(self, vs, es, fs):
        self.fs = fs
        self.es = es
        self.vs = vs

def find_loop(v, f):
    l = None
    for l2 in v.link_loops:
        if l2.vert == v and l2.face == f:
            l = l2
            break

    if l is None:
        print("failed", len(v.link_edges), f)
        raise RuntimeError("find_loop() failed")
        return None
    
    return l
    
def getdisk(v, f=None):
    if f is None:
        for f2 in v.link_faces:
            f = f2
            break
            
    val = len(v.link_edges)
    l = find_loop(v, f);
    
    startl = l
    
    vs = []
    es = []
    fs = []
    
    ret = Disk(vs, es, fs)
    
    first = l
    _i = 0
    while _i < 1000:
        l = l.link_loop_next;
        vs.append(l.vert)
        
        l = l.link_loop_next
        vs.append(l.vert)
        
        l = l.link_loop_next
        l = l.link_loop_radial_next;
        
        if l == first:
            break
            
        _i += 1
    
    vs.reverse()
    es.reverse()
    fs.reverse()
    
    offs = [
        0,
        0,
        0,
        1,
        3,
        5,
        7,
        9,
        0,
        0,
        0,
    ]
    
    vs2 = [startl.vert]
    es2 = [startl.vert]
    fs2 = [startl.vert]
    for i in range(len(vs)):
        i2 = (i + (2*(val-2)-1)) % len(vs)
        vs2.append(vs[i2])
        es2.append(vs[i2])
        fs2.append(vs[i2])
    vs[:] = vs2
    
    return ret

def opposite_quad(v, f=None):
    for l in f.loops:
        if l.vert == v:
            return l.link_loop_next.link_loop_next
            
    return None
    
def edgeloop_next_old(l):
    l = l.link_loop_radial_next.link_loop_next
    l = l.link_loop_radial_next
    return l
    
def edgeloop_next(l):
    l = l.link_loop_next.link_loop_radial_next.link_loop_next
    
    return l
    
def edgeloop_prev(l):
    l = l.link_loop_prev.link_loop_radial_prev.link_loop_prev
    
    #l = l.link_loop_radial_prev
    return l
    
def getPatchVs4(f):
    vs = []
    fl = list(f.loops)
    
    #vs = list(f.verts)
   
    l = edgeloop_prev(fl[1]).link_loop_prev

    l = edgeloop_prev(l)
    vs.append(l.vert)
    l = edgeloop_next(l)
    vs.append(l.vert)
    l = edgeloop_next(l)
    vs.append(l.vert)
    l = edgeloop_next(l)
    vs.append(l.vert)
    
    vs.append(edgeloop_prev(fl[0]).vert)
    vs.append(fl[0].vert)
    vs.append(fl[1].vert)
    vs.append(edgeloop_next(fl[0]).link_loop_next.vert)
    
    l8 = l = edgeloop_next(edgeloop_next(fl[2]))
    vs.append(l.vert)
    l = edgeloop_prev(l)
    vs.append(l.vert)
    l = edgeloop_prev(l)
    vs.append(l.vert)
    l = edgeloop_prev(l)
    vs.append(l.vert)
    
    l8 = l8.link_loop_next;
    
    l = edgeloop_prev(l8).link_loop_prev
    vs.append(l.vert)
    l = edgeloop_prev(l)
    vs.append(l.vert)
    l = edgeloop_prev(l)
    vs.append(l.vert)
    l = edgeloop_prev(l)
    vs.append(l.vert)
    
    return vs
    
visv = 13
visf = 79

def getPatchVs(v, f):
    global visv, visf
    
    disk = getdisk(v, f)
    vs = disk.vs
    
    """
    return list(f.verts)
    #"""
    
    v2 = None
    for l in f.loops:
        if l.vert == v:
            l2 = l.link_loop_next.link_loop_radial_next.link_loop_next
            l2 = l2.link_loop_prev.link_loop_prev.link_loop_radial_next
            l2 = l2.link_loop_next.link_loop_next
            
            #l2 = l2.link_loop_radial_next.link_loop_prev.link_loop_radial_next
            #l2 = l2.link_loop_next
            l = l2
            break
            
    if l is None:
        print("error l was undefined", v2, f)
        return vs
    if l.vert in vs:
        #print("error wrong vert", v2, f)
        #return vs
        pass
    
    vs.append(l.vert)

    vs.append(l.link_loop_prev.vert)
    l2 = l.link_loop_prev.link_loop_prev.link_loop_radial_next
    l2 = l2.link_loop_prev
    vs.append(l2.vert)
    l2 = l2.link_loop_prev.link_loop_radial_next
    l2 = l2.link_loop_prev
    vs.append(l2.vert)

    vs.append(l.link_loop_next.vert)
    l2 = edgeloop_next_old(l.link_loop_next)
    vs.append(l2.vert)
    l2 = edgeloop_next_old(l2).link_loop_next.link_loop_next
    vs.append(l2.vert)
    

    #if v.index == visv and f.index == visf:
    #    l2.edge.hide = 1
    
    if v.index == visv and f.index == visf:
        for v2 in vs:
            v2.hide = 1
            pass
        
    return vs
    
    
#// returns integer log2(1.0/x)
def invilog2_roundup(x):
    if (x==0.0): return 33
    xi = int(1.0/x)
    lg2 = 0
    
    while xi > 0:
        lg2 += 2
        xi = xi>>1
    
    return lg2

def mypow(x, y):
    if x < 0 or ((x == 0.0 and y == 0.0)):
        return 1.0
    return pow(x, y)
    
class Patch:
    def __init__(self, bmv, bmf, n):
        global visv, visf
        global matData
        
        self.n = len(bmv.link_edges)
        self.v = bmv
        self.f = bmf
        
        if self.n != 4:
            self.ev = matData[n]
            
            self.__ = bmv.co * 1
            
            self.control = getPatchVs(bmv, bmf)
            self.proj = [Vector() for x in range(len(self.control))]
            
            if bmv.index == visv and bmf.index ==  visf:
                for i, v2 in enumerate(self.control):
                    no = Vector(v2.normal)
                    no.normalize()
                    no *= 0.025;
                    getText(str(i), v2.co + no)
        else:
            def ci(x, y):
                return y*4 + x
            
            self.control = [None for v in range(16)]
            cs = self.control
            
            vs2 = getPatchVs4(bmf)
            for i, v in enumerate(vs2):
                cs[i] = Vector(v.co)
                
            if bmf.index == visf:
                for i, v in enumerate(vs2):
                    getText(str(i), v.co)
            
                
            self.patch = BezPatch()
            vi = 0
            
            for i in range(4):
                for j in range(4):
                    self.patch[j][i] = cs[vi]
                    vi += 1
            
    def evaluate_all(self, u, v):
        if self.n != 4:
            ret = [
                None,
                Vector(), #du
                Vector(), #dv
                Vector(), #du2
                Vector(), #dv2
            ]
            
            ret[0] = self.evaluate_ext(u, v, ret[1], ret[2], ret[3], ret[4])
            return ret
        else:
            return [self.evaluate(u, v)] + self.calcDvs(u, v)
            
    def evaluate(self, u, v):
        if self.n == 4:
            return self.patch.evaluate(u, v)
        else:
            return self.evaluate_ext(u, v)
    
    def calcDvs(self, u, v):
        if self.n == 4:
            return [
                self.patch.derivative_u(u, v),
                self.patch.derivative_v(u, v),
                self.patch.derivative_u2(u, v),
                self.patch.derivative_v2(u, v)
            ]
        else:
            ret = [
                Vector(), #du
                Vector(), #dv
                Vector(), #du2
                Vector(), #dv2
            ]
            
            self.evaluate_ext(u, v, ret[0], ret[1], ret[2], ret[3])
            return ret
    
    def evaluate_ext(self, u, v, dv_u=None, dv_v=None, dv_u2=None, dv_v2=None):
        if u == 0.0:
            u = 0.0000001
        if v == 0.0:
            v = 0.0000001
        
        """
        p = Vector()
        p[0] = u
        p[1] = v
        p[2] = 1
        return p
        #"""
        
        ev = self.ev
        ev.projectPoints(self.control, self.proj)

        cn = self.control
        cp = self.proj
        
        #print(len(self.control), len(self.proj), self.ev.k)
        #print(u, v)
        
        logu = -log(u, 2.0)
        logv = -log(v, 2.0)
        
        #logu = 0 if u == 0 else -log(u, 2)
        #logv = 0 if v == 0 else -log(v, 2)
        
        n = floor(min(logu, logv)+1)
        #n = invilog2_roundup(max(u, v))
        #print(u, v, n, logu, logv, invilog2_roundup(max(u, v)))
        
        pow2 = pow(2.0, n-1);
        u *= pow2; v *= pow2
        
        if v < 0.5:
            k=0; u=2*u-1; v=2*v
        elif u < 0.5:
            k=2; u=2*u; v=2*v-1
        else:
            k=1; u=2*u-1; v=2*v-1
        
        p = Vector()
        wsum = 0
        
        if dv_u is not None:
            dv_u[0] = dv_u[1] = dv_u[2] = 0.0
        if dv_v is not None:
            dv_v[0] = dv_v[1] = dv_v[2] = 0.0
            
        for i in range(0, 2*self.ev.n + 8): #self.ev.k):
            p2 = EvalSpline(ev.phi[k], i*16, u, v)
                
            if ev.val[i] == 0.0:
                ev.val[i] = 0.0000001
                pass
                continue
                
            w = pow(ev.val[i], n-1)
            
            dw = 1#pow(ev.val[i], n-1)
            #for higher-order derives: pow(2, n*dvorder)
            dvw1 = 2#pow(2, n); 
            dvw2 = 8#pow(2, 2*n);
            
            if dv_u:
                du = EvalSpline_du(ev.phi[k], i*16, u, v)*dw*dvw1
                du2 = EvalSpline_du2(ev.phi[k], i*16, u, v)*dw*dvw2
            if dv_v:
                dv = EvalSpline_dv(ev.phi[k], i*16, u, v)*dw*dvw1
                dv2 = EvalSpline_dv2(ev.phi[k], i*16, u, v)*dw*dvw2
            
            #print(w)
            #print(p2[k])
            
            for j in range(3):
                #p[j] += cn[i].co[j]*w
                p[j] += cp[i][j]*p2*w
                wsum += w
                
                if dv_u:
                    dv_u[j] += du*cp[i][j]
                    dv_u2[j] += du2*cp[i][j]
                if dv_v:
                    dv_v[j] += dv*cp[i][j]
                    dv_v2[j] += dv2*cp[i][j]
        
        #print(wsum, n, self.ev.k, len(self.control))
        
        #wsum=0        
        #if wsum != 0:
        #    p *= 1.0 / wsum
        
        return p
        
class EVALSTRUCT (Structure):
    _fields_ : [
        ("val" , POINTER(c_double)),
        ("vecI" , POINTER(c_double)),
        ("phi" , POINTER(POINTER(c_double))),
    ]
    
    n = None #valence
    ps = None #internal patches
        
    def projectPoints(self, pin, pout):
        n = self.n
        K = self.k #2*n + 8
        cp = pout

        for i in range(K):
            cp[i][0] = 0;
            cp[i][1] = 0;
            cp[i][2] = 0;
            
            inv = self.vecI
            
            for j in range(K):
                idx = j*K + i
                for m in range(3):
                    try:
                      cp[i][m] += inv[idx] * pin[j].co[m]
                    except ReferenceError:
                      if random.random() < 0.00002: 
                        print("Missing verts", pin[j]);
                    

    #def evaluate(
matData = None

def getMatData(path=None):
    global matData
    
    if matData is None:
        if hasattr(globals, "matData"):
            matData = globals.matData
        else:
            loadMatData(path)
            globals.matData = matData
            
    return matData
    

def loadMatData(path=None):
    global matData
    print("loading subsurf matrix data. . .")
    
    if path is None:
        basepath = os.path.split(bpy.data.filepath)[0]
        path = basepath + "/subdiv_matrices/ccdata50NT.dat"
    
    file = open(path, "rb");
    buf = [file.read()]
    file.close()
    _i = [0]
    
    def readint():
        s = buf[0][:4]
        buf[0] = buf[0][4:]
        return struct.unpack("<i", s)[0]
        
    def readdoubles(out, n):
        memmove(out, cadd + _i[0], n*8)
        _i[0] += n*8
      
    
    #don't include 4-byte size argument in buf2, that'll mess up 8-byte alignment
    nmax = readint()

    buf2 = bytearray(buf[0]);
    
    cbuf = (c_char*(len(buf2)))()
    cbuf[:] = buf2
    cadd = addressof(cbuf)
    
    #print(cbuf)
    #return
    evs = [None, None, None];
    
    for i in range(nmax-2):
      N = i+3;
      K = 2*N+8;
      
      ev = EVALSTRUCT()
      ev.n = N
      ev.k = K
      
      evs.append(ev)
      
      ev.val = (c_double*K)()
      ev.vecI = (c_double*(K*K))()
      ev.phi = (POINTER(c_double)*3)()
      
      #print("i", i+1, "of", nmax-2)
      #print(K, N)
      
      ev.phi[0] = (c_double*(K*16))();
      ev.phi[1] = (c_double*(K*16))();
      ev.phi[2] = (c_double*(K*16))();
    
      readdoubles(ev.val, K)
      readdoubles(ev.vecI, K*K)
      readdoubles(ev.phi[0], K*16)
      readdoubles(ev.phi[1], K*16)
      readdoubles(ev.phi[2], K*16)
      
    matData = evs
    makeBasis(matData)
      
    print("done", _i[0], "of", len(buf2));
    
    return evs;
    
class SSVert:
    def __init__(self, bmv):
        self.v = bmv
        self.patches = []
    
    #will align u to go along spoke edge 
    def calcDvs(self, p):
        if p.n == 4:
            vi = 0
            
            uvs = [
              [0, 0],
              [1, 0],
              [1, 1],
              [0, 1]
            ]
            
            for v in p.f.verts:
                if v == self.v:
                    break
                vi += 1
            
            flip = vi & 1
            
            if flip:
                v, u = uvs[vi]
            else:
                u, v = uvs[vi]
            
            return p.calcDvs(u, v)
        elif p.v == self.v:
            return p.calcDvs(0.000001, 0.00001)
        
def copyBM(bm):
  #bm2 = bm.copy() #doesn't work!!!
  bm2 = bmesh.new()
  
  bm2.verts.index_update()
  vs = []
  
  for v in bm.verts:
    v2 = bm2.verts.new(v.co)
    v2.normal = v.normal
    v2.select = v.select
    v2.hide = v.hide
    vs.append(v2)
    
  for e in bm.edges:
    v1 = vs[e.verts[0].index]
    v2 = vs[e.verts[1].index]
    e2 = bm2.edges.new([v1, v2])
    e2.select = e.select
    e2.hide = e.hide
    
  bm.verts.index_update()
  
  for f in bm.faces:
    vs2 = [vs[v.index] for v in f.verts]
    f2 = bm2.faces.new(vs2)
    
    f2.select = f.select;
    f2.hide = f.hide
    f2.normal = f.normal
    
  print("sadasadasdasd", len(bm2.verts))
  
  bm2.verts.index_update()
  bm2.edges.index_update()
  bm2.faces.index_update()
  
  return bm2
  
class SSMesh:
    def __init__(self):
        self.patches = []   
        self.verts = []
        self.vset = set()
        self.pmap = {}
        
    def load(self, bm):
        getMatData()
        
        #bm = copyBM(bm);
        
        self.verts = [SSVert(v) for v in bm.verts]
        
        subdivide(bm)
        subdivide(bm)
        
        #subdivide2(bm)
        bm.verts.index_update()
        bm.edges.index_update()
        bm.faces.index_update()
        
        fset = set()
        pmap = self.pmap
        
        for v in bm.verts:
            es = list(v.link_edges)
            if len(es) < 3 or len(es) == 4: continue
            
            for f in v.link_faces:
                if f in fset: continue
                
                fset.add(f)

                p = Patch(v, f, len(es))
                self.patches.append(p)
                pmap[f] = p
        
        for f in bm.faces:
            if f in fset: continue
            fv = None
            
            fset.add(f)
            
            for v in f.verts:
                if len(v.link_edges) == 4:
                    fv = v
                    break;
                    
                    
            p = Patch(fv, f, 4)
            self.patches.append(p)
            pmap[f] = p
        
        
        for v in self.verts:
            for f in v.v.link_faces:
                v.patches.append(pmap[f])
            
    def render_patch(self, p, bm, steps=7, render_dvs=False):
        def edge(v1, v2):
            if type(v1) == Vector:
                v1 = bm.verts.new(v1)
            if type(v2) == Vector:
                v2 = bm.verts.new(v2)
            return bm.edges.new([v1, v2])
            
        ds = 0.999999 / (steps-1)
        grid = [[0 for x in range(steps)] for y in range(steps)]
        
        #steps=3
        
        u = 0
        for i in range(steps):
            v = 0.0
            for j in range(steps):
                if render_dvs:
                    ret = p.evaluate_all(u, v)
                    co = ret[0]
                    
                    if len(ret) > 3:
                        no = Vector(ret[1]).cross(ret[2])
                        #print(ret[1]-ret[2])
                        no.normalize()
                        
                        du1 = ret[1]
                        du2 = ret[3]
                        dv1 = ret[2]
                        dv2 = ret[4]
                        
                        ku = du2.length / du1.length
                        #ku *= -1 if du2.dot(no) >= 0 else 1

                        kv = dv2.length / dv1.length
                        #kv *= -1 if dv2.dot(no) >= 0 else 1
                        
                        tan = dv1
                        #tan.normalize();
                        
                        #edge(co, co + no*ku);
                        
                        co2 = co + no*0.5*kv
                        edge(co, co2)
                        edge(co2, co2 + tan*0.25)
                        #edge(co, co+ret[1]*0.25)
                        #edge(co, co+ret[2]*0.25)
                else:
                    co = p.evaluate(u, v)
                grid[i][j] = bm.verts.new(co)
                v += ds
            u += ds
        
        for i in range(steps-1):
            for j in range(steps-1):
                v1 = grid[i][j]; v2 = grid[i][j+1]
                v3 = grid[i+1][j+1]; v4 = grid[i+1][j]
                bm.faces.new([v4, v3, v2, v1])
                #bm.edges.new([v1, v2])
                #bm.edges.new([v2, v3])
        
        #bm.update_normals()
        #bm.normals_update()
    
    def renderVal4Cages(self, bm, selected_only=False):
        for p in self.patches:
            if selected_only and not p.f.select: continue
            
            if p.n == 4:
                p.patch.renderCage(bm)
    
    def render(self, bm, steps=7, render_dvs=False, selected_only=False):
        tot4 = 0
        for p in self.patches:
            if selected_only and not p.f.select:
                continue
                
            self.render_patch(p, bm, steps, render_dvs)
        
    
    def renderDerivatives(self, bm, steps=4):
        def edge(v1, v2):
            if type(v1) == Vector:
                v1 = bm.verts.new(v1)
            if type(v2) == Vector:
                v2 = bm.verts.new(v2)
            return bm.edges.new([v1, v2])
        
        for v in self.verts:
            for p in v.patches:
                dvs = v.calcDvs(p)
                if dvs is None:
                    print("dvs was None!")
                    continue
                
                edge(v.v.co, v.v.co+dvs[0])
                no = dvs[0].cross(dvs[1])
                no.normalize()
                
                if len(dvs) <= 2:
                    continue
                
                k1 = dvs[2].length / dvs[0].length
                k2 = dvs[3].length / dvs[1].length
                
                edge(v.v.co, v.v.co+no*k1)
                edge(v.v.co, v.v.co+no*k2)
                
                #edge(v.v.co, v.v.co+dvs[1])
                #if len(dvs) > 2:
                #    edge(v.v.co, v.v.co-dvs[2])
                
                
                
        
        