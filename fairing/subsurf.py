import bpy, bmesh, random, time
from math import *
from mathutils import *

import bpy, bmesh, time, random
from mathutils import *
from math import *
from .sym import sym
from .smatrix import getSMat, getVertDisk, mulVertsMatrix, \
                     getSplitSMat, getSplitEdgeVertDisk, getBoundFlags
import numpy as np

def edgeloop_next(l):
    l = l.link_loop_next.link_loop_radial_next.link_loop_next
    
    return l
    
def edgeloop_prev(l):
    l = l.link_loop_prev.link_loop_radial_prev.link_loop_prev
    
    #l = l.link_loop_radial_prev
    return l

class SymVector (list):
    def __init__(self, co=[0, 0, 0]):
        self.append(sym(co[0]))
        self.append(sym(co[1]))
        self.append(sym(co[2]))
    
    def copy(self):
        b = SymVector(self)
        return b
        
    def binop(self, b, op):
        ret = SymVector()
        
        if type(b) not in [list, Vector, SymVector]:
            b = SymVector([b, b, b])
            
        for i in range(3):
            ret[i] = sym(self[i], b[i], op)
        return ret
    
    def simp(self):
        self[0] = self[0].simp()
        self[1] = self[1].simp()
        self[2] = self[2].simp()
        
        return self
        
    def __radd__(self, b):
        return SymVector(b).binop(self, "+")
        
    def __iadd__(self, b):
        return self.binop(b, "+")
    def __add__(self, b):
        return self.binop(b, "+")
    def __sub__(self, b):
        return self.binop(b, "-")
    def __mul__(self, b):
        return self.binop(b, "*")
    def __truediv__(self, b):
        return self.binop(b, "/")
    def __pow__(self, b):
        return self.binop(b, "**")
    
    def dot(b):
        return self[0]*b[0] + self[1]*b[1] + self[2]*b[2]
        
    @property
    def length(self):
        return self.dot(self)**0.5
    
    def normalize(self):
        l = self.length
        self[0] /= l;
        self[1] /= l;
        self[2] /= l;
        return self
        
class SymElem:
    def __init__(self, orig, sm):
        self.sm = sm
        self.orig = orig
        self.index = -1
        
class SymVert (SymElem):
    def __init__(self, orig, sm):
        SymElem.__init__(self, orig, sm)
        self.co = SymVector(orig.co)
    
    @property
    def link_edges(self):
        def generator():
            for e in self.orig.link_edges:
                yield self.sm.get(e)
        
        return generator()
        
    @property
    def link_faces(self):
        def generator():
            for f in self.orig.link_faces:
                yield self.sm.get(f)
        
        return generator()
        
class SymLoop (SymElem):
    def __init__(self, orig, sm):
        SymElem.__init__(self, orig, sm)
        
class SymEdge (SymElem):
    def __init__(self, orig, sm):
        SymElem.__init__(self, orig, sm)
        
class SymFace (SymElem):
    def __init__(self, orig):
        SymElem.__init__(self, orig, sm)

class SymMesh:
    def __init__(self):
        pass
class Tags:
    CREASE = 1,
    DART = 2,
    SMOOTH = 4
    
def catmull_clark(bm):
    bm.verts.index_update()
    bm.edges.index_update()
    bm.faces.index_update()
    
    origvs = set(bm.verts)
    origfs = set(bm.faces)
    origes = set(bm.edges)
    splitvs = set()
    centvs = set()
    splitmap = {}
    
    vi = len(bm.verts)
    vcos = [Vector(v.co) for v in bm.verts]
    
    def has_edge(v1, v2):
        for e in v1.link_edges:
            if e.other_vert(v1) == v2: return True
    
    creases = {}
    tags = {}
    
    for e in bm.edges:
        if len(e.link_faces) < 2:
            creases[e.index] = 1.0
            tags[e.index] = Tags.CREASE
        else:
            tags[e.index] = Tags.SMOOTH
            
    for e in list(bm.edges):
        origv = [e.verts[0], e.verts[1]]
        
        newe, newv = bmesh.utils.edge_split(e, e.verts[0], 0.5) #split_edges(bm, edges=[e])

        newv.index = vi
        vi += 1
        vcos.append(Vector(newv.co))
        
        splitvs.add(newv)
        splitmap[newv] = origv
        
    for f in origfs:
        cent = Vector()
        tot = 0.0
        
        for i, v in enumerate(f.verts):
            if i % 2 == 0:
                cent += vcos[v.index]
                #cent += v.co
                
                tot += 1.0
        cent /= tot
        
        cent = bm.verts.new(cent)
        cent.index = vi
        vi += 1
        vcos.append(Vector(cent.co))
        
        centvs.add(cent)
        
        n = len(f.verts)
        fvs = list(f.verts)
        
        for _i in range(n//2):
            i = _i*2
            v = fvs[i]
            
            v1 = cent
            v2 = fvs[(i-1+n)%n]
            v3 = v
            v4 = fvs[(i+1)%n]
            
            f2 = bm.faces.new([v1, v2, v3, v4])
        
        bm.faces.remove(f)
        
    for v in origvs:
        val = len(v.link_edges)
        vs = getVertDisk(v)
        is_boundary = getBoundFlags(v)
        smat, ismat = getSMat(val, is_boundary)
        
        cos = [0, 0, 0]
        #print(len(vs), val)
        
        vcos2 = [Vector(vcos[v.index]) for v in vs]
            
        #for boundary cases, padd with zero'd vectors
        while len(vcos2) < 2*val + 1:
            vcos2.append(Vector())
        
        #print(smat)
        #smat = smat@smat
        mulVertsMatrix(vcos2, smat)
        
        vs[0].co = vcos2[0]
        
    for v in splitvs:
        is_boundary = len(v.link_faces) != 4
        #print(len(v.link_faces))
        
        #v.co = vcos[v.index]
        co = Vector()
        tot = 0.0
        
        if is_boundary:
            v1, v2 = splitmap[v]
            
            v.co = vcos[v1.index]*0.5 + vcos[v2.index]*0.5
            #v.co = v1.co*0.5 + v2.co*0.5
            continue
        
        #can't use matrix here, since we might be at top level
        e1 = None
        e2 = None
        w0 = 0
        if is_boundary:
            #continue
            for e in v.link_edges:
                if e1 is None:
                    e1 = e
                else:
                    e2 = e

            w0 = 0.0
            
            #v1, v2 = splitmap[v]
            #v.co = (v1.co + v2.co)*0.5
            
        dot = 1
        if e1 is not None and e2 is not None:
            t1 = e1.verts[1].co - e1.verts[0].co
            t2 = e2.verts[1].co - e2.verts[0].co
            
            if e1.verts[0] != v:
                t1 = -t1
            if e2.verts[0] != v:
                t2 = -t2
                
            t1.normalize()
            t2.normalize()
            
            if is_boundary:
                dot = t1.dot(t2)
                w0 = 0#*dot
        co += v.co*w0
        tot += w0
            
        fval = len(v.link_faces)
        for e in v.link_edges:
            v2 = e.other_vert(v)
            w = 3/8
            
            if 1: #not is_boundary:
                co += vcos[v2.index]*w
            else:
                co += v2.co*w
            tot += w
            
            for f in e.link_faces:
                cent = Vector()
                tot2 = 0.0
                for v3 in f.verts:
                    if 1: #not is_boundary:
                        cent += vcos[v3.index]
                    else:
                        cent += v3.co
                    tot2 += 1.0
                cent /= tot2
                
                if is_boundary:
                    w = 1/16
                else:
                    w = 100/16
                
                co += cent*w
                tot += w
                
        co /= tot
        v.co = co
        
        #can't use matrix here, since we might be at top level
        continue
        ov1, ov2 = splitmap[v]
        #vs = getVertDisk(v)
        vs = getSplitEdgeVertDisk(v, ov1, ov2)
        smat, ismat = getSplitSMat()

        #XXX deal with boundary cases
        while len(vs) < len(smat):
            vs.append(v)
        
        vcos2 = [Vector(vcos[v2.index]) for v2 in vs]
        
        mulVertsMatrix(vcos2, smat)
            
        vs[0].co = vcos2[0]
        
    bm.normal_update()
    
def subdivide2(bm, scos=None, use_scos=True, first=False, ws=None):
    bm.verts.index_update()
    bm.edges.index_update()
    bm.faces.index_update()
    
    w1 = 3
    w2 = 2
    w3 = 0.2
    w6 = 0.5
    w8 = 0.58
    w9 = -0.1
    
    #ws = [3, 2, 0.2, 0.5, 0.58, -0.1]
    if ws is not None:
       w1, w2, w3, w6, w8, w9 = ws[:7]
        
    def has_edge(v1, v2):
        for e in v1.link_edges:
            if e.other_vert(v1) == v2: return True
    
    if scos is None and use_scos:
        scos = []
    
        i = 0
        for v in bm.verts:
            id = "v" + str(i)
            scos.append(SymVector([sym(id + "_x"), sym(id + "_y"), sym(id + "_z")]))
            i += 1
    
    if use_scos:
        scos2 = [None for v in range(len(bm.verts))]
    
    origfs = list(bm.faces)
    origvs = set(bm.verts)
    origes = list(bm.edges)
    splitvs = set()
    
    cos = []
    
    cos = [v.co.copy() for v in bm.verts]
            
    for v in origvs:
        sum = Vector()
        if use_scos:
            ssum = SymVector()
        tot = 0
        n = 0
        vi = 0
        fi = 0
        val = len(list(v.link_edges))
        
        for f in v.link_faces:
            cent = Vector()
            vi = 0
            
            for v2 in f.verts:
                if v2 == v: continue
                #if v2 not in origvs: continue
                
                w = w1 if has_edge(v, v2) else w6
                
                if use_scos:
                    key = "v" + str(v.index) + "_" + str(fi) + "_" + str(vi)
                    
                    sco = SymVector(scos[v2.index])
                    
                    """
                    if first:
                        sco = SymVector([sym(key + "_x"), sym(key + "_y"), sym(key + "_z")])
                    else:
                        sco = SymVector(scos[v2.index])
                    #"""
                    
                    ssum += sco*w
                    
                sum += cos[v2.index]*w
                tot += w
                
                n += 1
                vi += 1
            fi += 1
            
        #if n*n <= tot: continue
        #n -= 1
        w0 = n*3 #n*n - tot
        w0 = val*val*val*w8
        
        if use_scos:
            ssum += SymVector(scos[v.index])*w0
            
        sum += cos[v.index]*w0
        tot += w0
        
        v.co = sum / tot
        
        if use_scos and tot != 0:
            scos2[v.index] = (ssum / tot).simp()
        elif use_scos:
            scos2[v.index] = sym("div_by_zero_error")
            
    centvs = set()
    if 1:
        for e in origes:
            if use_scos:
                s1 = scos[e.verts[0].index]
                s2 = scos[e.verts[1].index]
                s3 = (s1 + s2)*0.5
                scos.append(s3)
                scos2.append(s3)
            
            co = (cos[e.verts[0].index] + cos[e.verts[1].index])*0.5

            newv = bmesh.utils.edge_split(e, e.verts[0], 0.5)[1] #split_edges(bm, edges=[e])
            newv.co = co
            
            newv.index = len(cos)
            cos.append(co.copy())
            
            splitvs.add(newv)
            
        for f in origfs:
            sum = Vector()
            ssum = SymVector()
            
            tot = 0.0
            
            i = 0
            for l in f.loops:
                v = l.vert
                
                sum += cos[v.index]
                if use_scos:
                    ssum += scos[v.index]
                
                tot += 1
                
            if tot == 0.0: #eek! evil!
                print("evil!!")
                return
             
            if use_scos:
                scentv = ssum / tot
                scos2.append(scentv)
                scos.append(scentv)
                
            centv = bm.verts.new(sum/tot)
            if f.select:
                centv.select = True
            
            centv.index = len(cos)
            
            cos.append(centv.co.copy())
            centvs.add(centv)
            
            i = 0
            for l in f.loops:
                if i % 2 == 0:
                    f2 = bm.faces.new([l.vert, l.link_loop_next.vert, centv, l.link_loop_prev.vert]);
                    
                    if f.select:
                        f2.select = 1
                        
                        for e in f2.edges:
                            e.select = 1
                i += 1
        
                
        for f in origfs:
            bm.faces.remove(f)
        
        for v in splitvs:
            sum = Vector()
            tot = 0.0
            #continue
            for f in v.link_faces:
                #"""
                cv = Vector()
                tot2 = 0
                for v2 in f.verts:
                    cv += cos[v2.index]
                    tot2 += 1
                
                cv /= tot2
                
                sum += cv*w3
                tot += w3
                continue
                #"""
                
                for v2 in f.verts:
                    if v == v2: continue
                    #if v2 in splitvs: continue
                    
                    #w = w3 if has_edge(v, v2) and v2 not in splitvs else w5
                    w = w3# if has_edge(v, v2) else w5
                    #w = w3
                    
                    sum += cos[v2.index]*w
                    tot += w
                    
            w0 = w9
            
            sum += cos[v.index]*w0
            tot += w0
            
            if tot > 0:
                sum /= tot
                v.co = sum
            else:
                pass
                #print("EEEK!")
    
    bm.normal_update()
    if use_scos:
        return scos2
    
def subdivide3(bm, scos=None, use_scos=True, ws=None, first=True,  cos=None, lvl=0, simple=False):
    bm.verts.index_update()
    bm.edges.index_update()
    bm.faces.index_update()
    
    is_first = scos is None
    
    w1 = 6
    w2 = 1
    w3 = 7
    w6 = 1
    w7 = 1
    w8 = 1
    w9 = 1
    w22 = 1
    
    if ws is not None:
       w1, w2, w3, w4, w6, w8, w9 = ws[:7]
    
    
    def has_edge(v1, v2):
        for e in v1.link_edges:
            if e.other_vert(v1) == v2: return e
    
    if scos is None and use_scos:
        scos = []
    
        i = 0
        for v in bm.verts:
            id = "v" + str(i)
            scos.append(SymVector([sym(id + "_x"), sym(id + "_y"), sym(id + "_z")]))
            i += 1
    
    if use_scos:
        scos2 = [None for v in range(len(bm.verts))]
    
    origfs = list(bm.faces)
    origvs = set(bm.verts)
    origes = list(bm.edges)
    splitvs = set()
    
    if cos is None:
        cos = [v.co.copy() for v in bm.verts]
            
    centvs = set()
                
    for v in origvs:
        if simple:
            break
        if not first:
            break
        val = len(list(v.link_edges))
        
        if val < 2: continue
        
        sum = Vector()
        tot = 0

        if use_scos:
            ssum = SymVector()
            
        wa = 3.0/(2*val)
        wb = 1.0/(4*val)
            
        for f2 in v.link_faces:
            cent = Vector()
            vi = 0
            
            vs2 = list(f2.verts)
            for j in range(len(vs2)):
                v2 = vs2[j]
                v3 = vs2[(j+1) % len(vs2)]
                e2 = has_edge(v2, v3)
                if v2 == v: continue
                        
                #if v2 == v or has_edge(v, v2): continue
                
                #v2 = vs2[j]
                
                ww = wa/val*w1 if has_edge(v, v2) else wb/val*w2
                #ww = wb/val
                
                sum += cos[v2.index]*ww
                tot += ww
        
        w0 = (1.0-wa-wb)*w4
            
        sum += cos[v.index]*w0
        tot += w0
        
        dv = (sum / tot) - v.co
        v.co = sum / tot
        
        #no = v.normal.copy()
        #no.normalize()
        
        #v.co += no*dv.length/4

    if 1:
        for e in origes:
            if use_scos:
                s1 = scos2[e.verts[0].index]
                s2 = scos2[e.verts[1].index]
                s3 = (s1 + s2)*0.5
                scos2.append(s3)
            
            co = (cos[e.verts[0].index] + cos[e.verts[1].index])*0.5

            newv = bmesh.utils.edge_split(e, e.verts[0], 0.5)[1] #split_edges(bm, edges=[e])
            newv.co = co
            
            newv.index = len(cos)
            cos.append(co.copy())
            
            splitvs.add(newv)
        
        #cos = [v.co.copy() for v in bm.verts]
        for v in splitvs:
            if simple: break
            
            sum = Vector()
            tot = 0.0
            
            for f in v.link_faces:
                for v2 in f.verts:
                    if v2 == v: continue
                    if has_edge(v, v2): continue
            
                    ww = 1/16
                    sum += cos[v2.index]*ww

                    tot += ww
            
            for e in v.link_edges:
                v2 = e.other_vert(v)
                
                ww = 3/8
                sum += cos[v2.index]*ww
                tot += ww
            
            if tot > 0:
                sum /= tot
                v.co = sum
            else:
                pass
                #print("EEEK!")
            
        for f in origfs:
            sum = Vector()
            ssum = SymVector()
            
            tot = 0.0
            
            i = 0
            for l in f.loops:
                v = l.vert
                e = l.edge
                
                sum += cos[v.index]
                tot += 1.0
                
                if use_scos:
                    ssum += scos2[v.index]
                
                #"""
                l2 = l.link_loop_radial_next
                if l2 != l:
                    f2 = l2.face
                    
                    for l3 in f2.loops:
                        if l3.link_loop_radial_next == l2: continue
                        v2 = l2.vert
                        e2 = l2.edge
                        ww = w3
                        
                        #print("yay", w3)
                        sum += cos[v2.index]*ww
                        tot += ww
                #"""
                
            if tot == 0.0: #eek! evil!
                print("evil!!")
                return
             
            if use_scos:
                scentv = ssum / tot
                scos2.append(scentv)
                
            centv = bm.verts.new(sum/tot)
            if f.select:
                centv.select = True
            
            centv.index = len(cos)
            
            cos.append(centv.co.copy())
            centvs.add(centv)
            
            i = 0
            for l in f.loops:
                if i % 2 == 0:
                    f2 = bm.faces.new([l.vert, l.link_loop_next.vert, centv, l.link_loop_prev.vert]);
                    
                    if f.select:
                        f2.select = 1
                        
                        for e in f2.edges:
                            e.select = 1
                i += 1
        
        
        for v in bm.verts: #origvs:
            #if v not in origvs: continue
            cos[v.index] = v.co.copy()
            
            
        for f in origfs:
            bm.faces.remove(f)
    
    bm.normal_update()
    return cos
    if use_scos:
        return scos2

def subdivide(bm):
    bm.verts.index_update()
    bm.edges.index_update()
    bm.faces.index_update()
    
    def has_edge(v1, v2):
        for e in v1.link_edges:
            if e.other_vert(v1) == v2: return True

    origfs = list(bm.faces)
    origvs = set(bm.verts)
    origes = list(bm.edges)
    splitvs = set()
    
    cos = []
    
    cos = [v.co.copy() for v in bm.verts]
        
    for v in origvs:
        sum = Vector()
        tot = 0
        n = 0
        
        for f in v.link_faces:
            cent = Vector()
            for v2 in f.verts:
                if v2 == v: continue
                #if v2 not in origvs: continue
                
                w = 6 if has_edge(v, v2) else 1
                    
                sum += cos[v2.index]*w
                tot += w
            #    n += 1.0
                
            n += 2.73
        
        if n*n <= tot: continue
        #n -= 1
        w0 = n*n - tot
        
        sum += cos[v.index]*w0
        tot += w0
        
        if tot != 0:
            v.co = sum / tot

    centvs = set()
    if 1:
        for e in origes:
            newv = bmesh.utils.edge_split(e, e.verts[0], 0.5)[1] #split_edges(bm, edges=[e])
            
            newv.index = len(cos)
            cos.append(newv.co.copy())

            splitvs.add(newv)
        for f in origfs:
            sum = Vector()
            tot = 0.0
            
            i = 0
            for l in f.loops:
                v = l.vert
                w = 1 if v in origvs else 0
                
                sum += cos[v.index]*w
                tot += w
                
            if tot == 0.0: #eek! evil!
                print("evil!!")
                return
             
            centv = bm.verts.new(sum/tot)
            if f.select:
                centv.select = True
                
            centv.index = len(cos)
            cos.append(centv.co.copy())
            
            centvs.add(centv)
            
            i = 0
            for l in f.loops:
                if i % 2 == 0:
                    f2 = bm.faces.new([l.vert, l.link_loop_next.vert, centv, l.link_loop_prev.vert]);
                    
                    if f.select:
                        f2.select = 1
                        
                        for e in f2.edges:
                            e.select = 1
                i += 1
        
                
        for f in origfs:
            bm.faces.remove(f)
        
        for v in splitvs:
            sum = Vector()
            tot = 0.0
            
            for f in v.link_faces:
                for v2 in f.verts:
                    if v == v2: continue
                    if v2 in splitvs: continue
                    
                    w = 6 if has_edge(v, v2) and v2 not in splitvs else 1
                    
                    sum += cos[v2.index]*w
                    tot += w
                    
            w0 = 0
            
            sum += cos[v.index]*w0
            tot += w0
            
            if tot > 0:
                sum /= tot
                v.co = sum
            else:
                pass
                #print("EEEK!")
    
    bm.normal_update()
    
class Loop (list):
    pass

class Bez:
    def __init__(self, v1, v2, a, b, c, d):
        self.v1 = v1
        self.v2 = v2
        self.a = Vector(a)
        self.b = Vector(b)
        self.c = Vector(c)
        self.d = Vector(d)
        
        self.length = (self.d - self.a).length
    
    def render(self, bm, steps=32):
        ds = 1.0 / steps
        s = ds
        
        lastv = bm.verts.new(self.eval(0))
        edges = []
        
        for i in range(steps):
            v = bm.verts.new(self.eval(s))
            edges.append(bm.edges.new([lastv, v]))
            
            lastv = v
            s += ds
        
        return edges
    
    def derivative(self, s):
        df = 0.00001

        a = self.eval(s)
        b = self.eval(s + df)
        
        return (b - a) / df
        
    def eval(self, s):
        def cubic(k1, k2, k3, k4, s):
            return -(k1*s**3-3*k1*s**2+3*k1*s-k1-3*k2*s**3+6*k2*s**2-3*k2*s+3*
              k3*s**3-3*k3*s**2-k4*s**3)
        
        a = self.a
        b = self.b
        c = self.c
        d = self.d
        
        x = cubic(a[0], b[0], c[0], d[0], s)
        y = cubic(a[1], b[1], c[1], d[1], s)
        z = cubic(a[2], b[2], c[2], d[2], s)
        
        return Vector([x, y, z])
                 


ws = [
    1, 1, 1, 0.0, 0.1, #0-4
    1, 1, 0.05, 0.06,  #5-8
    0, 4, 3.5,  0, 0, -1, 0,#9-15
    0.3333, 0.3333, 1.0, 1.0, 0.25, -0.25 #16- used in gregory.py
    
    #1, 0.75, -1.0, 0.5, 0.333, 0.333,#16- used in bezpatch.py
    #1.0, 1.0, 0.25, -0.25, 1.0 #22- used in bezpatch.py
]

"""
ws = [1.7448060062265873, 0.8060642425045409, 1.5126530368899969, 0.3100176716003821, 0.3162245563845574, 0, 1, 0.011822678380041361, 0.127833456099555, -0.12206065641666866, 3.911385697038632, 3.513963954701816, 0.17596301927764527, 0.0035423339919746465, -1.1994979196923674, 0.0, 0.4959587290058114, 0.4905257506024106, 1.0019032037538935, 1.0363621323760153, -1.0119696542616317, 0.9540459902135334]
ws = [1, 1, 1.333, -0.05, 0.075, 0, 1, 0.011822678380041361, 0.127833456099555, -0.12206065641666866, 4, 3.5, 0.0, 0.0, -1, 0.0, 0.5, 0.5, 1.0, 1.0, 1.5, -1.0]
ws = [1.0, 1.0, 1.333, -0.05, 0.075, 0, 1, 0.012027207441322822, 0.12044288040777608, -0.0725755885921653, 3.9817362197016006, 3.47407993600624, -0.0010677957856597767, 0.043831980567041864, -1.0171079528695377, 0.0, 0.4937553389188164, 0.49375588017996663, 0.9979098054065706, 0.9972201708047582, 1.4922617127749722, -1.0217693143094295]
ws = [1.0, 1.0, 1.333, -0.05, 0.075, 0, 1, 0.012027207441322822, 0.12044288040777608, -0.0725755885921653, 3.9817362197016006, 3.47407993600624, -0.0010677957856597767, 0.043831980567041864, -1.0171079528695377, 0.0, 0.4937553389188164, 0.49375588017996663, 0.9979098054065706, 0.9972201708047582, 1.4922617127749722, -1.0217693143094295]
ws = [1.0, 1.0, 1.333, -0.05, 0.075, 0, 1, 0.012027207441322822, 0.12044288040777608, -0.0725755885921653, 3.9817362197016006, 3.47407993600624, -0.0010677957856597767, 0.043831980567041864, -1.0171079528695377, 0.0, 0.4937553389188164, 0.49375588017996663, 0.9979098054065706, 0.9972201708047582, 1.4922617127749722, -1.0217693143094295]

ws = [1.0, 1.0, 1.333, -0.05, 0.075, 0, 1, 0.166583107158488, 0.0278466482694525, -0.5836633313762172, 5.138786114616642, 4.9034702586792065, -0.11724608764162522, -0.3846027167958684, -1.6026291354889066, 0.26972222124923095, -0.8636692842022919, 0.4853958585097474, 3.5518154701879614, 1.4082610707348977, 0.6077282739867563, -2.437962243437]
#"""

ws = [1.0, 1.0, 1.333, -0.05, 0.075, 0, 1, 0.17346173411646343, 0.05945459858497371, -0.6958808302214978, 5.144608478995847, 5.850520958668186, -0.17811424543601495, -0.43083925128538625, -2.16779423087856, -0.5382275062474484, -2.2576717797831583, 0.783867962889622, 10.684628518786853, 1.367474639095884, 0.590830270129469, -2.470696670595178]
ws = [1, 1, 1.333, -0.05, -0.02500000000000001, 0, 1, 0.17346173411646343, 0.05945459858497371, -0.6958808302214978, 4, 3.5, -0.05, 0.1, 0, 0, 0.5, 0.5, 1.0, 0.2, 0.0, -0.0]

#ws = [1, 1, 1.333, -0.05, -0.02500000000000001, 0, 1, 0.1510908226696087, 0.0045972439578884005, -0.3664559929580165, 4.280376357194705, 3.939665778751547, -0.10093202212690522, 0.3405227145480962, 0.08348095076124129, 0.0, 0.46012206854075255, 0.4607819908876017, 1.1005332596911035, 0.0650770718198898, -0.005974030170334473, -0.23938821018371537]
ws = [1, 1, 1.333, -0.05, -0.02500000000000001, 0, 1, 0.14792178174257248, 0.00455922624174619, -0.3843950882730577, 4.258829481939746, 3.921287966867826, -0.14428567955625923, 0.31256483936632573, 0.0561006328418586, 0.0, 0.4601739828988156, 0.46083390524566475, 1.102352436718691, 0.12728774384717007, 0.008929575799427426, -0.21597877251232328]
ws = [1, 1, 1.333, -0.05, -0.02500000000000001, 0, 1, 0.1852153723043396, 0.08036582350234255, -0.5871281494464144, 4.048278754367854, 3.6920413249229216, -0.07582706357653493, 0.24559624319800283, -0.3804297303978761, 0.0, 0.5030674225413296, 0.503225626542239, 1.0143310052740904, -0.9819763586124661, -0.17917919576015626, -0.25511214989375486]

ws = [1, 1, 1.333, -0.05, 
      -0.025, 0, 1, 0.211, 0.0321, 
      -0.4710, 4.0808, 3.6905, -0.1325,
      0.4808, -0.6320, 0.0, 0.5443,
      0.6, -0.62, -1.0301, -0.08434, -0.11355
     ]

def find_subsurf_co(v):
    global ws
    
    elen = len(v.link_edges)
    
    if elen == 3:
        w0 = ws[0]
        w1 = ws[1]
        w2 = ws[2]
        w3 = ws[4]
    elif elen == 4:    
        w0 = ws[0]
        w1 = ws[1]
        w2 = ws[2]
        w3 = ws[3]
    else:
        w0 = ws[5]
        w1 = ws[6]
        w2 = ws[7]
        w3 = ws[8]
    
    co = Vector(v.co*w0)
    tot = w0
    
    #if len(list(v.link_edges)) == 2:
        
    for e in v.link_edges:
        v2 = e.other_vert(v)
        co += (v.co+v2.co)*0.5*w1
        tot += w1
    
    for f2 in v.link_faces:
            cent0 = f2.calc_center_median()
            
            co += cent0*w2
            tot += w2
            
            tot2 = 0
            elen = len(list(f2.edges))
            
            for e2 in f2.edges:
                #continue
                for f3 in e2.link_faces:
                    if f3 == f2: continue
                    
                    cent3 = f3.calc_center_median()
                    cent3 = (cent3 + cent0)*0.5
                    
                    co += cent3*w3
                    tot += w3
                
                #co += (e2.verts[0].co + e2.verts[1].co + cent0)*0.33333*w2
                #co += cent0*w2
                #tot += w2
        
    return co / tot

def find_subsurf_boundary(e):
    a = find_subsurf_co(e.verts[0])
    d = find_subsurf_co(e.verts[1])
    b = 0
    c = 0
    
    for i in range(2):
        v = e.verts[i]
        
        w0 = ws[9]
        w1 = ws[10]
        w2 = ws[11]
        w3 = ws[12]
        w4 = ws[13]
        w5 = ws[14] if len(list(e.link_faces)) > 0 else ws[15]
        
        """
        vv = v
        vv = e.verts[i^1]
        if len(list(vv.link_edges)) == 3:
            w0 = 0
            w1 = 0.5
            w2 = 1
            w3 = 0
            w4 = 0
            w5 = 0
            
        elif len(list(vv.link_edges)) > 4:
            w0 = 0
            w1 = 1
            w2 = 1
            w3 = 0
            w4 = 0
            w5 = 0
        #"""
        
        tot = 0.0
        sum = Vector()
        
        sum += e.other_vert(v).co*w5
        tot += w5
        
        sum += v.co*w0
        tot += w0
        
        v2 = e.other_vert(v)
        sum += (v.co + (v2.co - v.co)*0.3333)*w1
        tot += w1
        
        for e2 in v.link_edges:
            v2 = e2.other_vert(v)
            co = v2.co
            
            if e2 == e: continue
            if len(list(e2.link_faces)) == 0.0:
                continue
                
            co = v.co + (co - v.co)*0.5
            
            sum += co*w4;
            tot += w4
            
        for f2 in e.link_faces:
            #continue
            cent0 = f2.calc_center_median()
            
            sum += cent0*w2
            tot += w2
            
            cent = Vector()
            tot2 = 0
            
            for e2 in f2.edges:
                for f3 in e2.link_faces:
                    if f3 == f2: continue
                    
                    cent3 = f3.calc_center_median()
                    
                    sum += cent3*w3
                    tot += w3
                
                #sum += (e2.verts[0].co + e2.verts[1].co + cent0)*0.33333*w3
                #tot += w3
        
        
        if tot == 0:
            sum = Vector(v.co) #isolated vertex
        else:
            sum /= tot
        
        vec = sum - v.co
        #vec = -vec
        #vec *= 0.5
        #sum = vec + v.co
        
        if i == 0:
            b = sum
        else:
            c = sum
    
    #b = Vector(e.verts[0].co)
    #c = Vector(e.verts[1].co)
    
    bz = Bez(e.verts[0], e.verts[1], a, b, c, d)
    
    return bz

def get_edgeloop(bm, selset, v):
    vset = set()
    fset = set()
    
    loop = Loop()
    loop.closed = True
    print("yay")
    
    while v not in vset:
        loop.append(v)
        vset.add(v)
        
        startv = v
        print("yay")
        for e in v.link_edges:
            v2 = e.other_vert(v);
            
            if v2 in selset and v2 not in vset:
                v = v2
                break
        
        if startv == v:
            loop.closed = False
            break
        
    return loop

def get_edgeloops(bm):
    selset = set()
    
    bm.verts.index_update()
    bm.edges.index_update()
    bm.faces.index_update()
    print("yay")
    for v in bm.verts:
        if not v.select or v.hide: continue
        
        selset.add(v)
    
    loops = []
           
    for v in list(selset):
        if v not in selset: continue
        loop = get_edgeloop(bm, selset, v)
        loops.append(loop)
        
        print(loop)
        delset = set()
        
        for v in loop:
            i = 0
            
            for e in v.link_edges:
                v2 = e.other_vert(v)
                i += v2 in selset
            
            #are we not part of another loop?
            #in that case, deselect
            if i < 3:
                delset.add(v)
        for v in delset:
            selset.remove(v)
            
    return loops

def make_bez(loop):
    if len(loop) == 0:
        print("error in make_bez: empty loop, no verts!")
        return Loop()
    
    loop2 = Loop()
    loop2.closed = loop.closed
    
    if 1: #loop.closed:
        loop3 = loop[:]
        #loop3.closed = loop.closed
        loop3.append(loop[0])
        loop = loop3
    
    center = Vector()
    ctot = 0.0
    
    for i in range(0, len(loop)-1):
        i0 = i-1 if i > 0 else i
        i3 = i+2 if i < len(loop)-2 else i
        
        v0 = loop[i0]
        v1 = loop[i]
        v2 = loop[i+1]
        v3 = loop[i3]
        
        a = find_subsurf_co(v1)
        b = find_subsurf_co(v0)
        c = find_subsurf_co(v3)
        d = find_subsurf_co(v2)
        
        n1 = a - b 
        n2 = c - d
        
        n = (d - a)
        #n = (v2.co - v1.co)
        n *= 1
        sz = 0.333
        
        n1 += (n - n1)*0.5
        n2 += (n - n2)*0.5
        
        #n2 *= (v0.co-v3.co).length*0.5;
        n1 *= sz
        n2 *= sz
        
        b = a + n1
        c = d - n2
        
        bez = Bez(v1, v2, a, b, c, d)
        loop2.append(bez)
        
        center += v1.co
        ctot += 1
        
    
    if ctot > 0:
        center /= ctot
    loop2.center = center
    
    normal = Vector();
    
    for i in range(0, len(loop)-1):
        v1 = loop[i]
        v2 = loop[i+1]
        
        t1 = v1.co - center
        t2 = v2.co - v1.co
        n2 = t1.cross(t2)
        n2.normalize()
        
        normal += n2
    
    normal.normalize()
    print(normal)
    
    loop2.normal = normal
    
    return loop2    
        

def loopto2d(loop2):
    n = loop2.normal
    
    ax = abs(n[0])
    bx = abs(n[1])
    cx = abs(n[2])
    
    #find axis for up vector
    if ax < cx and bx < cx:
        axis = 0
    elif bx < ax and cx < ax:
        axis = 1
    else:
        axis = 2
    
    world = Vector()
    world[axis] = 1.0
    
    z = n
    y = z.cross(world)
    x = y.cross(z)
    
    mat = Matrix([x, y, z])
    mat.to_4x4()
    
    #mat.transpose()
    for bz in loop2:
        bz.a = mat * bz.a
        bz.b = mat * bz.b
        bz.c = mat * bz.c
        bz.d = mat * bz.d
        
    return mat
        
def faircurve(loop, origs):
    def minmax(l):
        vmin = Vector([1e17, 1e17, 1e17])
        vmax = Vector(-vmin)
        
        for v in l:
            for i in range(3):
                vmin[i] = min(vmin[i], v.co[i])
                vmax[i] = max(vmax[i], v.co[i])
                
        return (vmin, vmax)
    
    startbb = minmax(loop)

    if len(loop) == 0:
        print("warning, loop had no verts");
        return
    
    vset = set()
    startcos = {}
    
    for i, v in enumerate(loop):
        v.index = i    
        vset.add(v)
        startcos[v] = Vector(v.co)
        
    def open_prev(v):
        return loop[v.index-1] if v.index > 0 else v
    def open_next(v):
        return loop[v.index+1] if v.index < loop.length-1 else v
    def close_prev(v):
        return loop[(v.index-1+len(loop)) % len(loop)]
    def close_next(v):
        return loop[(v.index+1) % len(loop)]
    
    closed = False
    for e in loop[0].link_edges:
        if e.other_vert(loop[0]) == loop[len(loop)-1]:
            closed = True
    
    if closed:
        prev = close_prev
        next = close_next
    else:
        prev = open_prev
        next = open_next
    
    #finite difference width
    #is twice of tesselated edge length
    df = (loop[2].co - loop[1].co).length*1.01
    
    def walkback(v, dis):
        walklen = 0
        while walklen < dis and v != prev(v):
            tmp = v
            v = prev(v)
            walklen += (tmp.co - v.co).length
            
        return v
    
    def walknext(v, dis):
        walklen = 0
        while walklen < dis and v != next(v):
            tmp = v
            v = next(v)
            walklen += (tmp.co - v.co).length
            
        return v
    
    def derivative(v):
        v0 = walkback(v, df)
        v2 = walknext(v, df)
        
        return (v2.co - v0.co) / (2 * df)
    
    def derivative2(v):
        v1 = walkback(v, df).co
        v2 = walkback(v, df*2).co
        
        v3 = walknext(v, df).co
        v4 = walknext(v, df*2).co
                    
        dva = (v4 - v1) / (2 * df)
        dvb = (v3 - v2) / (2 * df)
        
        dv1 = dva
        dv2 = (dvb - dva) / df
        
        return dv2
    
    def curvature(v):
        v1 = walkback(v, df).co
        v2 = walknext(v, df).co
        
        vv = walknext(v, df)
        
        v3 = walkback(vv, df).co
        v4 = walknext(vv, df).co
                    
        dva = (v2 - v1) / (2 * df)
        dvb = (v3 - v2) / (2 * df)
        
        dv1 = dva
        dv2 = (dvb - dva) / df
        
        k = dv1[0]*dv2[1] - dv1[1]*dv2[0]
        d = dv1[0]*dv1[0] + dv1[1]*dv1[1]
        
        if d != 0.0:
            k /= pow(d, 3/2)
            
        return k
    
    def dcurvature(v):
        #v1 = walkback(v, df*2)
        #v2 = walknext(v, df*2)
        
        df2 = 0.5
            
        k1 = curvature(walkback(v, df2))
        k2 = curvature(v)
        k3 = curvature(walknext(v, df2))
        
        dv = (k2 - k1) / df2
        dv2 = (k3 - k2) / df2
        
        dv2 = (dv2 - dv) / df2
        
        return dv2
    
    elens = {}
    eset = set()
    
    for v in loop:
        v.normal = Vector()
        
        for e in v.link_edges:
            if e.other_vert(v) in vset:
                eset.add(e)
                elens[e] = (e.verts[1].co - e.verts[0].co).length
    
    for step in range(55):
        err = 0.0
        
        cos = [v.co.copy() for v in loop]
        #for i, v in enumerate(loop):
        for _i in range(len(loop)):
            #i = int(random.random()*len(loop)*0.99999)
            i = _i
            v = loop[i]
            
            gs = Vector()
            k1 = dcurvature(v)
            
            err += abs(k1)
                
            for j in range(2):
                df2 = 0.001
                orig = v.co[j]
                
                v.co[j] += df2
                k2 = dcurvature(v)
                gs[j] = (k2 - k1) / df2
                
                v.co[j] = orig
             
            totg = gs[0]*gs[0] + gs[1]*gs[1]
            if totg == 0.0: continue
        
            #k1 = 1.0 / sqrt(totg);
            k1 = k1 / totg
            
            for j in range(2):
                if gs[j] == 0: continue
                
                v.normal[j] += -k1*gs[j]*0.5
                v.co[j] += v.normal[j]
                v.normal[j] *= 0.0
                
            
            #go towards original positions a bit
            v.co += (startcos[v] - v.co)*0.01
            
            #if gs[0] != 0.0: print(k1/gs[0])
            #cos[i] += derivative2(v)*2
        
        for v in loop:
            #apply a spring constraint
            for e in v.link_edges:
                if e not in eset: continue
                #continue#
            
                elen = elens[e] #rest length
                elen2 = (e.verts[1].co-e.verts[0].co).length
                
                v1, v2 = e.verts
                
                elen = elen/elen2*0.5
                cent = (v1.co + v2.co)*0.5
                
                co2 = (v1.co-cent) * elen + cent
                co3 = (v2.co-cent) * elen + cent
                
                v1.co += (co2 - v1.co)*0.75
                v2.co += (co3 - v2.co)*0.75
                
        print("err:", err/len(loop))
    
    endbb = minmax(loop)
    v1 = (startbb[1] -  startbb[0])
    v2 = (endbb[1] - endbb[0])
    cent = (startbb[0] + startbb[1]) * 0.5
    cent2 = (endbb[0] + endbb[1]) * 0.5
    
    #v1[2] = v2[2] = 0.0
    ratio = Vector([v1[0] / v2[0], v1[1] / v2[1], 1.0])
    #ratio = (ratio[0]+ratio[1])*0.5
    ratio = max(ratio[0], ratio[1])
    
    for v in loop:
        #break
        for i in range(2):
            v.co[i] = (v.co[i] - cent2[i])*ratio + cent[i]
        v.co[2] += (cent2[2] - v.co[2])*0.25

"""        
print("YYYY")   
#print("\n\n")          
scos = [find_sco(v) for v in bm.verts]

for i, v in enumerate(bm.verts):
    pass 
    #v.co = scos[i]
    
loops = get_edgeloops(bm)
    
for l in loops:
    l2 = make_bez(l)
    endv = None
    firstv = None
    
    mat = loopto2d(l2)
    mat2 = mat.copy()
    mat.invert()
    
    l3 = Loop()
    l3.closed = l2.closed
    
    origs = []
    
    bm.verts.index_update()
        
    for i, bz in enumerate(l2):
        nextbz = l2[(i+1)%len(l2)]
        
        startv = outbm.verts.new(bz.eval(0.0))
        origs.append([bz.v1, startv, Vector(startv.co), bz.v1.index])
        
        lastv = startv
        
        l3.append(startv)
        
        if firstv is None:
            firstv = startv
        
        if endv is not None:
            outbm.edges.new([endv, startv]);
            pass
            
        steps = 8
        ds = 1.0 / float(steps+1)
        s = ds
        
        ds = 0.1;
        
        def nexts(co, s):
            s2 = s + ds/bz.length
            s2 = min(max(s2, 0.0), 1.0)
            
            start = s
            end = s + ds/bz.length*2.0
            
            for i in range(3):
                co2 = bz.eval(s2) 
                dis = (co2 - co).length
                
                if dis > ds:
                    end = s2
                elif dis < ds:
                    start = s2
                else:
                    break
                
                s2 = (start+end)*0.5

            s2 = min(max(s2, 0.0), 1.0)
            
            return s2
            
        while s < 1.0:
            co = bz.eval(s)
            v = outbm.verts.new(co)
            l3.append(v)
            
            if lastv is not None:
                e = outbm.edges.new([lastv, v])
                
            lastv = v
            s = nexts(co, s)
            
            if s >= 1.0:
                break
        
        endv = lastv
    
    if firstv is not None:
        outbm.edges.new([lastv, firstv])
    
    outbm.verts.index_update()
    outbm.edges.index_update()
    
    faircurve(l3, origs)

    for v, v2, v3, idx in origs:
        co2 = v2.co
        v.co = mat2 * v.co
        v4 = -(v3-v2.co)*0.5
        v.co += (v2.co - v.co)*0.5 #v4
        #v.co = v3

        v.co = mat * v.co
    
    for bz in l3:
        bz.co = mat * bz.co
          
    break

for e in list(outbm.edges):
    #outbm.edges.remove(e)
    pass
    
for v in bm.verts:
    sco = find_sco(v)
    #v2 = outbm.verts.new(sco)
"""