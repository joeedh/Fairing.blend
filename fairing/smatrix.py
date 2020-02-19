import bpy, bmesh, random, time
from math import *
from mathutils import *

import bpy, bmesh, time, random
from mathutils import *
from math import *
from .sym import sym

import numpy as np

def mulVertsMatrix(vs, smat):
    for i in range(3):
        cos = [v[i] for v in vs]
        cos = cos @ smat
        for j in range(len(vs)):
            vs[j][i] = cos[j]
     
     
def getVertDisk(v):
    vs = [v]
    for l in v.link_loops:
        if l.vert == v:
            vs.append(l.link_loop_next.vert)
            vs.append(l.link_loop_next.link_loop_next.vert)
        else:
            vs.append(l.vert)
            vs.append(l.link_loop_prev.vert)
    return vs

def getSplitEdgeVertDisk(newv, origv1, origv2):
    vs = [newv, origv1, origv2]
    vs2 = getVertDisk(newv)
    
    for v in vs2:
        if v not in vs:
            vs.append(v)
    
    return vs
    
def getSplitSMat(): #for [new] split vertices
    smat = np.identity(9)
    
    wa = 1/16
    wb = 3/8
    #wb=0
    
    smat[0][0] = 0
    smat[1][0] = 3/8
    smat[2][0] = 3/8
    smat[3][0] = wa
    smat[4][0] = wb
    smat[5][0] = wa
    smat[6][0] = wb
    smat[7][0] = wa
    smat[8][0] = wb
    
    #return smat, None
    
    n = len(smat)
    
    for i in range(n):
        wsum = 0
        #break
        for j in range(n):
            wsum += smat[j, i]
            
        for j in range(n):
            smat[j, i] /= wsum
    
    #ismat = np.linalg.inv(smat)
    return smat, None


def getBoundFlags(v):
    return len(v.link_edges) - len(v.link_loops)
        
#for original cage verts   
cache = {} 
def getSMat(val, boundary_flag=0):
    key = val + (boundary_flag<<10)
    
    global cache
    if key in cache:
        return cache[key]
    
    n = 2*val + 1
    smat = np.identity(n)
        
    wa = 3.0/(2*val)
    wb = 1.0/(4*val)
    w0 = (1.0-wa-wb)
    
    if boundary_flag:
        #w0 = 3/4
        wa = 0
        wb = 0
        pass
        
    wsum = 0.0
    
    for i in range(1, n):
        w = wa/val if i%2==0 else wb/val
        
        smat[i, 0] = w
        wsum += w
    
    smat[0][0] = w0
    
    for i in range(1, n):
        iprev = 1 + ((i - 2 + 2*val) % (2*val))
        inext = 1 + (i % (2*val))
        
        """
        smat[i][i] = w0
        smat[iprev][i] = wa
        smat[inext][i] = wa
        smat[0][i] = wb
        #"""
        
        """
        smat[i][i] = w0
        smat[i][iprev] = wa
        smat[i][inext] = wa
        smat[i][0] = wb
        #"""
    
    n2 = n - boundary_flag
    for i in range(n):
        wsum = 0
        for j in range(n2):
            wsum += smat[j, i]
            
        for j in range(n2):
            smat[j, i] /= wsum
            pass
    
    ismat = np.linalg.inv(smat)
    
    cache[key] = (smat, ismat)
    return cache[key]
    
# ############  IGNORE stuff below this line  #############

class SMatrix (list):
    def __init__(self, n):
        self[:] = [[0 for x in range(n)] for y in range(n)]

class Disk:
    def __init__(self, vs, es, fs):
        self.fs = fs
        self.es = es
        self.vs = vs
        
def getdisk(v):
    l = None
    for l2 in v.link_loops:
        if l2.vert == v:
            l = l2
            break

    if l is None:
        print("failed")
        return None
    
    vs = [l.vert]
    es = [l.edge]
    fs = [l.face]
    
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
    
    return ret
    
#assumes all-quad perfectly manifold mesh with no holes
def buildSMatrix(v):
    val = len(v.link_edges)
    smat = SMatrix(val)
    
    disk = getdisk(v)
    
    print(len(disk.vs), len(disk.es), len(disk.fs), val)
    
def test(bm,  visbm):
    subdivide(bm)
    
    for v in bm.verts:
        if not v.select: continue
        buildSMatrix(v)
    
        