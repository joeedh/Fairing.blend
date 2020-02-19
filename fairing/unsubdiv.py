import bpy, bmesh, random, time, struct
from math import *
from mathutils import *

import bpy, bmesh, time, random, os, os.path
from ctypes import *
from mathutils import *
from math import *
from .sym import sym
from .bezpatch import BezPatch
from .subsurf import subdivide2, subdivide, subdivide3, catmull_clark
from . import globals

from .smatrix import getSMat, getVertDisk, getBoundFlags
from .subsurf_evaluate import *
import numpy as np
import numpy.linalg
#from .smatrix import *

def edgeloop_next(l):
    l = l.link_loop_next.link_loop_radial_next.link_loop_next
    
    return l
    
def edgeloop_prev(l):
    l = l.link_loop_prev.link_loop_radial_prev.link_loop_prev
    
    #l = l.link_loop_radial_prev
    return l

def unsubdivide_vert(v, mul_smat=False):
    co = Vector(v.co)
    val = len(v.link_edges)
    
    vs = [v]
    for l in v.link_loops:
        if l.vert == v:
            vs.append(l.link_loop_next.vert)
            vs.append(l.link_loop_next.link_loop_next.vert)
        else:
            vs.append(l.vert)
            vs.append(l.link_loop_prev.vert)
    
    n = 2*val + 1
    is_boundary = getBoundFlags(v)
    
    while len(vs) < n:
        vs.append(v) #XXX do boundary case rules
        
    vcos = [Vector() for x in range(n)]
    cos = [0, 0, 0]
    smat, ismat = getSMat(val, is_boundary)
    
    """
    if mul_smat:
        smat2 = smat @ smat
        try:
            ismat = np.linalg.inv(smat2)
        except:
            pass
        pass
    #"""
    
    #""" row vector
    for i in range(3):
        #print(val, len(cos[i]), len(ismat))
        cos[i] = [v.co[i] for v in vs]
        cos[i] = cos[i] @ ismat
        
        for j in range(len(cos)):
            vcos[j][i] = cos[i][j]
    #"""
    
    return vcos[0]
    print(cos[0])
    #these three add to one
    wa = 3.0/(2*val)
    wb = 1.0/(4*val)
    w0 = (1.0-wa-wb)
    
    co = Vector()
    for i in range(1, len(vs)):
        w = wa if i % 2 == 0 else wb
        w /= val
        
        v2 = vs[i]
        co += -v2.co*w
    
    co += v.co
    co /= w0
    
    return co
 
def unsubdivide(bm, origvs, kill_old=True, mul_smat=False):
    bm.verts.index_update()
    cos = [Vector(v.co) for v in bm.verts]
    
    for v in origvs:
        cos[v.index] = unsubdivide_vert(v, mul_smat)
    
    for v in origvs:
        v.co = cos[v.index]
    
    if not kill_old:
        return
        
    quads = []
    for v in origvs:
        #if v.index != 55: continue
        
        for l in v.link_loops:
            if l.vert != v:
                continue
                l = l.link_loop_prev.link_loop_prev
            
            l1 = l
            l2 = edgeloop_next(l).link_loop_next
            l3 = edgeloop_next(l2).link_loop_next
            l4 = edgeloop_next(l3).link_loop_next
            
            v1 = l1.vert
            v2 = l2.vert
            v3 = l3.vert
            v4 = l4.vert
            
            quads.append([v1, v2, v3, v4])
    
    for f in list(bm.faces):
        bm.faces.remove(f)
        
    for q in quads:
        v1, v2, v3, v4 = q
        try:
            bm.faces.new([v1, v2, v3, v4])
        except ValueError: #face exists
            pass
        
    for e in list(bm.edges):
        ok = True
        for f in e.link_faces:
            ok = False
        if ok:
            bm.edges.remove(e)
            
    for v in list(bm.verts):
        ok = True
        for e in v.link_edges:
            ok = False
        if ok:
            bm.verts.remove(v)
    
def test(inob_name, outob_name, levels=1):
    ob = bpy.data.objects[inob_name]
    bm = bmesh.new()
    
    bm.from_mesh(ob.data)
    outbm = bm.copy()
    
    bm.verts.index_update()
    bm.edges.index_update()
    bm.faces.index_update()
    outbm.verts.ensure_lookup_table()
    
    origcos = [Vector(v.co) for v in bm.verts]
    
    vset = set(bm.verts)
    origvs = [vset]
    bm2 = bm.copy()
    bm2.verts.index_update()
    origvs2 = set(bm2.verts)
    
    for i in range(levels):
        catmull_clark(bm)
        
        if i == 0:
            bm2 = bm.copy()
            bm2.verts.index_update()
            bm2.verts.ensure_lookup_table()
            #origvs2 = set(bm2.verts)
            
            origvs2 = set()
            for v in vset:
                v2 = bm2.verts[v.index]
                origvs2.add(v2)
        
        origvs.append(set(bm.verts))
    
    ocos2 = [Vector(v.co) for v in bm.verts]
    for i in range(2):
        for v in bm.verts:
            v.co += v.normal*0.2
            pass
        
        for j in range(10):
            bmesh.ops.smooth_vert(bm, verts=bm.verts, factor=1.0, use_axis_x=True, use_axis_y=True, use_axis_z=True)
            pass
    
    """
    outob = bpy.data.objects[outob_name]
    bm.to_mesh(outob.data)
    outob.data.update()
    return
    #"""
    
    """
    bm2 = outbm.copy()
    bm2.verts.index_update()
    origvs2 = set(bm2.verts)
    catmull_clark(bm2)
    #"""
    
    bm.verts.index_update()
    bm2.verts.index_update()
    bm2.verts.ensure_lookup_table()
    bm.verts.ensure_lookup_table()
    
    #print(len(origvs[0]))
    
    for v in origvs[1]:
        #v = bm2.verts[i]
        bm2.verts[v.index].co = bm.verts[v.index].co
        pass
        
    unsubdivide(bm2, origvs2, False, True)
    
    """
    outob = bpy.data.objects[outob_name]
    bm2.to_mesh(outob.data)
    outob.data.update()
    return
    #"""
    
    #for i in range(levels):
    #    unsubdivide(bm, origvs[levels-1-i], i != levels-1) #i!=0)
    
    for v in vset:
        v = bm2.verts[v.index]
        dv = v.co - origcos[v.index]
        
        #if not v.select: continue
            
        outbm.verts[v.index].co += dv
        #outbm.verts[v.index].co = v.co
        
    outob = bpy.data.objects[outob_name]
    outbm.to_mesh(outob.data)
    outob.data.update()
    
    
    