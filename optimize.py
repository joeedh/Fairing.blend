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

class OptState:
    def __init__(self, bm):
        self.bm = bm
        self.sm = None
        
    def init(self):
        bm = self.bm
        sm = self.sm = SSMesh()
        sm.load(bm)
        

def test(inob_name, outob_name, steps=3):
    ob = bpy.data.objects[inob_name]
    bm = bmesh.new()
    
    bm.from_mesh(ob.data)
    outbm = bmesh.new()
    
    optstate = OptState(bm)
    optstate.init()
    
    optstate.sm.renderVal4Cages(outbm, False)
    optstate.sm.renderDerivatives(outbm)
    
    outob = bpy.data.objects[outob_name]
    outbm.to_mesh(outob.data)
    outob.data.update()
    
    
    
    
    