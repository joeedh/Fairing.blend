import bpy, bmesh
from mathutils import *
from math import *
import time, random, bmesh

def ob_force_update(ob):
    graph = bpy.context.evaluated_depsgraph_get()
    graph.update()
    obeval = ob.evaluated_get(graph)
    
    return obeval
    
barrier = 0

def ob_get_final_mesh(ob):
    global barrier
    
    if barrier:
        print("depsgraph called tool; loop detected")
        return
    
    barrier = 1
    obeval = ob_force_update(ob)

    #XXX this is crashing. . .
    ret = obeval.to_mesh()
    
    #ret = ob.to_mesh()
    
    barrier = 0
    return ret
    
def getShadow(ob, link_to_scene=True, wireframe=True, ctx=None):
    ctx = bpy.context if ctx is None else ctx
    
    key = "_" + ob.name + "_fshd"
    
    if key not in bpy.data.objects:
        me = bpy.data.meshes.new(key)
        ob2 = bpy.data.objects.new(key, me)
    else:
        ob2 = bpy.data.objects[key]
    
    if ob.mode == "EDIT":
        bm = bmesh.from_edit_mesh(ob.data)
    else:
        bm = bmesh.new()
        bm.from_mesh(ob.data)
        
    bm.to_mesh(ob2.data)
    ob2.data.update()
    
    ob2.location = ob.location
    ob2.scale = ob.scale
    ob2.rotation_euler = ob.rotation_euler
    
    if wireframe:
        ob2.display_type = "WIRE"
    else:
        ob2.display_type = "TEXTURED"
        
    if link_to_scene:
        view_layer = ctx.view_layer
        
        if ob2.name not in view_layer.active_layer_collection.collection.objects:
            view_layer.active_layer_collection.collection.objects.link(ob2)
            
    print(ob2)
    
    return ob2

def getSubSurf(ob):
    for mod in ob.modifiers:
        if mod.type == "SUBSURF" and mod.show_viewport:
            return mod
    
    return None
    
def getUnSubShadow(ob, link_to_scene=True, wireframe=True, ctx=None):
    ctx = bpy.context if ctx is None else ctx
    ob2 = getShadow(ob, link_to_scene, wireframe)
    
    me = ob2.data

    bm = bmesh.new()
    bm.from_mesh(me)
    
    bm.verts.index_update()
    bm.edges.index_update()
    
    cos = [Vector(v.co) for v in bm.verts]
    nos = [Vector(v.normal) for v in bm.verts]
    
    n1 = Vector()
    n2 = Vector()
    no = Vector()
    
    """
    for v in bm.verts:
        lens = 0
        tot = 0
        for e in v.link_edges:
            lens += (cos[e.verts[1].index] - cos[e.verts[0].index]).length
            tot += 1
        
        th = 0
        totth = 0
        
        for e1 in v.link_edges:
            for e2 in v.link_edges:
                if e1 == e2: continue
                
                n1[:] = cos[e1.verts[0].index] - cos[e1.verts[1].index]
                n2[:] = cos[e2.verts[0].index] - cos[e2.verts[1].index]
                
                n1.normalize()
                n2.normalize()
                
                d = (n1.dot(n2))
                
                th += d
                totth += 1
                
        if tot == 0 or totth == 0: continue
        #lens /= tot
        #th /= totth
        
        #th = 1.0 - th
        
        no[:] = nos[v.index]
        no.normalize()
        
        print(th)
        
        if tot == 3:
            no *= 2.0
        
        v.co += no * lens * 0.05
    #"""
    
    if "origindex" not in bm.faces.layers.int:
        bm.faces.layers.int.new("origindex")
    lay = bm.faces.layers.int.get("origindex")
    
    for i, f in enumerate(bm.faces):
        f[lay] = i
    
    vs = []
    for v in bm.verts:
        if v.select and not v.hide:
            vs.append(v)
    
    for v in vs:
        for f in v.link_faces:
            for v2 in f.verts:
                for f2 in v2.link_faces:
                    f2[lay] = -1
            
    for v in vs:
        break
        try:
            bmesh.utils.vert_dissolve(v)
        except ReferenceError:
            print("v removed", v)
            pass
            
    me = ob2.data
    bm.to_mesh(me)
    
    if "origindex" not in me.polygon_layers_int:
        me.polygon_layers_int.new(name="origindex")
    
    #layer = me.polygon_layers_int["origindex"].data
    #for i, p in enumerate(me.polygons):
    #    layer[i].value = i
        
    me.update()
    ob2.modifiers.clear()
    
    #mirror
    for mod in ob.modifiers:
        if mod.type == "MIRROR":
            mmod = ob2.modifiers.new("mirror", "MIRROR")
            mmod.mirror_object = mod.mirror_object
            
            for k in dir(mod):
                if k.startswith("_"): continue
                
                v = getattr(mod, k)
                if not hasattr(v, "__call__"): #type(v) in [bool, int, float, str]:
                    try:
                        setattr(mmod, k, v)
                    except AttributeError:
                        pass #read-only?
                
    #decimate?
    dmod = ob2.modifiers.new("decimate", "DECIMATE")
    dmod.decimate_type = "UNSUBDIV"
    dmod.iterations = 1

    smod = ob2.modifiers.new("subsurf", "SUBSURF")
    
    ss = getSubSurf(ob)
    if ss is not None:
        smod.levels = ss.levels + 1
    else:
        smod.levels = 3
        
    return ob2







