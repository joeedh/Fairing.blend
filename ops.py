from . import util, shadow
from . import simplemesh

import bpy, bmesh, time, random
from bpy_extras import view3d_utils
from mathutils import *
from math import *

from bpy.props import *

import bgl
import blf

class AXES:
    X = 1
    Y = 2
    Z = 4

colormap = [
    Vector([0, 0, 1, 1]),
    Vector([1, 0, 1, 1]),
    Vector([0, 1, 0, 1]),
    Vector([0, 1, 1, 1]),
    Vector([1, 1, 0, 1]),
    Vector([1, 0, 0, 1])
]

def fract(f):
    return f - floor(f)
def colormap_get(f):
    f = min(max(f, 0.0), 1.0)
    t = fract(f*len(colormap)*0.9999999)
    f = int(f*len(colormap)*0.9999999)
    
    if f >= len(colormap)-1:
        return colormap[f]
    else:
        a = colormap[f]
        b = colormap[f+1]
        return a + (b - a) * t
    
def handle_mirror_clipping(self, ob, bm, vcos):
    axes = 0
    limit = 0.0
    
    for mod in ob.modifiers:
        #print(mod.type)
        if mod.type != "MIRROR" or not mod.use_clip: continue
        
        if mod.use_axis[0]:
            axes |= AXES.X
        if mod.use_axis[1]:
            axes |= AXES.Y
        if mod.use_axis[2]:
            axes |= AXES.Z
        
        limit = max(limit, mod.merge_threshold)
    
    for i, v in enumerate(bm.verts):
        if not v.select or v.hide: continue
        
        for j in range(3):
            if not (axes & (1<<j)):
                continue
            
            d = abs(vcos[i][j])
            if d <= limit:
                v.co[j] = 0
                
def draw_callback_3dpx(self, context):
    if not hasattr(self, "sm"):
        print("no 3d draw data")
        return
        
    matrix = bpy.context.region_data.perspective_matrix
    sm = self.sm
    
    bgl.glDisable(bgl.GL_DEPTH_TEST)
    bgl.glEnable(bgl.GL_DEPTH_TEST)
    
    #bgl.glPolygonOffset(100000, 100000);
    
    #bgl.glDisable(bgl.GL_BLEND)
    bgl.glEnable(bgl.GL_BLEND)
    sm.draw({
        "uColor" : [0.7, 0.8, 1, 0.3],
        "viewProjectionMatrix" : matrix,
        "uPolyOffset" : 0.5
    })
    
    #bgl.glEnable(bgl.GL_BLEND)
    if self.sm2 is not None:
        sm2 = self.sm2
        sm2.draw({
            "uColor" : [1, 1, 1, 0.7],
            "viewProjectionMatrix" : matrix,
            "uPolyOffset" : 0.5
        })
        
    if self.sm3 is not None:
        self.sm3.draw({
            "uColor" : [1, 1, 1, 0.5],
            "viewProjectionMatrix" : matrix,
            "uPolyOffset" : 0.5
        })
    
    bgl.glDisable(bgl.GL_DEPTH_TEST)
    
    bgl.glEnable(bgl.GL_DEPTH_TEST)
    
def draw_callback_px(self, context):
    #print("mouse points", len(self.mouse_path))

    font_id = 0  # XXX, need to find out how best to get this.

    area = context.area
    w = area.width
    h = area.height
    
    self.text = "yay"
    # draw some text
    if self.text != "":
        blf.position(font_id, 15, 30, 0)
        blf.size(font_id, 20, 72)
        blf.draw(font_id, self.text)
    
    sm = simplemesh.SimpleMesh()
    d = 100
    #sm.tri([0,0,0], [0, d, 0], [d, d, 0])
    
    for l in self._lines:
        v1 = [(l[0][0]-w*0.5)/w*2.0, (l[0][1]-h*0.5)/h*2.0, 0]
        v2 = [(l[1][0]-w*0.5)/w*2.0, (l[1][1]-h*0.5)/h*2.0, 0]
        v1[0] = v1[1] = 0
        #print(v1, v2)
        sm.line(v1, v2)
    
    #sm.line([0, 0, 0], [d, d, 0])
    
    sm.draw({
        "uColor" : [1, 1, 1, 0.75]
    })
    
    #print(dir(bgl))
    return
    # 50% alpha, 2 pixel width line
    bgl.glEnable(bgl.GL_BLEND)
    bgl.glColor4f(0.0, 0.0, 0.0, 0.5)
    bgl.glLineWidth(2)

    bgl.glBegin(bgl.GL_LINE_STRIP)
    for x, y in self.mouse_path:
        bgl.glVertex2i(x, y)

    bgl.glEnd()

    # restore opengl defaults
    bgl.glLineWidth(1)
    bgl.glDisable(bgl.GL_BLEND)
    bgl.glColor4f(0.0, 0.0, 0.0, 1.0)

class ProjectToUnSubD(bpy.types.Operator):
    """Modal object selection with a ray cast"""
    bl_idname = "mesh.project_to_unsubd"
    bl_label = "UnSubD Project"
    bl_options = {'REGISTER', 'UNDO'}
    
    factor: bpy.props.FloatProperty(name="factor")

    sm2 = None
    sm3 = None
    
    def main(self, context):
        # get the context arguments
        scene = context.scene
        region = context.region
        rv3d = context.region_data
        
        LT = simplemesh.LayerTypes
        lf = LT.VERTEX | LT.NORMALS | LT.UVS | LT.COLORS
        self.sm3 = simplemesh.SimpleMesh(shader=simplemesh.SimpleShader3D, layerflag=lf)
        sm3 = self.sm3

        dl = self.factor
        
        def obj_ray_cast(obj, matrix_inv, ray_origin, ray_target):
            """Wrapper for ray casting that moves the ray into object space"""

            #ray_origin_obj = matrix_inv @ ray_origin
            #ray_target_obj = matrix_inv @ ray_target
            ray_direction = ray_target - ray_origin

            # cast the ray
            success, location, normal, face_index = obj.ray_cast(ray_origin, ray_direction)
            
            dist = (ray_origin - location).length
            
            if success:
                return location, dist, normal, face_index
            else:
                return None
        
        ob = context.active_object
        if ob is None or ob.type != "MESH" or ob.name.startswith("_") or ob.mode != "EDIT":
            print("invalid object", ob)
            return
        
        
        ob2 = self._ob2
       
        self.sm2 = simplemesh.SimpleMesh(shader=simplemesh.SimpleShader3D)
        sm2 = self.sm2
        
        bm = bmesh.from_edit_mesh(ob.data)
        
        cos = self.cos
        nos = self.nos
        
        for i, v in enumerate(bm.verts):
            v.co = cos[i]
         
        # get the ray relative to the object
        matrix_inv = ob.matrix_world.inverted()
        
        dav = 0
        dtot = 0
        matrix = ob.matrix_world
            
        for i, v in enumerate(bm.verts):
            if not v.select or v.hide: continue
            
            no = v.normal
            target = v.co + no*1000
            
            ret = obj_ray_cast(ob2, matrix_inv, v.co, target)

            target = v.co - no*1000
            ret2 = obj_ray_cast(ob2, matrix_inv, v.co, target)
            
            if ret is None and ret2 is None: continue
            elif ret is not None and ret2 is not None:
                if ret2[1] < ret[1]:
                    ret = ret2
            elif ret is None and ret2 is not None:
                ret = ret2
            
            no = Vector(v.normal)
            no.normalize()
            
            v.co = cos[i] + (ret[0] - cos[i]) * dl
            dist = (v.co - cos[i]).length
            dav += dist
            dtot += 1
            
            sm2.line(matrix @ v.co, matrix @ Vector(ret[0]))
            sm2.line(matrix @ v.co, matrix @ Vector(ret[0]))
        
        for e in bm.edges:
            ok = not e.verts[0].hide and e.verts[0].select
            ok = ok or (not e.verts[1].hide and e.verts[1].select)
            
            if not ok: continue
            
            v1, v2 = e.verts
            #sm3.line(matrix @ v1.co, matrix @ v2.co)
            
        if dtot > 1:
            dav /= dtot
            
            for i, v in enumerate(bm.verts):
                if not v.select or v.hide: continue
                
                no = Vector(nos[i])
                no.normalize()                
                
                sl = -1 if dl < 0 else 1
                v.co += no*sl*dav
            
                
        handle_mirror_clipping(self, ob, bm, self.cos)
        bmesh.update_edit_mesh(ob.data, destructive=False)

    @classmethod
    def poll(cls, context):
        return (context.mode == 'EDIT_MESH')

    def execute(self, context):
        self._ob2 = shadow.getUnSubShadow(context.active_object, ctx=context)
        self._lines = []
        
        self.main(context)
        self.stop()
        return {'FINISHED'}
        
    def makeDrawData(self, ob2):
        me = shadow.ob_get_final_mesh(ob2)
        
        LT = simplemesh.LayerTypes
        lf = LT.VERTEX | LT.NORMALS  | LT.UVS | LT.COLORS
        
        drawbm = bmesh.new()
        self.sm = simplemesh.SimpleMesh(shader=simplemesh.SimpleShader3D, layerflag=lf)

        fset = set()
        sm = self.sm
        layer = me.polygon_layers_int["origindex"].data
        
        for i, p in enumerate(me.polygons):
            i2 = layer[i].value
            
            if i2 == -1: #i2 in fset:
                li = p.loop_start
                vs = []
                
                for j in range(p.loop_total):
                    vi = me.loops[li].vertex_index
                    v = me.vertices[vi]
                    vs.append(drawbm.verts.new(Vector(v.co)))
                    li += 1
                
                drawbm.faces.new(vs)
                
        matrix = ob2.matrix_world
        for v in drawbm.verts:
            v.co = matrix @ v.co
        
        drawbm.normal_update()
        
        c = [1, 1, 1, 1.0];
        
        for f in drawbm.faces:
            #c = colormap_get(0.9) #random.random()*0.15 + 0.15)
            
            if len(f.verts) == 3:
                v1, v2, v3 = f.verts
                
                t = sm.tri(v1.co, v2.co, v3.co)
                t.normals(v1.normal, v2.normal, v3.normal)
                t.colors(c, c, c);
                #t.uvs([0, 0], [0, 1], [1, 1])
            elif len(f.verts) == 4:
                v1, v2, v3, v4 = f.verts
                q = sm.quad(v1.co, v2.co, v3.co, v4.co)
                q.normals(v1.normal, v2.normal, v3.normal, v4.normal)
                q.colors(c, c, c, c);
                #q.uvs([0, 0], [0, 1], [1, 1], [1, 0])
            else:
                print("error; should have gotten subd surface with all quads");
                
        ob2.to_mesh_clear()
        
    def modal(self, context, event):
        coord = event.mouse_region_x, event.mouse_region_y
        
        dx = coord[0] - self.start_mpos[0]
        dy = coord[1] - self.start_mpos[1]
        #dl = sqrt(dx*dx + dy*dy) / 250
        #print(dl)
        self.factor = -dy / 250
        
        self._lines = [
            [self.start_mpos, coord]
        ]
        
        #print(event.type, dir(event), event.value, event.oskey, event.tilt)
        
        if event.type in {'MIDDLEMOUSE', 'WHEELUPMOUSE', 'WHEELDOWNMOUSE'}:
            context.area.tag_redraw()
            # allow navigation
            return {'PASS_THROUGH'}
        elif event.type == 'MOUSEMOVE':
            context.area.tag_redraw()
            self.main(context)
            return {'RUNNING_MODAL'}
        elif event.type == "LEFTMOUSE" or (event.type == "RET" and event.value != "RELEASE"):
            self.stop()
            return {'FINISHED'}
        elif event.type in {'RIGHTMOUSE', 'ESC'}:
            self.stop()
            context.area.tag_redraw()
            ob = context.active_object
            
            #can't rely on aborting earlier if this is false (in invoke) cause of dumb crash
            if ob is not None and ob.type == "MESH" and not ob.name.startswith("_"):
                bm = bmesh.from_edit_mesh(ob.data)
                for i, v in enumerate(bm.verts):
                    v.co = self.cos[i]
                bmesh.update_edit_mesh(ob.data)
            
            return {'CANCELLED'}

        return {'RUNNING_MODAL'}
    
    def stop(self):
        if hasattr(self, "_handle") and self._handle is not None:
            bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
            self._handle = None
            
        if hasattr(self, "_handle3d") and self._handle3d is not None:
            bpy.types.SpaceView3D.draw_handler_remove(self._handle3d, 'WINDOW')
            self._handle3d = None
        
        self._ob2 = None
        self._me = None
            
    def invoke(self, context, event):
        self._lines = []
        
        if context.space_data.type == 'VIEW_3D':
            self._ob2 = shadow.getUnSubShadow(context.active_object, ctx=context)
            
            self.makeDrawData(self._ob2)
            
            #print(event, dir(event))
            self.start_mpos = event.mouse_region_x, event.mouse_region_y
            self.mouse_path = []
            
            args = (self, context)
            self._handle = bpy.types.SpaceView3D.draw_handler_add(draw_callback_px, args, 'WINDOW', 'POST_PIXEL')
            self._handle3d = bpy.types.SpaceView3D.draw_handler_add(draw_callback_3dpx, args, 'WINDOW', 'POST_VIEW')
            self.text = ""
            
            context.window_manager.modal_handler_add(self)
            ob = context.active_object
            
            if ob.mode == "EDIT":
                ob.update_from_editmode()
            
            if ob is not None and ob.type == "MESH" and not ob.name.startswith("_"):
                #bpy.ops.ed.undo_push()
                
                self.cos = [Vector(v.co) for v in ob.data.vertices]
                self.nos = [Vector(v.normal) for v in ob.data.vertices]
                self.wos = [0.0 for v in ob.data.vertices]
                
                return {'RUNNING_MODAL'}
            else:
                return {'RUNNING_MODAL'} #XXX below line is crashing!
                #return {'CANCELLED'}
        else:
            self.report({'WARNING'}, "Active space must be a View3d")
            return {'CANCELLED'}

from . import simple_loop_optimize, solver
import imp
    
class FairLoopOperator(bpy.types.Operator):
    """UV Operator description"""
    bl_idname = "mesh.fair_loop"
    bl_label = "Fair Loop"
    bl_options = {'UNDO', 'REGISTER', 'PRESET'}
    
    factor = FloatProperty(name="factor", default=1.0)
    
    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj and obj.type == 'MESH' and obj.mode == 'EDIT'

    def execute(self, context):
      obj = context.active_object
      me = obj.data
      bm = bmesh.from_edit_mesh(me)
      
      imp.reload(solver)
      imp.reload(simple_loop_optimize)
      
      simple_loop_optimize.fairLoops(bm, obj.matrix_world, self.factor)
      
      bmesh.update_edit_mesh(me)
      return {'FINISHED'}

registrar = util.Registrar([
    ProjectToUnSubD,
    FairLoopOperator
]);

