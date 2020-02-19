from mathutils import *
from bgl import *
import gpu
from gpu_extras.batch import batch_for_shader
shader = gpu.shader.from_builtin('3D_UNIFORM_COLOR')

vertex_shader = '''
    uniform mat4 viewProjectionMatrix;
    uniform vec4 uColor;
    
    in vec3 position;
    in vec2 uv;
    in vec4 color;
    
    out vec3 pos;
    out vec2 in_uv;
   
    void main()
    {
        pos = position;
        in_uv = uv;
        gl_Position = vec4(position, 1.0f);
    }
'''

fragment_shader = '''
    uniform float brightness;
    uniform vec4 uColor;

    in vec3 pos;
    in vec2 in_uv;
    
    void main()
    {
        vec4 c = uColor;
        
        c[0] = in_uv[0];
        
        gl_FragColor = c;
    }
'''

SimpleShader = gpu.types.GPUShader(vertex_shader, fragment_shader);


vertex_shader3d = '''
    uniform mat4 viewProjectionMatrix;
    uniform vec4 uColor;
    uniform float uPolyOffset;
    
    in vec3 position;
    in vec2 uv;
    in vec3 normal;
    in vec4 color;
    
    centroid out vec3 pos;
    centroid out vec2 f_uv;
    centroid out vec3 f_normal;
    centroid out vec3 f_color;
    
    void main()
    {
        pos = position;
        f_uv = uv;
        f_normal = normalize((viewProjectionMatrix * vec4(normal, 0.0f)).xyz);
        f_color = color.rgb;
        
        vec4 p = viewProjectionMatrix * vec4(position, 1.0f);
        p[2] -= uPolyOffset / p[3];
        
        gl_Position = p;
    }
'''

fragment_shader3d = '''
    uniform float brightness;
    uniform vec4 uColor;
    
    centroid in vec3 pos;
    centroid in vec2 f_uv;
    centroid in vec3 f_normal;
    centroid in vec3 f_color;
    
    void main()
    {
        vec4 c = uColor * vec4(f_color, 1.0);
        
        //dumb line to keep uv from getting optimized out, stupid gpu
        //shader system is buggy that way and apparently can't live without
        //uv data, bleh
        
        c *= f_uv[0]*0.01 + c[0]*0.99;
        
        vec3 normal = normalize(f_normal);
        
        float f = pow(abs(normal[2]), 3.0);
        float f2 = pow(abs(normal[1]), 3.0);
        float f3 = pow(abs(normal[0]), 3.0);
        
        f = f*0.7 + f2*0.05 + f3*0.25;
        
        c.rgb *= f*1.4;
        
        gl_FragColor = c;
    }
'''
SimpleShader3D = gpu.types.GPUShader(vertex_shader3d, fragment_shader3d);

#class Shader:
#    def __init__(self, v, f, attrs, uniforms):
#        pass
        
class LayerTypes:
    VERTEX = 1
    UVS = 2
    COLORS = 4
    NORMALS = 8
    IDS = 16

LayerTypeList = [1, 2, 4, 8, 16]

ElemSizes = {}
ElemSizes[LayerTypes.VERTEX] = 3
ElemSizes[LayerTypes.UVS] = 2
ElemSizes[LayerTypes.COLORS] = 4
ElemSizes[LayerTypes.NORMALS] = 3
ElemSizes[LayerTypes.IDS] = 1

AttrNames = {}
AttrNames[LayerTypes.VERTEX] = "position"
AttrNames[LayerTypes.UVS] = "uv"
AttrNames[LayerTypes.COLORS] = "color"
AttrNames[LayerTypes.NORMALS] = "normal"
AttrNames[LayerTypes.IDS] = "id"

DefaultLayers = LayerTypes.VERTEX #| LayerTypes.COLORS | LayerTypes.UVS

class VBuffer:
    def __init__(self, type):
        self.data = []
        self.type = type
        self.esize = ElemSizes[type]
    
    def gldestroy(self):
        pass
    
    def getBuffer(self):
        buf = bgl.Buffer(bgl.FLOAT32, [len(self.data)])
        for i, d in enumerate(self.data):
            buf[i] = d
        
        self._buf = buf
        return buf
        
class BufferList (list):
    def __init__(self):
        self.map = {}
        self.totelem = 0
    
    def gldestroy(self):
        for b in self:
            b.gldestroy()
    
    def newElems(self, count):
        ret = self.newElem()
        
        for i in range(count-1):
            self.newElem()
            
        return ret
        
    def newElem(self):
        i = self.totelem
        
        for buf in self:
            for j in range(buf.esize):
                buf.data.append(0)
                
        self.totelem += 1
        return i
        
    def create(self, type):
        buf = VBuffer(type)
        self.append(buf)
        self.map[type] = buf
        return buf
        
    def get(self, type):
        #if type not in self.map:
        #    return self.create(type)
        
        if type not in self.map:
            raise RuntimeError("buffer " + str(type) + " not in map")
            
        return self.map[type]
    
    def has(self, type):
        return type in self.map
        
class TriEditor:
    def __init__(self, sm):
        self.sm = sm
        self.i = sm.tribufs.newElems(3)
        
    def vs(self, v1, v2, v3):
        sm = self.sm
        i = self.i*3
        
        buf = sm.tribufs.get(LayerTypes.VERTEX)
        d = buf.data
        
        d[i]   = v1[0]; d[i+1] = v1[1]; d[i+2] = v1[2];
        d[i+3] = v2[0]; d[i+4] = v2[1]; d[i+5] = v2[2];
        d[i+6] = v3[0]; d[i+7] = v3[1]; d[i+8] = v3[2];
        
        return self
            
    def uvs(self, u1, u2, u3):
        sm = self.sm
        i = self.i*2
        buf = sm.tribufs.get(LayerTypes.UVS)
        d = buf.data
        
        d[i]   = u1[0]; d[i+1] = u1[1];
        d[i+2] = u2[0]; d[i+3] = u2[1];
        d[i+4] = u3[0]; d[i+5] = u3[1];
        
        return self
        
    def normals(self, v1, v2, v3):
        sm = self.sm
        i = self.i*3
        
        buf = sm.tribufs.get(LayerTypes.NORMALS)
        d = buf.data
        
        d[i]   = v1[0]; d[i+1] = v1[1]; d[i+2] = v1[2];
        d[i+3] = v2[0]; d[i+4] = v2[1]; d[i+5] = v2[2];
        d[i+6] = v3[0]; d[i+7] = v3[1]; d[i+8] = v3[2];
        
        return self
        
    def colors(self, c1, c2, c3):
        sm = self.sm
        i = self.i*4
        
        buf = sm.tribufs.get(LayerTypes.COLORS)
        d = buf.data
        
        d[i] = c1[0]; d[i+1] = c1[1]; d[i+2] = c1[2]; d[i+3] = c1[3];
        d[i+4] = c2[0]; d[i+5] = c2[1]; d[i+6] = c2[2]; d[i+7] = c2[3];
        d[i+8] = c3[0]; d[i+9] = c3[1]; d[i+10] = c3[2]; d[i+11] = c3[3];
        
        return self

class QuadEditor:
    def __init__(self, sm, v1, v2, v3, v4):
        self.t1 = sm.tri(v1, v2, v3)
        self.t2 = sm.tri(v1, v3, v4)
    
    def vs(self, v1, v2, v3, v4):
        self.t1.vs(v1, v2, v3)
        self.t2.vs(v1, v3, v4)
        return self
        
    def uvs(self, v1, v2, v3, v4):
        self.t1.uvs(v1, v2, v3)
        self.t2.uvs(v1, v3, v4)
        return self
        
    def normals(self, v1, v2, v3, v4):
        self.t1.normals(v1, v2, v3)
        self.t2.normals(v1, v3, v4)
        return self
        
    def colors(self, v1, v2, v3, v4):
        self.t1.colors(v1, v2, v3)
        self.t2.colors(v1, v3, v4)
        return self
        
class LineEditor:
    def __init__(self, sm):
        self.sm = sm
        self.i = sm.linebufs.newElems(2)
        
    def vs(self, v1, v2):
        sm = self.sm
        i = self.i*3
        
        buf = sm.linebufs.get(LayerTypes.VERTEX)
        d = buf.data
        
        d[i]   = v1[0]; d[i+1] = v1[1]; d[i+2] = v1[2];
        d[i+3] = v2[0]; d[i+4] = v2[1]; d[i+5] = v2[2];
        
        return self
            
    def uvs(self, u1, u2):
        sm = self.sm
        i = self.i*2
        buf = sm.linebufs.get(LayerTypes.UVS)
        d = buf.data
        
        d[i]   = v1[0]; d[i+1] = v1[1];
        d[i+2] = v2[0]; d[i+3] = v2[1];
        
        return self
    
class SimpleMesh:
    def __init__(self, layerflag=None, shader=None):
        self.layerflag = DefaultLayers if layerflag is None else layerflag
        shader = SimpleShader if shader is None else shader
        
        self.shader = shader
        self.tribufs = BufferList()
        self.linebufs = BufferList()
        self._regen = True
        self.tribatch = None
        self.linebatch = None
        
        for type in LayerTypeList:
            if self.layerflag & type:
                self.tribufs.create(type)
                self.linebufs.create(type)
    
    def tri(self, v1, v2, v3):
        self._regen = True
        return TriEditor(self).vs(v1, v2, v3)
    
    def quad(self, v1, v2, v3, v4):
        return QuadEditor(self, v1, v2, v3, v4)
        
    def line(self, v1, v2):
        self._regen = True
        return LineEditor(self).vs(v1, v2)
        
    def regen(self):
        shader = self.shader
        
        self.linebatch = None
        self.tribatch = None
        
        bufs2 = {}
        for step in range(2):
            bufs = self.linebufs if step else self.tribufs
            if bufs.totelem == 0:
                continue
                
            for b in bufs:
                data = []
                name = AttrNames[b.type]
                bufs2[name] = data
                
                #print("L::", len(b.data) / b.esize, b.esize)
                
                for i in range(len(b.data) // b.esize):
                    l = []
                    for j in range(b.esize):
                        l.append(b.data[i*b.esize + j])
                    l = tuple(l)
                    
                    data.append(l)
                
            #attr = shader.attr_from_name(AttrNames[b.type])
            #b.getBuffer()
            #b.glid = attr
            
            batch = batch_for_shader(shader, 'TRIS' if not step else 'LINES', bufs2)
            if step:
                self.linebatch = batch
            else:
                self.tribatch = batch
        
                
    def destroy():
        pass
        
    def draw(self, uniforms={}):
        if self._regen:
            self._regen = False
            self.regen()
        
        if self.tribatch is None and self.linebatch is None:
            #print("failed to generate gl data")
            return
            
        shader = self.shader
        
        shader.bind()
        for k in uniforms.keys():
            v = uniforms[k]
            if type(v) == float:
                shader.uniform_float(k, v)
            elif type(v) == int:
                shader.uniform_int(k, v)
            elif type(v) in [Vector, list, Matrix]:
                shader.uniform_float(k, v)
        
        if self.tribatch is not None:
            self.tribatch.draw(shader)
        
        if self.linebatch is not None:
            self.linebatch.draw(shader)
        #matrix = bpy.context.region_data.perspective_matrix
        #shader.uniform_float("viewProjectionMatrix", matrix)
        
        
        
        
        
        
        
        
        
        
        
        
        
        