import bpy
from ctypes import c_double
from ctypes import addressof
from ctypes import *

#degree raising/lowering
"""
on factor;
off period;

procedure bez(a, b);
    a + (b - a)*s;

clear lin, quad, cubic;

comment: degree raising;

lin := bez(k1, k2);
quad := bez(lin, sub(k2=k3, k1=k2, lin));
cubic := bez(quad, sub(k3=k4, k2=k3, k1=k2, quad));

lin := sub(k1=oldk1, k2=oldk2, lin);

f1 := sub(s=0, lin) - sub(s=0, quad);
f2 := sub(s=1, lin) - sub(s=1, quad);
f3 := sub(s=0.5, lin) - sub(s=0.5, quad);

deg0to1 := list(k1 = oldk1, k2 = oldk1);
deg1to2 := part(solve({f1, f2, f3}, {k1, k2, k3}), 1);

quad := sub(k1=oldk1, k2=oldk2, k3=oldk3, quad);

f1 := sub(s=0, quad) - sub(s=0, cubic);
f2 := sub(s=1, quad) - sub(s=1, cubic);
f3 := sub(s=0, df(quad, s)) - sub(s=0, df(cubic, s));
f4 := sub(s=1, df(quad, s)) - sub(s=1, df(cubic, s));

deg2to3 := part(solve({f1, f2, f3, f4}, {k1, k2, k3, k4}), 1);

comment: degree lowering;

clear lin, quad, cubic;

lin := bez(k1, k2);
quad := bez(lin, sub(k2=k3, k1=k2, lin));
cubic := bez(quad, sub(k3=k4, k2=k3, k1=k2, quad));

cubic := sub(k1=oldk1, k2=oldk2, k3=oldk3, k4=oldk4, cubic);

f1 := sub(s=0, quad) - sub(s=0, cubic);
f2 := sub(s=1, quad) - sub(s=1, cubic);
comment: f3 := sub(s=0.5, quad) - sub(s=0.5, cubic);

f3 := sub(s=0.0, df(quad, s)) - sub(s=0.0, df(cubic, s));
deg3to2a := solve({f1, f2, f3}, {k1, k2, k3});

f3 := sub(s=1.0, df(quad, s)) - sub(s=1.0, df(cubic, s));
deg3to2b := solve({f1, f2, f3}, {k1, k2, k3});

deg3to2 := list(oldk1, part(deg3to2a, 1, 2)*0.5 + part(deg3to2b, 1, 2)*0.5, oldk4);

quad := sub(k1=oldk1, k2=oldk2, k3=oldk3, quad);

deg2to1 := list(k1 = oldk1, k2 = oldk3);
deg1to0 := list(k1 = (oldk1+oldk2)*0.5);

comment: degree 0 to 1;
part(deg0to1, 1);
part(deg0to1, 2);

comment: degree 1 to 2;
part(deg1to2, 1);
part(deg1to2, 2);
part(deg1to2, 3);

comment: degree 2 to 3;
part(deg2to3, 1);
part(deg2to3, 2);
part(deg2to3, 3);
part(deg2to3, 4);

comment: degree 3 to 2;
part(deg3to2, 1);
part(deg3to2, 2);
part(deg3to2, 3);

comment: degree 2 to 1;
part(deg2to1, 1);
part(deg2to1, 2);

comment: degree 1 to 0;
part(deg1to0, 1);
"""    

def bez_raise(bz):
    if len(bz) == 1:
        oldk1 = bz[0]
        return [
            oldk1, oldk1
        ]
    elif len(bz) == 2:
        oldk1 = bz[0]
        oldk2 = bz[1]
        return [
            oldk1, (oldk1+oldk2)*0.5, oldk2
        ]
    elif len(bz) == 3:
        oldk1 = bz[0]
        oldk2 = bz[1]
        oldk3 = bz[2]
        return [
            oldk1, (oldk1+2*oldk2)/3, (2*oldk2+oldk3)/3, oldk3
        ]
    
    raise RuntimeError("unsupported degree " + str(len(bz)-1) + " (" + str(len(bz)) + " points)")
    
def bez_lower(bz):
    if len(bz) == 2:
        oldk1 = bz[0]
        oldk2 = bz[1]
        
        return [
            (oldk1 + oldk2)*0.5
        ]
    elif len(bz) == 3:
        oldk1 = bz[0]
        oldk2 = bz[1]
        oldk3 = bz[2]
        
        return [
            oldk1, oldk3
        ]
    elif len(bz) == 4:
        oldk1 = bz[0]
        oldk2 = bz[1]
        oldk3 = bz[2]
        oldk4 = bz[3]
        
        return [
            oldk1, 
            (3*oldk3-oldk4+3*oldk2-oldk1)/4,
            oldk4
        ]
        
    raise RuntimeError("unsupported degree " + str(len(bz)-1) + " (" + str(len(bz)) + " points)")
        
def degree_change(bez, n):
    while len(bez) > n+1:
        bez = bez_lower(bez)
    while len(bez) < n+1:
        bez = bez_raise(bez)
    return bez
    
def sortlist(lst, cb):
    class Item:
        def __init__(self, item):
            self.item = item
        
        def __eq__(self, b):
            return cb(self.item, b.item) == 0
        def __lt__(self, b):
            return cb(self.item, b.item) < 0
        def __gt__(self, b):
            return cb(self.item, b.item) > 0
            
    for i in range(len(lst)):
        lst[i] = Item(lst[i])
    lst.sort()
    
    for i in range(len(lst)):
        lst[i] = lst[i].item

_bidgen = 0

def _buf_id(buf):
    global _bidgen
    
    buf._bid = _bidgen
    _bidgen += 1
    
    def hash():
        return buf._bid
        
    buf.__hash__ = hash
        
def makeBuffer(size):
    ret = (c_double*size)()
    _buf_id(ret);
    
    return ret
    
def sliceBuffer(buf, a, b):
    rettype = POINTER(c_double*(b-a))
    
    ptr2 = cast(buf, c_void_p)
    ptr2.value += sizeof(c_double)*a
    
    ret = cast((ptr2), rettype)[0]
    _buf_id(ret)
    
    return ret 

def testBuffers():
    a = makeBuffer(5)
    a[0] = 1.1
    a[1] = 4.7
    a[3] = -2.2
    a[4] = 88.234
    
    b = sliceBuffer(a, 1, 3)
    
    print(len(a), len(b))
    
    b[1] = 55.0;
    print(b[1], a[2])
    
    a[2] = 0.5;
    print(b[1], a[2])
    
    
class cachering (list):
    def __init__(self, size, cls):
        #list.__init__(self)
        
        self.cur = 0
        for i in range(size):
            self.append(cls())
    
    def next(self):
        ret = self[self.cur]
        
        self.cur = (self.cur + 1) % len(self)
        return ret

class Reg:
    def reg(self):
        pass
    def unreg(self):
        pass
        
class CustomReg (Reg):
    def __init__(self, regcb, unregcb, data):
        self.rcb = rcb
        self.ucb = ucb
        self.data = data
    
    def reg(self):
        self.rcb(self.data)
        
    def unreg(self):
        self.ucb(self.data)
        
class BpyReg (Reg):
    def __init__(self, data):
        self.data = data
        
    def reg(self):
        bpy.utils.register_class(self.data)
    
    def unreg(self):
        bpy.utils.unregister_class(self.data)
        
class Registrar (list, Reg):
    def __init__(self, classes):
        for c in classes:
            if not isinstance(c, Reg):
                c = BpyReg(c)
            self.append(c)
            
        self.registered = False
    
    def reg(self):
        if self.registered:
            self.unreg()
            
        for cls in self:
            cls.reg()
        
        self.registered = True
        
    def unreg(self):
        if not self.registered:
            return
        
        for cls in self:
            cls.unreg()
            
        self.registered = False
    

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
    