ops = set([
    "const", "+", "-", "*", "/", "**", "&", "|", "&&", "||", "==", ">=", "<=", "var"
])

prec = {
    "**" : 0,
    "*" : 1,
    "/" : 1,
    "+" : 2,
    "-" : 2
}

class sym:
    def __init__(self, a=0, b=None, op=None):
        self.parent = None
        self.prec = -1
        
        if isinstance(a, sym) and b is None:
            self.load(a)
            return
        elif type(a) == str and b is None and (op is None or op == "var"):
            self.op = "var"
            self.a = a
            self.b = None
            return
        elif type(a) in [int, float] and b is None and (op is None or op == "const"):
            self.op = "const"
            self.a = a
            self.b = None
            return
        
        if not isinstance(a, sym):
            a = sym(a)
            a.parent = self
            
        if b is not None and not isinstance(b, sym):
            b = sym(b)
            b.parent = self
            
        self.prec = -1 if op is None or op not in prec else prec[op]
        self.op = op
        
        self.a = a
        self.b = b
        
        self.a.parent = self
        if self.b is not None:
            self.b.parent = self
        
    def __str__(self):
        if self.op == "var":
            return self.a
        elif self.op == "const":
            return str(self.a)
        else:
            ret = str(self.a) + " " + self.op + " " + str(self.b)
            if self.parent is not None and self.prec > self.parent.prec:
                ret = "(" + ret + ")"
                
            return ret
        
    def __repr__(self):
        return str(self)
        
    def load(self, b):
        self.op = b.op
        self.prec = b.prec
        
        if b.op == "const" or b.op == "var":
            self.a = b.a
            self.b = b.b
        else:
            #print(b, b.op, type(b.a), type(b.b))
            self.a = b.a.copy()
            self.b = b.b.copy()
            
            self.a.parent = self
            self.b.parent = self
            
        return self
        
    def copy(self):
        s = sym()
        s.load(self)
        return s
    
    def binop(self, b, op):
        ret = sym(a=self, b=sym(b), op=op)
        self.parent = ret
        
        ret.a.parent = ret
        ret.b.parent = ret
        
        return ret
        
    def __add__(self, b):
        return self.binop(b, "+")
    def __radd__(self, b):
        return sym(b).binop(self, "+")
    
    def __sub__(self, b):
        return self.binop(b, "-")
    def __rsub__(self, b):
        return sym(b).binop(self, "-")

    def __mul__(self, b):
        return self.binop(b, "*")
    def __rmul__(self, b):
        return sym(b).binop(self, "*")

    def __truediv__(self, b):
        return self.binop(b, "/")
    def __rtruediv__(self, b):
        return sym(b).binop(self, "/")
    
    def simp(self):
        return self
        