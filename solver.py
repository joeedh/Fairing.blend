from math import *

class Constraint:
  def __init__(self, name, klst, mlst, func, params, wlst=None):
    self.name = name
    self.k = 1.0
    
    if wlst is None:
      wlst = []
      for k in klst:
        wlst.append([1 for x in range(len(k))])
    
    glst = []
    for k in klst:
      glst.append([0 for x in range(len(k))])
    
    self.mlst = mlst
    self.df = 0.0001
    self.threshold = 0.00001
    self.func = func
    self.params = params
    
    self.wlst = wlst
    self.glst = glst
    self.klst = klst
    
  def evaluate(self, calc_dvs=True):
    r1 = self.func(self.params)
    
    if abs(r1) < self.threshold: return 0.0
    
    df = self.df
    for i in range(len(self.klst)):
      ks = self.klst[i]
      gs = self.glst[i]
      
      for j in range(len(ks)):
        orig = ks[j]
        ks[j] += df
        r2 = self.func(self.params)
        
        gs[j] = (r2 - r1) / df
        
        ks[j] = orig
    
    return r1

class Solver:
  def __init__(self):
    self.cons = [];
    self.gk = 1.0
    self.threshold = 0.0005
  
  def solve(self, steps=5, gk=1.0):
    error = 0.0
    
    for i in range(steps):
      error = self.solve_intern(gk)
        
      print(i+1, "%.5f" % error)
      if error < self.threshold:
        break
      
    return error
  
  def add(self, con):
    self.cons.append(con)
    
  def solve_intern(self, gk=1.0):
    gk *= self.gk
    error = 0.0
    
    for con in self.cons:
      r1 = con.evaluate()
      
      if abs(r1) == 0.0: continue
      
      error += abs(r1)
      
      totg = 0.0
      totmass = 0.0
      
      for i, ks in enumerate(con.klst):
        gs = con.glst[i]
        
        for g in gs:
          totg += g*g
          
        for m in con.mlst:
          #totmass += 1.0 / m
          totmass = max(1.0 / m, totmass)
          
      totmass = 1.0 / totmass
      
      if totg == 0.0: continue
      
      r1 /= totg
      k = con.k * gk
      
      for i, ks in enumerate(con.klst):
        gs = con.glst[i]
        ws = con.wlst[i]
        ms = con.mlst[i]
        mw = (1.0 / ms)
        
        for j in range(len(ks)):
          ks[j] += -gs[j]*ws[j]*k*r1*mw
          
    return error
      


  