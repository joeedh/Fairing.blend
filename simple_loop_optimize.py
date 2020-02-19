import bpy
import bmesh
from mathutils import *
from math import *
from bpy.props import *
from .solver import Constraint, Solver

def getLoop(v, visit):
  loop = [v]
  visit.add(v.index)
  
  startv = v
  
  while 1:
    found = False
    
    for e in v.link_edges:
      v2 = e.other_vert(v)
      if v2.hide or not v2.select or v2.index in visit: continue
      
      visit.add(v2.index)
      loop.append(v2)
      
      v = v2
      found = True
      break
    
    if not found:
      break
    
    if v == startv:
      break
      
  return loop
  
def safe_acos(f):
  f = min(max(f, -0.99999), 0.999999)
  return acos(f)
  
def fairLoop(bm, matrix, loop, closed, factor=1.0):
  print("FACTOR", factor)
  totlen = 0
  cos = [matrix @ v.co for v in loop]
  masses = []
  
  for i in range(len(loop)):
    if not closed and (i == 0 or i == len(loop)-1):
      masses.append(100000)
    else:
      masses.append(1)
  
  knots = []
  
  for i in range(len(loop)-1):
    v1 = cos[i]
    v2 = cos[i+1]
    
    l = (v2 - v1).length
    knots.append(totlen)
    
    totlen += l
  
  knots.append(knots[len(knots)-1])
  
  origcos = [co.copy() for co in cos]
  
  def evaluate(s):
    start = 0
    end = len(loop)-1
    mid = (start + end) // 2
    
    while end - start > 1:
      if s > knots[mid]:
        start = mid
      elif s < knots[mid]:
        end = mid
      elif s == knots[mid]:
        break
      
      mid = (start + end) // 2
    
    #print(s, knots[mid], knots[mid+1], totlen)
    ks = knots[mid]
    ds = s - ks
    
    if (mid < len(knots)-1):
      dmid = abs(knots[mid+1] - knots[mid])
    else:
      dmid = 0.0
      
    #print(mid, s, knots[mid], ds, dmid)
    #return origcos[mid]
    
    if ds > 0.0 and mid < len(knots)-1 and dmid > 0.0:
      t = ds / dmid
      
      v1 = origcos[mid]
      v2 = origcos[mid+1]
      
      return v1 + (v2 - v1)*t
    else:
      return origcos[mid]
    
    return origcos[mid]
  
  bpy._evaluate_test = evaluate
  #return
  
  #temporarily make verts evenly spaced
  for i in range(len(cos)):
    if not closed and (i == 0 or i == len(cos)-1):
      continue
      
    #print("s", i * totlen / len(cos), totlen)
    cos[i] = evaluate(float(i) * totlen / (len(cos)))
    
    
  def goal_c(params):
    v1, goal = params
    
    return (v1 - goal).length
    
  def len_c(params):
    v1, v2, restlen = params
    
    return (v2 - v1).length - restlen
  
  def curv_c(params):
    v0, v1, v2, v3, v4 = params
    
    l1 = (v0 - v1).length + 0.00001
    l2 = (v1 - v2).length + 0.00001
    l3 = (v3 - v2).length + 0.00001
    l4 = (v4 - v3).length + 0.00001
    
    t1 = (v1 - v0) / l1
    t2 = (v2 - v1) / l2

    t3 = (v2 - v1) / l2
    t4 = (v3 - v2) / l3

    t5 = (v3 - v2) / l3
    t6 = (v4 - v3) / l4
    
    d1 = t2 - t1
    d2 = t4 - t3
    d3 = t6 - t5
    
    bi1 = t1.cross(d1)
    bi2 = t3.cross(d2)
    bi3 = t5.cross(d3)
    
    d4 = d3 - d1
    
    return d4.length

  def bitan_c(params):
    v0, v1, v2, v3, v4 = params
    
    l1 = (v0 - v1).length + 0.00001
    l2 = (v1 - v2).length + 0.00001
    l3 = (v3 - v2).length + 0.00001
    l4 = (v4 - v3).length + 0.00001
    
    t1 = (v1 - v0) / l1
    t2 = (v2 - v1) / l2

    t3 = (v2 - v1) / l2
    t4 = (v3 - v2) / l3

    t5 = (v3 - v2) / l3
    t6 = (v4 - v3) / l4
    
    d1 = t2 - t1
    d2 = t4 - t3
    d3 = t6 - t5
    
    d1.normalize()
    d2.normalize()
    d3.normalize()
    
    bi1 = t1.cross(d1)
    bi2 = t3.cross(d2)
    bi3 = t5.cross(d3)
    
    #bi1.normalize()
    #bi2.normalize()
    #bi3.normalize()
    
    return (bi3 - bi1).length
    
  def bend_c(params):
    v0, v1, v2 = params
    
    t1 = v1 - v0
    t2 = v2 - v1
    
    t1.normalize()
    t2.normalize()
    
    return safe_acos(t1.dot(t2))
    
    
  restlen = totlen / len(loop)
  slv = Solver()
  
  for i, v in enumerate(loop):
    #if (i == 0 or i == len(loop)-1) and not closed:
    #  continue
    
    #"""
    if closed or i > 1 or i < len(loop)-2:
      i0 = (i - 2 + len(loop)) % len(loop)
      i1 = (i - 1 + len(loop)) % len(loop)
      i2 = i
      i3 = (i + 1) % len(loop)
      i4 = (i + 2) % len(loop)
      
      con = Constraint(
        "curv_c", 
        [cos[i0], cos[i1], cos[i2], cos[i3], cos[i4]], 
        [masses[i0], masses[i1], masses[i2], masses[i3], masses[i4]], 
        curv_c, 
        [cos[i0], cos[i1], cos[i2], cos[i3], cos[i4]]
      )
      con.k = 0.5
      slv.add(con)

      con = Constraint(
        "bitangent", 
        [cos[i0], cos[i1], cos[i2], cos[i3], cos[i4]], 
        [masses[i0], masses[i1], masses[i2], masses[i3], masses[i4]], 
        bitan_c, 
        [cos[i0], cos[i1], cos[i2], cos[i3], cos[i4]]
      )
      con.k = 0.5
      slv.add(con)
      
      #"""
      con = Constraint(
        "bend_c", 
        [cos[i1], cos[i2], cos[i3]], 
        [masses[i1], masses[i2], masses[i3]], 
        bend_c, 
        [cos[i1], cos[i2], cos[i3]]
      )
      con.k = 0.1
      #slv.add(con)
      #"""
    
    #Constraint(name, klst, mlst, func, params, wlst=None):
    v2 = loop[(i+1) % len(loop)]
    i2 = (i + 1) % len(loop)
    
    #restlen = (cos[i2] - cos[i]).length
    restlen = totlen / len(loop)/2.0
    
    con = Constraint("goal_c", [cos[i]], [masses[i]], goal_c, [cos[i], Vector(cos[i])])
    con.k = 0.01
    slv.add(con)
    
    con = Constraint("len_c", [cos[i], cos[i2]], [masses[i], masses[i2]], len_c, [cos[i], cos[i2], restlen])
    slv.add(con)

  slv.solve(55)
  
  mat = Matrix(matrix)
  mat.invert()
  
  for i in range(len(knots)):
    knots[i] /= totlen
    
  totlen = 0;
  for i in range(len(cos)-1):
    v1 = cos[i]
    v2 = cos[i+1]
    totlen += (v2 - v1).length
  
  for i in range(len(knots)):
    knots[i] *= totlen
    
  origcos = cos
  for i in range(len(cos)):
    co = Vector(cos[i])
    #co = evaluate(i * totlen / len(cos)) #knots[i])
    co = evaluate(knots[i])
    
    if not closed and (i == 0 or i == len(loop)-1):
      continue
      
    co = mat @ co
    loop[i].co += (co - loop[i].co) * factor
    
def fairLoops(bm, matrix, factor=1.0):
  bm.verts.index_update()
  
  loops = []
  visit = set()
  
  #first find non-closed loops
  for v in bm.verts:
    if v.hide or not v.select or v.index in visit: continue
    
    tot = 0
    for e in v.link_edges:
      v2 = e.other_vert(v)
      
      if v2.select and not v2.hide and v2.index not in visit:
        tot += 1
    
    if tot == 1:
      loop = getLoop(v, visit)
      loops.append([loop, False])
    
  for l in loops:
    print("found loop:", len(l))
    fairLoop(bm, matrix, l[0], l[1], factor)
