bl_info = {
    "name": "Fairing Tools",
    "author": "Joseph Eagar",
    "version": (0, 1),
    "blender": (2, 80, 0),
    "description": "Tools for fairing meshes",
    "warning": "",
    "wiki_url": "",
    "category": "Mesh",
}

import imp
 
__all__ = [
    "globals", "util", "ops", "props", "ui", "timer", "timer_oldstyle",
    "shadow", "simplemesh", "subsurf", "sym", "smatrix", "subsurf_evaluate",
    "bezpatch", "unsubdiv", "optimize", "simple_loop_optimize", "solver"
]

from . import globals, util, ops, props, ui, timer,  \
              timer_oldstyle, shadow, simplemesh, \
              subsurf, sym, smatrix, bezpatch, subsurf_evaluate, \
              unsubdiv, optimize, simple_loop_optimize, solver

timer.clearTimers()

props.registrar.unreg()
ops.registrar.unreg()
ui.registrar.unreg()
timer_oldstyle.registrar.unreg()

imp.reload(util)
imp.reload(sym)
imp.reload(solver)
imp.reload(simplemesh)
imp.reload(timer_oldstyle)
imp.reload(timer)
imp.reload(bezpatch)
imp.reload(smatrix)
imp.reload(subsurf)
imp.reload(subsurf_evaluate)
imp.reload(unsubdiv)
imp.reload(optimize)
imp.reload(props)
imp.reload(shadow)
imp.reload(ui)
imp.reload(simple_loop_optimize)
imp.reload(ops)

registrar = util.Registrar([
    ops.registrar,
    props.registrar,
    ui.registrar,
    timer_oldstyle.registrar
])

def register():
    registrar.reg()

def unregister():
    registrar.unreg()

