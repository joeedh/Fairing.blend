from . import globals
import bpy

OLD_STYLE = True

def clearTimers():
    for t in globals.timers[:]:
        if not OLD_STYLE and bpy.app.timers.is_registered(t):
            bpy.app.timers.unregister(t)
        globals.timers.remove(t)
    
def addTimer(func, first_ival=0.01):
    if not globals._modal_timer_running:
        globals._modal_timer_running = True
        bpy.ops.wm._fairing2_util_timer()
        
    globals.timers.append(func)
    if not OLD_STYLE:
        bpy.app.timers.register(func, first_interval=first_ival)

    
def remTimer(func):
    if not OLD_STYLE and bpy.app.timers.is_registered(func):
        bpy.app.timers.unregister(func)
    
    if func in globals:
        globals.timers.remove(func)
