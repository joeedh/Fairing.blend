import bpy
import time

from . import util

class ModalTimerOperator(bpy.types.Operator):
    """Operator which runs its self from a timer"""
    bl_idname = "wm._fairing2_util_timer"
    bl_label = "my timer runner"

    _timer = None
    last_time = 0
    
    def modal(self, context, event):
        from . import globals
        
        if event.type in {'ESC'} or not globals._modal_timer_running:
            self.cancel(context)
            print(globals._modal_timer_running)
            globals._modal_timer_running = False
            return {'CANCELLED'}

        if event.type == 'TIMER' and time.time() - self.last_time > 0.5:
            for t in globals.timers[:]:
                ret = t()
                if ret is None:
                    globals.timers.remove(t)
                    
            self.last_time = time.time()
            

        return {'PASS_THROUGH'}

    def execute(self, context):
        from . import globals

        globals._modal_timer_running = True
        
        wm = context.window_manager
        self._timer = wm.event_timer_add(0.1, window=context.window)
        wm.modal_handler_add(self)
        
        return {'RUNNING_MODAL'}

    def cancel(self, context):
        from . import globals

        globals._modal_timer_running = False
        
        wm = context.window_manager
        wm.event_timer_remove(self._timer)

from . import globals
globals._modal_timer_running = False

registrar = util.Registrar([
    ModalTimerOperator
]);
