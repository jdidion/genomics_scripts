try:
    from progressbar import ProgressBar, Widget, Percentage, Bar, ETA
    available = True
except ImportError:
    available = False
    
class MyProgressBar(ProgressBar):
    def increment(self, n=1):
        self.update(self.currval + n)

class Status(Widget):
    def __init__(self, value=None):
        self.value = value
        
    def update(self, pbar):
        return self.value or ""

def default_progress(name=None, maxval=None, status=False):
    widgets = []
    if status:
        status = Status()
    if name and status:
        widgets.extend((name, " (", status, "): "))
    elif name:
        widgets.extend((name, ": "))
    elif status:
        widgets.extend((status, ": "))
    widgets.extend((Percentage(), " ", Bar(), " ", ETA()))
    
    progress = MyProgressBar(maxval=maxval, widgets=widgets)
    progress.status = status
    return progress
