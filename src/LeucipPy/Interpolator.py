
try:
    from LeucipPy import Interpolator
except:
    import Interpolator

'''
Abstract class so we can use the same functions
'''

class Interpolator:
    def __init__(self):
        self.Method = "Generic"

    def getValue(self):
        return 0  #use whatever data has been loaded into the class and generically get value so we can differentiate