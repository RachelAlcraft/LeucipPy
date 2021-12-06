try:
    from LeucipPy import Interpolator
except:
    import Interpolator

'''
Give it some data and get interp of given degree and dimensions
Give it 1,2 or 3d data and get back a multivariate equation
'''

class BSpline(Interpolator):
    def __init__(self):
        self.Method = "BS"

    def getValue(self):
        return 0