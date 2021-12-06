try:
    from LeucipPy import Interpolator
except:
    import Interpolator

'''
Give it some data and get interp of given degree and dimensions
Give it 1,2 or 3d data and get back a multivariate equation
'''

class Multivariate(Interpolator):
    def __init__(self):
        self.Method = "MV"

    def getValue(self):
        return 0

    def getMultivariateDescription(self,npdata):
        return ""
