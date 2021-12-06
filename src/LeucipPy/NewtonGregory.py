
try:
    from LeucipPy import Interpolator
except:
    import Interpolator
'''
Give it some data and get poly, bipoly or tripoly interpolation back
Give it a sequence and get back the polynomial
'''


class NewtonGregory(Interpolator):
    def __init__(self):
        self.Method = "NG"

    def getValue(self):
        return 0

    def getPolynomialDescription(self,sequence):
        return ""
