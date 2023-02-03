
import math


def getMagnitude(x,y,z):
    return math.sqrt(x*x + y*y + z*z)

def getDistance(x1, y1, z1, x2, y2, z2):
    """
    param: x1,y1,z1= the coordinates of the first atom
    param: x2,y2,z2= the coordinates of the second atom
    returns: distance between atoms and details of the calculation
    """

    xd = float(x2) - float(x1)
    yd = float(y2) - float(y1)
    zd = float(z2) - float(z1)
    sumsqu = (xd * xd) + (yd * yd) + (zd * zd)
    mag = math.sqrt(sumsqu)
    #if mag > 10:
    #    print("mag",mag)
    return (mag)

def getAngle(x1, y1, z1, x2, y2, z2, x3, y3, z3):
    """
    param: x1,y1,z1= the coordinates of the first atom
    param: x2,y2,z2= the coordinates of the second atom
    param: x3,y3,z3= the coordinates of the third atom
    returns: angle between atoms and details of the calculation
    """

    xd = float(x2) - float(x1)
    yd = float(y2) - float(y1)
    zd = float(z2) - float(z1)
    xe = float(x2) - float(x3)
    ye = float(y2) - float(y3)
    ze = float(z2) - float(z3)
    dot = (xd * xe) + (yd * ye) + (zd * ze)
    magA = math.sqrt((xd * xd) + (yd * yd) + (zd * zd))
    magB = math.sqrt((xe * xe) + (ye * ye) + (ze * ze))
    cos_theta = dot / (magA * magB)
    theta = math.acos(cos_theta)
    theta_deg = (theta / 3.141592653589793238463) * 180
    return (round(theta_deg, 3))

def crossProduct(A, B):
    x = (A[1] * B[2]) - (A[2] * B[1])
    y = (A[2] * B[0]) - (A[0] * B[2])
    z = (A[0] * B[1]) - (A[1] * B[0])
    return (x, y, z)

def getDihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4):
    """
    param: x1,y1,z1= the coordinates of the first atom
    param: x2,y2,z2= the coordinates of the second atom
    param: x3,y3,z3= the coordinates of the third atom
    param: x4,y4,z4= the coordinates of the fourth atom
    returns: angle between atoms and details of the calculation
    """

    xA = float(x2) - float(x1)
    yA = float(y2) - float(y1)
    zA = float(z2) - float(z1)

    xB = float(x2) - float(x3)
    yB = float(y2) - float(y3)
    zB = float(z2) - float(z3)

    xC = float(x4) - float(x3)
    yC = float(y4) - float(y3)
    zC = float(z4) - float(z3)

    crossAB = crossProduct((xA, yA, zA), (xB, yB, zB))
    crossBC = crossProduct((xB, yB, zB), (xC, yC, zC))
    dot = (crossAB[0] * crossBC[0]) + (crossAB[1] * crossBC[1]) + (crossAB[2] * crossBC[2])

    ABsq = (crossAB[0] ** 2) + (crossAB[1] ** 2) + (crossAB[2] ** 2)
    BCsq = (crossBC[0] ** 2) + (crossBC[1] ** 2) + (crossBC[2] ** 2)
    magAB = math.sqrt(ABsq)
    magBC = math.sqrt(BCsq)

    try:
        cos_theta = dot / (magAB * magBC)
        theta = math.acos(cos_theta)
        theta_deg = (theta / 3.141592653589793238463) * 180
        cross = crossProduct(crossAB, crossBC)
        dotB = (cross[0] * xB) + (cross[1] * yB) + (cross[2] * zB)
        if dotB > 0:
            theta_deg *= -1

        return (round(theta_deg, 3))
    except:
        return -180 #there is a 0 max for cos theta