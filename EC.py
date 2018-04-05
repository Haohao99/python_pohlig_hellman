######################################################################
# Matt Watson
# Directed Reading Program Project
# Elliptic Curve Implementation
######################################################################

# Class for representing an elliptic curve
class Elliptic_Curve:
    def __init__(self, b, c):
        self.b = b
        self.c = c
    
    # Returns boolean value for if the point is on the curve
    def is_point(self, pt):
        if pt.is_inf:
            return True
        return (pow(pt.y, 2) == pow(pt.x, 3) + self.b * pt.x + self.c)
    
    # adds two points on the curve together
    def add(self, pt1, pt2):
        # Check that points are on curve
        if not self.is_point(pt1) or not self.is_point(pt2):
            print('Given points are invalid')
            return Point(None, None, is_infinty=True)
        
        # Consider infinity cases 
        if pt1.is_inf and pt2.is_inf:
            return Point(None, None, is_infinity=True)
        if pt1.is_inf:
            return Point(pt2.x, pt2.y)
        if pt2.is_inf:
            return Point(pt1.x, pt1.y)
        
        # Add points otherwise
        if pt1.x == pt2.x and pt1.y == pt2.y:
            if pt1.y == 0:
                return Point(None, None, is_infinity=True)
            m = (3 * pow(pt1.x, 2) + self.b) / (2 * pt1.y)
        else:
            m = (pt2.y - pt1.y) / (pt2.x - pt1.x)
        x_3 = pow(m, 2) - pt1.x - pt2.x
        y_3 = (m * (pt1.x - x_3)) - pt1.y
        return Point(x_3, y_3) 
        
# Class for a point. Supports infinity
class Point:
    def __init__(self, x, y, is_infinity=False):
        self.x = x
        self.y = y
        self.is_inf = is_infinity
        if is_infinity:
            self.coordinate = None
    
    def __str__(self):
        if self.is_inf:
            return 'inf'
        else:
            return '(' + str(self.x) + ', ' + str(self.y) + ')'

if __name__ == '__main__':
    ec = Elliptic_Curve(0, 73)
    pt1 = Point(2, 9)
    print(ec.is_point(pt1))
    pt2 = Point(3, 10)
    pt3 = ec.add(pt1, pt2)
    print(pt3)
    print(ec.add(pt3, pt3))
    inf = Point(None, None, is_infinity=True)
    print(ec.add(pt1, inf))

