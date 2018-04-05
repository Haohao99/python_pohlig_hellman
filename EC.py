######################################################################
# Matt Watson
# Directed Reading Program Project
# Elliptic Curve Implementation
######################################################################

import DRP

# Class for representing an elliptic curve
class Elliptic_Curve:
    def __init__(self, b, c, p=None):
        self.b = b
        self.c = c
        self.p = p
    
    # Returns boolean value for if the point is on the curve
    def is_point(self, pt):
        if pt.is_inf:
            return True
        if self.p is None:
            return (pow(pt.y, 2) == pow(pt.x, 3) + self.b * pt.x + self.c)
        else:
            return (pow(pt.y, 2, self.p) == 
                    (pow(pt.x, 3) + self.b * pt.x + self.c) % self.p)
    
    # adds two points on the curve together
    def add(self, pt1, pt2):
        # Check that points are on curve
        if not self.is_point(pt1) or not self.is_point(pt2):
            raise ValueError('The points being added are invalid.')
        
        # Consider infinity cases 
        if pt1.is_inf and pt2.is_inf:
            return Point(None, None, is_infinity=True)
        if pt1.is_inf:
            return Point(pt2.x, pt2.y)
        if pt2.is_inf:
            return Point(pt1.x, pt1.y)
        
        # Add points otherwise
        if self.p is None:
            if pt1.x == pt2.x and pt1.y == pt2.y:
                if pt1.y == 0:
                    return Point(None, None, is_infinity=True)
                m = (3 * pow(pt1.x, 2) + self.b) / (2 * pt1.y)
            else:
                m = (pt2.y - pt1.y) / (pt2.x - pt1.x)
            x_3 = pow(m, 2) - pt1.x - pt2.x
            y_3 = (m * (pt1.x - x_3)) - pt1.y
            return Point(x_3, y_3)
        else:   # Mod p
            if pt1.x == pt2.x and pt1.y == pt2.y:
                if pt1.y == 0:
                    return Point(None, None, is_infinity=True)
                m = (3 * pow(pt1.x, 2) + self.b) 
                m = (m * DRP.inverse((2 * pt1.y), self.p)) % self.p
            else:
                m = (pt2.y - pt1.y)
                m = (m * DRP.inverse((pt2.x - pt1.x), self.p)) % self.p
            x_3 = (pow(m, 2) - pt1.x - pt2.x) % self.p
            y_3 = ((m * (pt1.x - x_3)) - pt1.y) % self.p
            return Point(x_3, y_3)

    # Pohlig-Hellman analog
    def pohlig_hellman(self, A, B):
        if not self.is_point(A) or not self.is_point(B):
            raise ValueError('The points are not valid.')
        # First, find smallest integer n so that nA = inf
        n = 1
        while True:
            A_aux = A
            if A_aux.is_inf:
                break
            A_aux = add(A, A_aux)
            n = n + 1
        ######## This is where I made it to
    
    # Check a k value by returning point B s.t. B = kA
    def check_discrete_log_ans(self, k, A):
        original_A, ans, i = Point(A.x, A.y, A.is_inf), A, 0
        while i < k:
            ans = self.add(ans, A)
        return ans
    
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
    # Working on regular elliptic curve
    # Example p. 350
    ec = Elliptic_Curve(0, 73)
    pt1 = Point(2, 9)
    print(ec.is_point(pt1))
    pt2 = Point(3, 10)
    pt3 = ec.add(pt1, pt2)
    print(pt3)
    print(ec.add(pt3, pt3))
    inf = Point(None, None, True)
    print(ec.add(pt1, inf))
    
    # Working on elliptic curve mod 5
    # Found on p. 352
    ec_mod = Elliptic_Curve(4, 4, 5)
    pt1 = Point(1, 2)
    pt2 = Point(4, 3)
    pt3 = Point(3, 2)
    print(ec_mod.is_point(pt1)) # Should be on
    print(ec_mod.is_point(pt2)) # Should be on
    print(ec_mod.is_point(pt3)) # Should not be on
    print(ec_mod.is_point(inf)) # Should be on
    print(ec_mod.add(pt1, pt2))

    # Example p. 353
    ec_mod = Elliptic_Curve(4, 4, 2773)
    pt = Point(1, 3)
    print(ec_mod.add(pt, pt))

    




    
    

