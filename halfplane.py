from sage.all import *

def inequality_to_geodesic(a, b, c):
    g = HyperbolicPlane().UHP().get_geodesic

    if a == 0:
        if b == 0:
            raise ValueError("invalid inequality with a == b == 0")
        else:
            if b < 0:
                return g(-c / b, oo)
            else:
                return g(oo, -c / b)
    else:
        discriminant = b**2 - 4 * a * c
        if bool(discriminant <= 0):
            raise ValueError("discriminant was non-positive")
        else:
            center = -b / (2 * a)
            # TODO: taking sqrt this way forces floating point. Is that bad?
            left_root = center - (sqrt(discriminant) / (2 * a))
            right_root = center + (sqrt(discriminant)  / (2 * a))
            if a * center**2 + b * center + c > 0:
                return g(right_root, left_root)
            else:
                return g(left_root, right_root)

def is_trivial(ineq):
    return bool(ineq[1]**2 - 4*ineq[0]*ineq[2] <= 0)

def plot(geodesic):
    start = geodesic.start()
    end = geodesic.end()
    if not (start == oo or end == oo):
        if start < end:
            return geodesic.plot(axes=True, color="green")
        else:
            return geodesic.plot(axes=True, color="red")
    elif start == oo:
        return geodesic.plot(axes=True, color="red")
    else:
        return geodesic.plot(axes=True, color="green")


