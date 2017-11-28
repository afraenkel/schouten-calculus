
import sympy
from multivector import Mv


class stdSympl(object):

    def __init__(self, N):
        self.dim = N
        self.x = sympy.symbols(['x%s' % x for x in range(N)])
        self.y = sympy.symbols(['y%s' % x for x in range(N)])
        self.Mv = Mv(self.x + self.y, {(s, t): 1 for (s, t) in zip(self.x, self.y)})
        
    def brac(self, f, g):
        return self.Mv.brac(f, g)


class sl2(object):

    def __init__(self):
        x, y, z = sympy.symbols('x y z')
        self.Mv = Mv([x, y, z], {(x, y): z, (y, z): x, (z, x): y})

    def brac(self, f, g):
        return self.Mv.brac(f, g)
