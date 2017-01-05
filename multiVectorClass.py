#!/usr/bin/env python

import sympy
from collections import defaultdict
from copy import deepcopy
from functools import reduce
from operator import mul


class Mv(object):

    def __init__(self, gens, P=defaultdict(int)):
        self.gens = gens
        self._P = P

    def __str__(self):
        return str(self._to_poly())

    def __repr__(self):
        return repr(self._to_poly())

    def __add__(self, other):
        P = deepcopy(self._P)
        for v, c in other._P.items():
            P[v] += c
        return Mv(self.gens, P)

    def __sub__(self, other):
        return (self + (-1) * other)

    def __rmul__(self, symb):
        P = deepcopy(self._P)
        for v in P:
            P[v] *= symb
        return Mv(self.gens, P)
    
    def __mul__(self, other):
        P = defaultdict(int)
        for v1, c1 in other._P.items():
            for v2, c2 in self._P.items():
                P[v1 + v2] += (c1 * c2)
        return Mv(self.gens, P)

    def brac(self, *F):
        s = list(self.sdeg().keys())
        if (len(s) != 1):
            raise ValueError("brackets only allowable from homogeneous tensors")
        if (s[0] != len(F)):
            raise ValueError("not enough input functions")
        
        mvs = [Mv(self.gens, {(): f}) for f in F]
        A = sBr(self, mvs.pop())
        while mvs:
            A = sBr(A, mvs.pop())

        return A._P[()]

    def deg(self):
        '''decompose the multivector w/r/t the polynomial degree of the coefficients'''
        degdict = defaultdict(lambda: defaultdict(int))
        for v, c in self._P.items():
            try:
                monoms = zip(map(sum, c.as_poly().monoms()), c.as_coeff_add()[1])
                for d, coeff in monoms:
                    degdict[d][v] += coeff
            except:
                degdict[0][v] += c

        return {n: Mv(self.gens, P) for (n, P) in degdict.items()}

    def diff(self, x):
        '''differentiate the coefficients w/r/t x'''
        P = defaultdict(int)
        for v, c in self._P.items():
            try:
                P[v] = c.diff(x)
            except AttributeError:
                pass
        return Mv(self.gens, P)

    def i(self, other):
        '''contraction of multivector :self: by :other:'''
        raise NotImplementedError
    
    def sBr(self, other):
        return sBr(self, other)
    
    def sdeg(self):
        '''decompose the multivector w/r/t the multivector degree'''
        degdict = defaultdict(lambda: defaultdict(int))
        for v, c in self._P.items():
            degdict[len(v)][v] += c

        return {n: Mv(self.gens, P) for (n, P) in degdict.items()}

    def sdiff(self, x):
        '''differentiate the vector basis w/r/t x'''
        P = defaultdict(int)
        for v, c in self._P.items():
            try:
                idx = v.index(x)
                w = v[:idx] + v[idx + 1:]
                P[w] = (-1)**(len(self.gens) - idx) * c
            except ValueError:
                pass
        return Mv(self.gens, P)
    
    def sort(self):
        return Mv(self.gens, sort(self.gens, self._P))
    
    def _to_poly(self):
        s = 0
        for (v, c) in self._P.items():
            vec = map(lambda x: sympy.symbols('d%s' % x.__str__(), commutative=False), v)
            s += c * reduce(mul, vec, 1)
        return s


def sBr(A, B):
    '''take the schouten bracket of A and B'''
    gens = tuple(set(A.gens) | set(B.gens))
    out = Mv(gens)

    for i, Ai in A.sdeg().items():
        for j, Bj in B.sdeg().items():
            for var in gens:
                sgn = (-1)**((i - 1) * (j - 1) + 1)
                out += (Ai.sdiff(var) * Bj.diff(var) + sgn * Bj.sdiff(var) * Ai.diff(var))

    return out


def sort(gens, P):
    out = defaultdict(int)
    for vals, coeff in P.items():
        s, parity = _sort(vals, gens)
        out[s] += parity * coeff
    return out


def _sort(values, gens):
    '''
    (c0,c3,c2), (c0, c1, c2, c3) => ((c0,c2,c3),-1)
    '''
    idxs = [gens.index(x) for x in values]

    N = len(idxs)
    num_swaps = 0
    for i in range(N - 1):
        for j in range(i + 1, N):
            if idxs[i] == idxs[j]:
                return values, 0
            if idxs[i] > idxs[j]:
                idxs[i], idxs[j] = idxs[j], idxs[i]
                num_swaps += 1
    return tuple(gens[i] for i in idxs), (-1)**(num_swaps % 2)
        

# ---------------------------------------------------------------------
# misc
# ---------------------------------------------------------------------

class stdSympl(object):

    def __init__(self, N):
        self.dim = N
        self.x = sympy.symbols(['x%s' % x for x in range(N)])
        self.y = sympy.symbols(['y%s' % x for x in range(N)])
        self.Mv = Mv(self.x + self.y, {(s, t): 1 for (s, t) in zip(self.x, self.y)})
        
    def brac(self, f, g):
        return self.Mv.brac(f, g)


# ---------------------------------------------------------------------
# lifting coords
# ---------------------------------------------------------------------

    
def lift(expr, gb, invHom, limit=10):
    '''
    returns a lifted expression, where expr is
    rewritten using the symbols in embRing.
    '''
    if not limit:
        print("limit reached before completion")
        return expr
    
    p = gb.reduce(expr)
    inverse = {x: y for (y, x) in invHom.items()}
    embRing = [inverse[x.as_expr()] for x in gb.args[0]]
    
    subbed = sum(x * y for (x, y) in zip(p[0], embRing)) + p[1]

    if p[1]:
        raise
    orig_gens_left = subbed.atoms() & (set(gb.gens) - set(embRing))

    if not len(orig_gens_left):
        return subbed
    else:
        inv = list(map(lambda x: x.as_expr(), gb.args[0]))
        gens = set(gb.gens) | set(embRing)
        gb = sympy.polys.groebner(inv, *gens)
        return lift(subbed, gb, invHom, limit=limit - 1)

    
def lift_from_quotient(invHom):

    gb = sympy.polys.groebner(invHom.values())
    embRing = invHom.keys()

    brac = stdSympl(len(embRing)).brac
    
    P = {}
    for k1, i1 in enumerate(embRing):
        for k2, i2 in enumerate(embRing):
            if k1 < k2:
                p = brac(invHom[i1], invHom[i2])
                P[(i1, i2)] = lift(p.as_expr(), gb, invHom)
    return P

'''
def poissBr(P):
    Pfull = deepcopy(P)
    Pfull.update({(k2, k1): -1 * v for ((k1, k2), v) in P.items()})

    def bracket(f, g):
        return sum(
            f.diff(xi) * g.diff(yj) * Pfull[(xi, yj)] for (xi, yj) in Pfull.keys()
        )

    return bracket
'''
