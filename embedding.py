
import sympy
from examples import stdSympl


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
