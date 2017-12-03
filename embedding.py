
import sympy


def lift(rmap, domain, codomain):
    '''
    returns a function lifting elements of the codomain of
    a ring homomorphism rmap to the domain.
    '''

    inv_rmap = {y: x for (x, y) in rmap.items()}
    gb = sympy.polys.groebner(
        rmap.values(),
        *(tuple(domain) + tuple(codomain))
    )

    coordlist = [inv_rmap[x] for x in gb.exprs]
    
    def closure(expr, limit=10):
        if not limit:
            print('limit reached before completion')
            return expr

        b, r = gb.reduce(expr)

        if r != 0:
            print("groebner reduction has non-zero remainder")
            print(b, r)
            return r

        subbed = sum(x * y for (x, y) in zip(b, coordlist))
        try:
            orig_gens_left = subbed.atoms() & set(codomain)
        except:
            orig_gens_left = set()

        if not orig_gens_left:
            return subbed
        else:
            return closure(subbed, limit=(limit - 1))

    return closure

