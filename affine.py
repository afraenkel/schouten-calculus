
import numpy as np


class Polynomial(object):
    
    def __init__(self, monomials):
        self.terms = monomials
        
    def __repr__(self):
        return self.terms.__repr__()
    
    def _repr_latex_(self):
        return _repr_latex_(self.terms)

    def __add__(self, other):
        return Polynomial(np.vstack([self.terms, other.terms]))

    def __sub__(self, other):
        return self + -1 * other

    def __pow__(self, N):
        if N == 0:
            return Polynomial(np.array([[1, 0]]))
        elif N == 1:
            return self
        else:
            return self * (self ** (N - 1))

    def __rmul__(self, other):
        new = self.terms.copy()
        new[:, 0] = new[:, 0] * other
        return Polynomial(new)
    
    def __mul__(self, other):
        
        p, q = self.terms, other.terms
        coeffs = np.outer(p[:, 0], q[:, 0]).flatten()
        a, b = q[:, 1:], p[:, 1:]
        r, s = np.meshgrid(a, b)
        N, M = a.shape[1], min(a.shape[0], b.shape[0])

        idxs = np.indices((r+s).shape).T
        idxs = np.diff(idxs, axis=-1) % N == 0
        cc = (r+s)[idxs.T.squeeze()].reshape(-1, M).T.reshape(-1, N)
        cc = np.hstack(np.split(cc, M)).reshape(-1, N)
        
        out = np.hstack([coeffs.reshape(-1, 1), cc])
        return Polynomial(out)


def _repr_latex_(arr):
    """
    latex representation of a multivariate polynomial
    """
    coeffs, monomials = arr[:, 0], arr[:, 1:]
    toprint = []
    for monomial, coeff in zip(monomials, coeffs):
        terms = []
        for k, x in enumerate(monomial):
            if x == 1:
                terms.append('x_{%d}' % k)
            elif x > 1:
                terms.append('x_{%d}^{%d}' % (k, x))

        if coeff == 0:
            toprint.append('0')
        elif coeff == 1:
            toprint.append(''.join(terms))
        else:
            toprint.append(str(coeff) + ''.join(terms))

    return '$' + '+'.join(toprint) + '$'
