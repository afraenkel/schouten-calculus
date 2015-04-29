"""
Implements a class for multivector fields over the ring of sympy symbols.
Includes wedge product, schouten bracket, etc...
"""

from sympy import Symbol, diff


#---------------------------------------------------------------------
# Globals / setting coordinates
#---------------------------------------------------------------------

COORDS = list()

def setCoordinates(dim,symb='x'):
    """Sets the dimension and variable symbol (right now, only 'x' works). """
    global COORDS
    COORDS = [ Symbol('%s_%d' %(symb,k)) for k in range(dim) ]
    pass


#---------------------------------------------------------------------
# MultiVector Class definition
#---------------------------------------------------------------------

class MultiVect(object):
    """
    Class of multivectors over a ring (any class with well defined +/*).
    a MultiVect object is a list of monomial multivectors, accessed via self.monomials.
    """
    
    def __init__(self, X=[] ):
        # TODO: instatiate with single elt / single tuple...
        mons = list()
        for x in X:
            if type(x) is not _MultiVectMonomial:
                mons.append( _MultiVectMonomial(x) )
            else:
                mons.append( x )
        if len(mons) == 0:
            mons = [ _MultiVectMonomial(0) ]
        self.monomials = mons

    def __str__(self):
        s = ''
        # sort the monomials by degree and print them separated by ' + '.
        for monom in sorted(self.order().monomials, key=lambda x:x.degree() ):
            if str(monom):
                monomStr = "%s + " %str( monom )
                s += monomStr
        if s:
            # strip the last ' + ' sign.
            return s[:-3]
        else:
            return '0'

    def __add__(self,other):
        return MultiVect( self.monomials + other.monomials ).order()
    
    def __sub__(self,other):
        return (self + (-1)*other).order()

    def __mul__(self,other):
        mulVect = list()
        for mon1 in self.monomials:
            for mon2 in other.monomials:
                mulVect.append( mon1*mon2 )
        return MultiVect(mulVect).order()
    
    def __rmul__(self,other):
        """Multiplication of a multivector and a functions."""
        rmulVect = list()
        for mon in self.monomials:
            rmulVect.append(_MultiVectMonomial(other)*mon)
        return MultiVect(rmulVect).order()

    def sBr(self,other):
        """ Schouten Bracket of self and other. """
        return sBr(self,other)
    
    def degree(self):
        """ decomposes a multivector into a list of unique degrees """
        return sorted(set([ mon.degree() for mon in self.monomials ]))

    def order(self):
        """
        Orders the anticommuting multivectors into ascending order 
        and collects the summands. 
        """
        ordered = dict()
        for mon in self.monomials:
            sgn,sortedMonomial = _signedSort( mon.mvect )
            # if the term isn't 0 either add it to sum, 
            # or combine with a previous term.
            if sgn != 0:
                if sortedMonomial in ordered.keys():
                    ordered[ sortedMonomial ] += sgn*mon.coeff
                else:
                    ordered[ sortedMonomial ] = sgn*mon.coeff
        # Delete leading zeros:
        try:
            if ordered[()] == 0:
                del ordered[()]
        except KeyError:
            pass
        return MultiVect(ordered.items())

#---------------------------------------------------------------------
# multivector monomial class
#---------------------------------------------------------------------

class _MultiVectMonomial(object):
    """
    A class implementing monomial multivectors. Purely a helper class for 
    defining the class MultiVector.  Monomials are represented as 
    tuples (mvect,coeff)  where mvect is a tuple representing the anticommuting
    part of the multvector.  e.g. x1*x3*dx1*dx4*d2 <----> ( (1,4,2), x1*x3 )
    """
    
    def __init__(self,d=((),0) ):
        try:
            self.coeff = d[1]
            self.mvect = d[0]
        except:
            self.coeff = d
            self.mvect = ()

    def __str__(self):
        if self.mvect == ():
            return str( self.coeff )
        else:
            if self.coeff == 0:
                return ''
            else:
                d = '*'.join([ "dx%d" %i for i in self.mvect ])
                return "(%s)*%s" %( str( self.coeff), d )

    def __mul__(self,other):
        """multiplication of multivectors. """
        mulVect =  _MultiVectMonomial()
        mulVect.mvect = self.mvect + other.mvect
        mulVect.coeff = self.coeff * other.coeff
        return mulVect.order()

    def __rmul__(self,other):
        """scalar multiplication """
        scalVect =  _MultiVectMonomial()
        scalVect.mvect = self.mvect
        scalVect.coeff = self.coeff * other
        return scalVect
                                
    def degree(self):
        """ decomposes a multivector into a dict of degrees """
        if self == ( (), 0 ):
            return 'minus infinity'
        else:
            return len( self.mvect)

    def diff(self,var,deg=0):
        deriv = _MultiVectMonomial()
        if deg == 0:
            deriv.mvect = self.mvect
            deriv.coeff = self.coeff.diff(var)
        else:
            try:
                ind = self.mvect.index(var)
            except ValueError:
                return _MultiVectMonomial(0)
            deriv.mvect = self.mvect[:ind] + self.mvect[ind+1:]
            deriv.coeff = (-1)**(len(self.mvect)-ind)*self.coeff
        return deriv

    def order(self):
        sgn,sortedMonomial = _signedSort( self.mvect )
        if sgn != 0 and self.coeff != 0 :
            ordered = _MultiVectMonomial()
            ordered.mvect = sortedMonomial
            ordered.coeff = sgn*self.coeff
            return ordered
        else:
            return _MultiVectMonomial( 0 )

#---------------------------------------------------------------------
# Aux. Functions for MultiVector / MultiVectorMonomial classes
#---------------------------------------------------------------------

def _signedSort(L):
    """
    Returns a tuple (sgn,T) where T is a sorted version of L and
    sgn is {-1,0,1}, the sign of the permutation taking L to T.
    """

    if len(L)<2:
        return (1,tuple(L))
    
    L1,S = list(L[:]),sorted(L)
    for k in range(len(S)-1):
        if S[k] == S[k+1]:
            return (0,tuple(S))
        
    sgn = 1
    while L1:
        k = L1.pop(0)
        ind = S.index(k)
        if ind % 2:
            sgn *= -1
        del S[ind]
    return (sgn,tuple(sorted(L)))

#---------------------------------------------------------------------
# Schouten Bracket Functions
#---------------------------------------------------------------------

def sBr(X,Y):
    """Schouten Bracket [X,Y]. """
    V = MultiVect([])
    for mon1 in X.monomials:
        for mon2 in Y.monomials:
            V = V + _sBrMonomials(mon1,mon2)
    return V.order()

def _sBrMonomials(mon1,mon2):
    """ Schouten Bracket between two monomial multivectors. """
    outVect = MultiVect([])
    for k,var in enumerate(COORDS):
        term1 = mon1.diff(k,deg=1)*mon2.diff(var)
        term2 = mon2.diff(k,deg=1)*mon1.diff(var)
        sgn = (-1)**( (mon1.degree() - 1)*(mon2.degree() - 1) + 1)
        outVect +=  MultiVect( [ term1, sgn*term2 ] ).order()
    return outVect

#---------------------------------------------------------------------
# Other functions
#---------------------------------------------------------------------

# Poisson Bracket

# Hamiltonian Vector Field

# Poisson Algebra Class (?)

# differential for Poisson (Co)-Homology  (may be useful -- use sympy
# solve resulting equations.)

# Contraction / inner product

# Curl Operator! (easy!)

# canonical examples (Volume Form, std Poisson str, sl2, etc...)
