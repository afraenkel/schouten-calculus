# schouten-calculus
An implementation of the Schouten algebra of multi-vector fields over the ring of sympy symbols.

This isn't supposed to be blazing fast (it is sympy-based).
Motivation: easy interaction / ability to derive algrebraic equation that can be solved using other means (e.g. using Macaulay2).

## Currently uses a naive 'list of monomials' implementation to lay down the structure.
* TODO: Sparse dict-based version would be a fairly easy tweak.
* TODO: A numpy array-based version would also be a fairly easy tweak.

## Features to be added:

* Poisson Bracket
* Hamiltonian Vector Field
* Poisson Algebra Class (?)
* Differential for Poisson (Co)-Homology  (may be useful -- use sympy / M2 to solve resulting equations.)
* Contraction / inner product
* Curl Operator! (easy!)
* Canonical examples (Volume Form, std Poisson str, sl2, etc...)
