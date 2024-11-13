from sympy import symbols, div, mod_inverse, gcd, Poly, primefactors, FiniteField, GF
from sympy.polys.polytools import invert
from sympy.abc import x
import itertools
import math

class FiniteFieldElement:
    def __init__(self, coeffs, p, k, irreducible_poly):
        self.poly = Poly(coeffs, x, domain=GF(p))
        self.p = p
        self.k = k
        self.irreducible_poly = Poly(irreducible_poly, x, domain=GF(p))

    def __eq__(self,other):
        return (self.poly==other.poly) and (self.p==other.p) and (self.k==other.k) and (self.irreducible_poly==other.irreducible_poly)

    def __neg__(self):
        result_coeffs = (-self.poly) % self.irreducible_poly
        return FiniteFieldElement(result_coeffs, self.p, self.k, self.irreducible_poly)
        
    def __add__(self, other):
        if isinstance(other, int):
            return self.__radd__(other)
        
        result_coeffs = (self.poly + other.poly) % self.irreducible_poly
        return FiniteFieldElement(result_coeffs, self.p, self.k, self.irreducible_poly)
    
    def __sub__(self, other):
        if isinstance(other, int):
            return -self.__rsub__(other)
        
        result_coeffs = (self.poly - other.poly) % self.irreducible_poly
        return FiniteFieldElement(result_coeffs, self.p, self.k, self.irreducible_poly)
    
    def __mul__(self, other):
        if isinstance(other, int):
            return self.__rmul__(other)
        
        result_coeffs = (self.poly * other.poly) % self.irreducible_poly
        return FiniteFieldElement(result_coeffs, self.p, self.k, self.irreducible_poly)

    def __radd__(self, other):
        result_coeffs = (self.poly + other)
        return FiniteFieldElement(result_coeffs, self.p, self.k, self.irreducible_poly)

    def __rsub__(self, other):
        result_coeffs = (other - self.poly)
        return FiniteFieldElement(result_coeffs, self.p, self.k, self.irreducible_poly)

    def __rmul__(self, other):
        result_coeffs = (self.poly * other)
        return FiniteFieldElement(result_coeffs, self.p, self.k, self.irreducible_poly)
    
    def __truediv__(self, other):
        inv_poly = invert(other.poly, self.irreducible_poly)
        result_coeffs = (self.poly * inv_poly) % self.irreducible_poly
        return FiniteFieldElement(result_coeffs, self.p, self.k, self.irreducible_poly)

    def __pow__(self, other):
        assert other>=0
        if other==0:
             return FiniteFieldElement(1, self.p, self.k, self.irreducible_poly)
        if other%2==0:
            return (self*self)**(other//2)
        else:
            return (self*self)**(other//2)*self
    def frob(self):
        return self**self.p
    def trace(self):
        result = FiniteFieldElement(0, self.p, self.k, self.irreducible_poly)
        result+=self
        for i in range(self.k-1):
            result = result.frob()
            result+=self
        return result.poly()
        
        
    def __str__(self):
        return str(self.poly.as_expr())

# finds an irreducible polynomial of degree k in F_p[x], give this as input to FiniteFieldElement
def find_irreducible_poly(p, k):
    x = symbols('x')
    if k==1:
        return Poly([1,0],x,domain=GF(p))
    for coeffs in itertools.product(range(p), repeat=k):
        if coeffs[-1] == 0:  # Skip if the leading coefficient is zero
            continue
        poly = Poly((1,)+coeffs, x, domain=GF(p))
        if is_irreducible(poly, p, k):
            return poly
    return None

def is_irreducible(poly, p, k):
    x = symbols('x')
    q = p**k
    # Check if poly divides x^q - x in F_p

    a=FiniteFieldElement([1,0],p,k,poly)
    if a**q!=a:
        return False
    for d in primefactors(k):
        if k==1:
           continue
        if gcd(poly.as_expr(), (a**(p**(k//d)) - a).poly, domain=GF(p)) != 1:
            return False
    return True
