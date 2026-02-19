from __future__ import annotations
import numpy as np
from math import factorial
from collections import Counter


def isNumeric(num) -> bool:
    return (isinstance(num, float) or isinstance(num, int)
            or isinstance(num, complex))


def isReal(num) -> bool: return isinstance(num, float) or isinstance(num, int)


def isInt(thing) -> bool:
    return (type(thing) is float and thing.is_integer()) or type(thing) is int


def isPoly(other) -> bool:
    """ Checks whether the input is a term (deg 0 poly) or a polynomial"""
    return isinstance(other, Term) or isinstance(other, Polynomial)


def roundComplex(x, a: int) -> bool:
    """Rounds complex number to specified digits"""
    if isinstance(x, complex):
        return round(x.real, a) + round(x.imag, a)*1j
    return round(x, a)


def arrayToTerm(arr: list[float]) -> list[Term]:
    deg = len(arr) - 1
    output = [0] * (deg + 1)
    for i in range(deg + 1):
        output[i] = Term(arr[i], deg-i)
    return output


class Function():
    """
    Abstract Function Class. This should not be instantiated.
    """

    def __repr__(self) -> str:
        return str(self)

    def __rmul__(self, other: Function | int | float | complex) -> Function:
        return self * other

    def __sub__(self, other: Function | int | float | complex) -> Function:
        return self + -1*other

    def __rsub__(self, other: Function | int | float | complex) -> Function:
        return -1*self + other

    def __radd__(self, other: Function | int | float | complex) -> Function:
        return self + other

    def evaluate(self, a: int | float | complex) -> int | float | complex:
        pass

    def getPoles(self) -> list[np.float64 | np.complex64]:
        return []

    def limit(self, a: complex) -> complex:
        return self.evaluate(a)


class Term(Function):
    """
    A mathematical object defined by ax^n.
    """

    def __init__(self, coefficient: int | float | complex, power: int) -> None:
        """
        Constructs a Term of a polynomial with a numeric coefficient and integer power.
        """
        if not isInt(power) or power < 0:
            raise Exception("Invalid Power")
        if not isNumeric(coefficient):
            raise Exception("Coefficient is not a number")
        self._coefficient = coefficient
        self._power = power

    def getPower(self) -> int:
        return self._power

    def getCoefficient(self) -> int | float | complex:
        return self._coefficient

    def __str__(self) -> str:
        coeff = self.getCoefficient()
        sign = "" if ((type(coeff) is complex and coeff.real > 0)
                      or (isReal(coeff) and coeff > 0)) else "-"
        if coeff == 0:
            return "0"
        elif self.getPower() == 0:
            return f"{self._coefficient}"
        elif abs(coeff) != 1 and self.getPower() != 1:
            return f"{sign}{abs(coeff)}x^{self._power}"
        elif abs(coeff) == 1 and self.getPower() != 1:
            return f"{sign}x^{self._power}"
        elif abs(coeff) != 1 and self.getPower() == 1:
            return f"{sign}{abs(coeff)}x"
        else:
            return f"{sign}x"

    def __lt__(self, other: Term) -> bool:
        return self.getPower() < other.getPower()

    def __eq__(self, other) -> bool:
        if self.isZero() and isNumeric(other):
            return other == 0
        if isNumeric(other) and self.getPower() == 0:
            return other == self.getCoefficient()
        if isinstance(other, Term):
            return self.getPower() == other.getPower()
        return False

    def __gt__(self, other: Term) -> bool:
        return self.getPower() > other.getPower()

    def __add__(self, other) -> Polynomial:
        p1 = Polynomial(self.toArray())
        if isNumeric(other):
            p2 = Polynomial([other])
        else:
            p2 = Polynomial(other.toArray())
        return p1 + p2

    def __mul__(self, other) -> Term:
        if isinstance(other, Term):
            return Term(self.getCoefficient() * other.getCoefficient(),
                        self.getPower() + other.getPower())
        elif isNumeric(other):
            return Term(other * self.getCoefficient(), self.getPower())

    def derivative(self) -> Term:
        if self.getPower() > 0:
            return Term(self.getPower() * self.getCoefficient(), self.getPower() - 1)
        if self.getPower() == 0:
            return Term(0, 0)
        else:
            raise Exception("no")

    def __pow__(self, pow: int) -> Term:
        return Term(self.getCoefficient()**pow, self.getPower() * pow)

    def evaluate(self, a: int | float | complex) -> int | float | complex:
        return self.getCoefficient() * (a)**(self.getPower())

    def compose(self, p) -> int | float | complex | Function:
        if isNumeric(p):
            return self.evaluate(p)
        if isinstance(p, Term) and p.getPower() == 1 and p.getCoefficient() == 1:
            return self
        elif isPoly(p):
            return self.getCoefficient() * p ** self.getPower()
        elif isinstance(p, Rational):
            return Rational(self.getCoefficient() * p.getNum()**self.getPower(),
                            p.getDenom()**self.getPower())

    def toArray(self) -> list[complex]:
        return [self.getCoefficient()] + [0] * self.getPower()

    def getRoots(self) -> list[int]:
        return [0]

    def real(self) -> Term:
        c = self.getCoefficient()
        return Term(c.real, self.getPower()) if isinstance(c, complex) else self

    def imag(self) -> Term:
        c = self.getCoefficient()
        return Term(c.imag, self.getPower()) if isinstance(c, complex) else 0

    def isZero(self) -> bool:
        return self.getCoefficient() == 0

    def __neg__(self) -> bool:
        return Term(-self.getCoefficient(), self.getPower())


class Polynomial(Function):
    """
    A polynomial in the form of ax^n + bx^n-1 + ... + nx + o
    """

    def __init__(self,
                 vals: list[Term | int | float | complex] = [Term(0, 0)]) -> None:
        """
        Constructs a polynomial via either a list of Term objects
            or descending coefficients ie [a, b,..., n, o] 
        """
        deg = len(vals) - 1
        for index, val in enumerate(vals):
            if isinstance(val, Term):
                continue
            elif isNumeric(val):
                vals[index] = Term(val, deg - index)

        self._terms: list[Term] = sorted(vals, reverse=True)

    def getPower(self) -> int:
        return len(self._terms) - 1

    def __add__(self, term: float | int | Term | Polynomial) -> Polynomial:
        if isNumeric(term):
            arr = self.toArray()
            t = arr[-1]
            arr[-1] = t + term
            return Polynomial(arr)

        elif isPoly(term):
            new = []
            if term.getPower() > self.getPower():
                diff = term.getPower() - self.getPower()

                new = term.toArray()
                old = [0] * diff + self.toArray()

                new = [x + y for (x, y) in zip(new, old)]

            else:
                diff = self.getPower() - term.getPower()
                new = [0] * diff + term.toArray()
                new = [x + y for (x, y) in zip(new, self.toArray())]

            return Polynomial(arrayToTerm(new))

    def __mul__(self, other: Term | Polynomial) -> Polynomial:

        if isNumeric(other):
            arr = self.toArray()
            for t in range(len(arr)):
                arr[t] = arr[t]*other
            return Polynomial(arrayToTerm(arr))

        if isinstance(other, Term):
            arr = self.toArray()
            for i in range(len(arr)):
                arr[i] *= other.getCoefficient()
            arr += [0] * other.getPower()
            return Polynomial(arrayToTerm(arr))

        elif isinstance(other, Polynomial):

            output = Polynomial()
            for term in other._terms:
                output += self * term
            return output

    def __pow__(self, num) -> Polynomial:
        if isinstance(num, int):
            if num == 0:
                return Polynomial([Term(1, 0)])
            output = self
            for _ in range(num - 1):
                output = self * output
            return output
        else:
            raise Exception(f"Cannot raise polynomial to {num} power")

    def __truediv__(self, other) -> Polynomial | Rational:
        if isNumeric(other) and other != 0:
            return self * (1/other)
        if isinstance(other, Term):
            return Rational(self, Polynomial([other]))
        if isinstance(other, Polynomial):
            return Rational(self, other)

    def __neg__(self) -> Polynomial:
        out = Polynomial()
        for term in self._terms:
            out -= term.getCoefficient()
        return out

    def __str__(self) -> str:
        if max(self._terms) == 0 and min(self._terms) == 0:
            return "0"
        f = 0
        i = 0
        while f == 0:
            f = self._terms[i]
            output = str(f)
            i += 1
        for term in (self._terms[i:]):
            if term.getCoefficient() != 0:
                if isinstance(term.getCoefficient(), complex):
                    sign = "+" if (term.getCoefficient().real > 0 or term.getCoefficient(
                    ).real == 0 and term.getCoefficient().imag > 0) else ""
                else:
                    sign = "+" if (term.getCoefficient() > 0) else ""
                output += f"{sign}{str(term)}"
        return output

    def __lt__(self, other: Term | Polynomial) -> bool:
        if isPoly(other):
            return self.getPower() < other.getPower()

    def __eq__(self, other: Term | Polynomial) -> bool:
        if isPoly(other):
            return self.toArray() == other.toArray()

    def __gt__(self, other: Term | Polynomial) -> bool:
        if isPoly(other):
            return self.getPower() > other.getPower()

    def toArray(self) -> list[complex]:
        pad = (self._terms[0].getPower()-len(self._terms)+1)
        return [term.getCoefficient() for term in self._terms] + pad * [0]

    def evaluate(self, a: int | float | complex) -> float:
        return sum(map(lambda i: i.evaluate(a), self._terms))

    def derivative(self) -> Polynomial:
        return Polynomial(list(map(Term.derivative, self._terms[:-1])))

    def getRoots(self) -> list[np.float64 | np.complex64]:
        return list(map(lambda x: roundComplex(x, 3), np.roots(self.toArray())))

    def compose(self, p) -> int | float | complex | Polynomial:
        if isNumeric(p):
            return self.evaluate(p)
        if isinstance(p, Term) and p.getPower() == 1 and p.getCoefficient() == 1:
            return self
        elif isPoly(p):
            out = Polynomial()
            for term in self._terms:
                out += term.getCoefficient() * p ** term.getPower()
            return out
        elif isinstance(p, Rational):
            raise "Cannot compose Polynomials with Rationals"

    def real(self) -> Polynomial:
        p = Polynomial()
        for t in self._terms:
            p += t.real()
        return p

    def imag(self) -> Polynomial:
        p = Polynomial()
        for t in self._terms:
            p += t.imag()
        return p

    def isZero(self) -> bool:
        return self.toArray() == [0]


class Rational(Function):
    """
    Defines a rational function f(x) = h(x) / g(x) where g(x) != 0
    """

    def __init__(self, num: Polynomial = Polynomial([1]),
                 denom: Polynomial = Polynomial([1])) -> None:
        """
        Constructs a rational function using two polynomial objects
        """
        if isNumeric(num):
            num = Polynomial([num])
        if isNumeric(denom):
            denom = Polynomial([denom])
        self._num = num
        self._denom = denom
        if denom == Polynomial([0]):
            raise Exception("Cannot divide by 0")

    def getDenom(self) -> Polynomial:
        return self._denom

    def getNum(self) -> Polynomial:
        return self._num

    def setNum(self, newNum) -> None:
        self._num = newNum

    def setDenom(self, newDenom) -> None:
        self._denom = newDenom

    def __add__(self, poly) -> Rational:
        if isPoly(poly) or isNumeric(poly):
            return Rational(self.getNum() + poly * self.getDenom(), self.getDenom())
        elif isinstance(poly, Rational) and self.getDenom() == poly.getDenom():
            return Rational(self.getNum() + poly.getNum(), self.getDenom())
        elif isinstance(poly, Rational) and self.getDenom() != poly.getDenom():
            return Rational(self.getNum() * poly.getDenom() +
                            poly.getNum() * self.getDenom(),
                            self.getDenom() * poly.getDenom())

    def __mul__(self, poly) -> Rational:
        if isPoly(poly):
            return Rational(self.getNum() * poly, self.getDenom())
        elif isinstance(poly, Rational):
            return Rational((self.getNum() * poly.getNum()),
                            self.getDenom() * poly.getDenom())
        elif isNumeric(poly):
            return Rational(self.getNum() * poly, self.getDenom())

    def __pow__(self, p: int) -> Function:
        if p == 0:
            return Polynomial([1])
        return Rational(self.getNum()**p, self.getDenom()**p)

    def evaluate(self, a) -> float | complex:
        if self.getDenom().getRoots() != 0:
            return self.getNum().evaluate(a)/self.getDenom().evaluate(a)
        else:
            raise Exception("Cannot divide by 0")

    def compose(self, p) -> float | complex | Rational:
        if isPoly(p):
            return Rational(self.getNum().compose(p), self.getDenom().compose(p))
        elif isNumeric(p):
            return self.evaluate(p)

    def getPoles(self) -> list[np.float64 | np.complex64]:
        c = Counter(self.getDenom().getRoots())
        c.subtract(Counter(self.getNum().getRoots()))
        return list(c.elements())

    def __str__(self) -> str:
        return f"({self.getNum()})/({self.getDenom()})"

    def derivative(self) -> Rational:
        num = self.getNum()
        den = self.getDenom()
        return Rational((den * num.derivative() - num * den.derivative()), den**2)

    def __truediv__(self, other) -> Rational:
        if isNumeric(other) and other != 0:
            return self * (1/other)
        if isinstance(other, Term):
            return Rational(self.getNum(), self.getDenom()*Polynomial([other]))
        if isinstance(other, Polynomial):
            return Rational(self.getNum(), self.getDenom()*other)

    def simplify(self) -> Polynomial | Rational:
        if self.getNum().isZero() and not self.getDenom().isZero():
            return Polynomial()

        if self.getDenom().getPower() == 0 and self.getDenom().evaluate(1) != 0:
            return self.getNum() * 1/self.getDenom().evaluate(1)

        numRoots = self.getNum().getRoots()
        uniqueNumRoots = list(set(numRoots))
        denomRoots = self.getDenom().getRoots()
        uniqueDenRoots = list(set(denomRoots))
        intersection = [x for x in numRoots if x in denomRoots]
        if intersection != []:
            p1 = Polynomial([1])
            p2 = Polynomial([1])
            for root in uniqueNumRoots:
                p1 *= Polynomial([1, -root])**(numRoots.count(root) -
                                               denomRoots.count(root))
            for root in uniqueDenRoots:
                p2 *= Polynomial([1, -root]
                                 )**(denomRoots.count(root) - numRoots.count(root))
            return Rational(p1, p2)
        return self

    def real(self) -> Rational:
        a, b = self.getNum().real(), self.getNum().imag()
        c, d = self.getDenom().real(), self.getDenom().imag()
        if d.isZero():
            return Rational(a, c)
        return Rational(a*c-b*d, c**2 + d**2)

    def imag(self) -> Rational:
        a, b = self.getNum().real(), self.getNum().imag()
        c, d = self.getDenom().real(), self.getDenom().imag()
        if c.isZero():
            return Rational(b, d)
        return Rational(b*c - a*d, c**2 + d**2)

    def limit(self, a: complex) -> complex:
        if self.getDenom().evaluate(a) != 0:
            return self.evaluate(a)
        elif (abs(self.getDenom().evaluate(a)) <= 10**(-10)
                and abs(self.getNum().evaluate(a)) <= 10**(-10)):
            return Rational(self.getNum().derivative(),
                            self.getDenom().derivative()).limit(a)
        else:
            raise Exception(f"Division by 0 {self.getNum().evaluate(a)}/0")


sine = Polynomial([(-1)**((i % 4) + ((i+1)/2))/factorial(i) if i %
                   2 == 1 else 0 for i in range(101, -1, -1)])
cosine = Polynomial([(-1)**((i % 4) + (i/2))/factorial(i) if i %
                    2 == 0 else 0 for i in range(101, -1, -1)])
exp = Polynomial([1/factorial(i) for i in range(101, -1, -1)])

x = Term(1, 1)
sinh = (exp - exp.compose(-1*x)) / 2
cosh = (exp + exp.compose(-1*x)) / 2

eix: Polynomial = cosine + 1j * sine
