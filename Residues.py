from Functions import Function, Polynomial, Rational
from math import factorial


def getOrder(f: Rational) -> dict[complex, int]:
    order = dict()
    for pole in f.getPoles():
        order[pole] = 1 if not order.get(pole) else order.get(pole) + 1
    return order


def getResidue(f: Function, a: complex) -> complex:
    if isinstance(f, Rational):
        info = getOrder(f)
        if a not in info.keys():
            return 0
        order = info[a]
        p = (f/factorial(order - 1)) * Polynomial([1, -a])**order
        [p := p.derivative() for _ in range(order - 1)]
        p = p.simplify()
        return p.limit(a)
    return 0
