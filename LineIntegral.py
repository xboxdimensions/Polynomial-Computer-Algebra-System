import numpy as np
import scipy.integrate
from VectorField import VectorField


def deriv(f: callable, x):
    h = 1.49*10**(-8)
    return (f(x+h)-f(x))/h


# def r(t): return np.array([t, 4*t])


# def f(t):
#     x, y = t
#     return x**2 * y

# Integrate f over a curve


def trapezoidRule(t: list[float], y: list[float]) -> float:
    return sum([(t[n]-t[n-1]) * (y[n]+y[n-1])/2 for n in range(1, len(t))])


def scalarInt(a: float, b: float, f: callable, r: callable, size: int) -> float:
    t: list[float] = np.linspace(a, b, size, dtype=np.float64)
    mag = np.sqrt(sum(map(lambda x: x**2, deriv(r, t))))
    y: list[float] = list(map(lambda x: x[0]*x[1], zip(f(r(t)), mag)))
    return scipy.integrate.simpson(y, t)


# print(scalarInt(0, 1, f, r, 10000))


def lineInt(a: float, b: float, F: VectorField, r: np.array, size: int) -> float:
    t: list[float] = np.linspace(a, b, size)
    direction = deriv(r, t).T
    force = F.generatePoints(r(t))
    y: list[float] = [np.dot(force[i], direction[i]) for i in range(size)]
    return scipy.integrate.simpson(y, t)
