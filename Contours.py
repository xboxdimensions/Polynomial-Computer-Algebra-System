from __future__ import annotations
import numpy as np
import scipy.integrate
from LineIntegral import lineInt
from Functions import Function
from VectorField import VectorField
from matplotlib.path import Path
import matplotlib.pyplot as plt
from Residues import getOrder, getResidue


class Contour():
    """
    A region (or just the bounds) of the curve provided
    """

    def __init__(self, Trajectory: list[callable], endPoints: list[float]) -> None:
        """
        Initialises a piecewise contour defined by the Trajectory.
        - Precondition: endPoints must be sorted ascending
        """
        self._trajectory = Trajectory
        self._endpoints = endPoints
        self._start = endPoints[0]
        self._end = endPoints[-1]

    def getEndPoints(self) -> list[float]:
        return self._endpoints

    def _getTrajectory(self, t: float | np.ndarray) -> list[list[float]]:
        endpoints = self.getEndPoints()
        if isinstance(t, np.ndarray):
            res = []
            for ts in t:
                for i in range(len(endpoints)-1):
                    if endpoints[i] <= ts and endpoints[i+1] >= ts:
                        res.append(self._trajectory[i](ts))
                        break
            return np.array(res)
        else:
            for i in range(len(endpoints)-1):
                if endpoints[i] <= t and endpoints[i+1] >= t:
                    return self._trajectory[i](t)

    def computeContourIntegral(self, f: Function, force: VectorField) -> float:
        return self.computeCurveIntegral(force) + self.computeLineIntegral(f)

    def computeLineIntegral(self, f: Function) -> float:
        start = self._getTrajectory(self._start)
        end = self._getTrajectory(self._end)

        a1, b1 = start
        c1, d1 = end

        # It is just a point (or it's closed)
        if ((a1-c1)**2 + (b1-d1)**2) < 10**(-8):
            return 0

        # Suppose that both the start and end points are on the Real Line.

        xs = np.linspace(c1, a1)
        ys = list(map(f.evaluate, xs))
        return scipy.integrate.simpson(ys, xs)

    def computeCurveIntegral(self, force: VectorField, num: int = 100) -> float:
        endpoints = self.getEndPoints()
        sum = 0
        for i in range(len(endpoints) - 1):
            start = endpoints[i]
            end = endpoints[i+1]
            int = lineInt(start, end, force,
                          self._trajectory[i], round(num/len(endpoints)))
            sum += int
        return sum

    def sumRes(self, f: Function) -> float | complex:
        """
        Computes the contour integral via residues - via Cauchy's Residue Theorem
        """
        sum = 0
        poles = getOrder(f)
        for pole in poles:
            sum += getResidue(f, pole) if self.contains(pole) else 0
        return np.pi*2j*sum

    def collectPoints(self, points: list[list[float]]) -> list[tuple[float, float]]:
        return [(points[i][0], points[i][1]) for i in range(len(points))]

    def contains(self, a: complex) -> bool:
        shape = Path(self.collectPoints(self.generatePoints(100)))
        return shape.contains_point((a.real, a.imag))

    def generatePoints(self, num: int) -> list[list[float]]:
        endpoints = self.getEndPoints()
        total = np.array([])
        for i in range(len(endpoints)-1):
            ts = np.linspace(endpoints[i], endpoints[i+1],
                             num//len(self._trajectory), dtype=np.float64)
            total.resize(len(total) + len(ts))
            total[-len(ts):] = ts
        return [self._getTrajectory(t) for t in total]

    def plotContour(self, f: Function = None) -> None:
        points = self.collectPoints(self.generatePoints(100))
        re, im = [x for (x, _) in points], [y for (_, y) in points]
        plt.plot(re, im)
        a, b = self._getTrajectory(self._start)
        c, d = self._getTrajectory(self._end)
        if np.sqrt((a-c)**2 + (b-d)**2) > 10**(-5):
            plt.plot([a, c], [b, d])
        if f:
            real = []
            imag = []
            poles = set(map(complex, f.getPoles()))
            for pole in poles:
                if self.contains(pole):
                    real.append(pole.real)
                    imag.append(pole.imag)
            plt.scatter(real, imag)
        plt.xlabel("Re(z)")
        plt.ylabel("Im(z)")
        plt.show()
