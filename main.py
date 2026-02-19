import numpy as np
from Contours import Contour
from Residues import *
from Functions import *
from LineIntegral import *
from VectorField import VectorField
import matplotlib.pyplot as plt

R = 300

endPoints = [0, np.pi]
components = [lambda t: np.array(
    [R*np.cos(t), R*np.sin(t)])]

semiCircle = Contour(components, endPoints)
f = Rational(1, Polynomial([1, 0, 1]))
def f1(x, y): return f.evaluate(x+y*1j).real
def f2(x, y): return f.evaluate(x+y*1j).imag


vec = VectorField(f1, f2)
print(semiCircle.sumRes(f))
semiCircle.plotContour(f)