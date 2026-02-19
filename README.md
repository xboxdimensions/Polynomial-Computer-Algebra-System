I built a Computer Algebra System to handle polynomials and rational functions (as I didn't want to use existing libraries ie [SymPy](https://www.sympy.org/en/index.html])). This was for a numerical methods course, yet ironically - I tried to use as much analytic means to avoid as much numerical approximations as possible.

Complex Analysis was used - specifically Cauchy's Residue Theorem - $\oint_{C} f(z) dz = 2\pi i \sum Res(f,a)$ - which allowed for computation of improper integrals (the main purpose of this project). 

Furthermore, we exploit the fact that we have the following isomorphism between $ℂ$ and $ℝ^{2}$: 

<p align="center">
 $x+yi \to (x,y)$
</p>
which allows us to use Vector Calculus methods (ie Line Integrals) when integrating. This was done via NumPy and it's ability to vectorise arrays.

`main.py` has an example of parameterising a 2D contour and integrating $f(z) = \frac{1}{z^2 + 1}$ around the upper semi-circle of the Complex Plane.

Currently, I only support 2D vector fields, however `VectorField.py` can be easily extended to higher dimensions if you so choose.

The functions $z$, $\sin{z}$, $\cos{z}$, $e^{z}$, $\sinh{z}$, $\cosh{z}$ and $e^{iz}$ are available and predefined as per `Functions.py` 
given they are entire (converges to their power series over all of $ℂ$). You can probably make other functions but I wanted to integrate over the entire Real line so I required my series to converge everywhere. Note, that these approximations are best at $z=0$ so your mileage may vary when $|z|$ increases.
If you choose to add more terms for the transcendental functions - modify the starting term (default is 101).
