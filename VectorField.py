from numpy import array


class VectorField():
    """
    A 2 Dimensional Vector Field
    """

    def __init__(self, f1: callable, f2: callable) -> None:
        """
        Constructs F(x,y) = (f1(x,y), f2(x,y))
        - Precondition: f1, f2 must be functions that take in exactly 2 inputs
        """
        self._f1 = f1
        self._f2 = f2

    def getF1(self, x, y) -> callable:
        return self._f1(x, y)

    def getF2(self, x, y) -> callable:
        return self._f2(x, y)

    def generatePoints(self, rs: list[list[float]]) -> list[list[float]]:
        """
        Maps the input of a 2D array by the Vector Function F
        """
        vals = []
        for i in range(len(rs[0])):
            x, y = rs[0][i], rs[1][i]

            val = [self.getF1(x, y), self.getF2(x, y)]

            vals.append(val)
        return array(vals)
