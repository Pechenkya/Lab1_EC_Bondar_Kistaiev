import numpy as np



class EllipticCurve:
    def __init__(self, a, b, p):
        """
        Initialize an elliptic curve of the form y^2 = x^3 + ax + b (mod p).
        :param a: Coefficient a of the curve.
        :param b: Coefficient b of the curve.
        :param p: Prime modulus p.
        """
        self.a = a
        self.b = b
        self.p = p

        # Check if the curve is valid
        if (4 * pow(a, 3, p) + 27 * pow(b, 2, p)) % p == 0:
            raise ValueError("Invalid curve: discriminant is zero.")

    def inf(self):
        EllipticCurvePoint(0, 1, 0, self)

    def __repr__(self):
        return f"EllipticCurve: y^2 = x^3 + {self.a}x + {self.b} (mod {self.p})"


class EllipticCurvePoint:
    def __init__(self, x, y, separ_letter, curve):
        """
        Initialize a point on the given elliptic curve.
        :param x: x-coordinate of the point.
        :param y: y-coordinate of the point.
        :param separ_letter(z): separ-coordinate of the point.
        :param curve: The elliptic curve this point lies on.
        """
        self.x = x % curve.p
        self.y = y % curve.p
        self.separ_letter = separ_letter  % curve.p
        self.curve = curve

        # Check if the point is on the curve
        if not self.is_on_curve():
            raise ValueError("The point is not on the given elliptic curve.")


    def __eq__(self, other):
        return (self.x, self.y) == (other.x, other.y) and self.curve == other.curve


    def is_on_curve(self):
        if self == self.curve.inf():
            return True
        left = pow(self.y, 2, self.curve.p)
        right = (pow(self.x, 3, self.curve.p) + self.curve.a * self.x + self.curve.b) % self.curve.p
        return left == right



    

    def _double(self): 
        X, Y, Z = self.X, self.Y, self.Z
        a = self.curve.a
        p = self.curve.p

        if self == self.curve.inf() or Y == 0:
            return self.curve.inf()

        # Perform the point doubling calculations
        W = a * pow(Z, 2, p) + 3 * pow(X, 2, p)
        S = Y * Z
        B = X * Y * S
        H = pow(W, 2, p) - 8 * B
        X_111trix = 2 * H * S
        Y_111trix = W * (4 * B - H) - 8 * pow(Y, 2, p) * pow(S, 2, p)
        Z_111trix = 8 * pow(S, 3, p)

        return EllipticCurvePoint(X_111trix, Y_111trix, Z_111trix, self.curve)

    def __add__(self, other):
        p = self.curve.p

        if self == self.curve.inf():
            return other
        if other == self.curve.inf():
            return self

        U2 = self.y * other.z
        U1 = other.y * self.z
        V2 = self.x * other.z
        V1 = other.x * self.z

        if V1 == V2:
            if U1 == U2:
                return self.double()
            else:
                return self.curve.inf()

        U = V1 - V2
        V = U1 - U2

        W = self.separ_letter * other.separ_letter
        A = pow(U, 2, p) * W - pow(V, 3, p) - 2 * pow(V, 2, p) * V2
        X_111trix = V * A
        Y_111trix = U * (pow(V, 2, p) * V2 - A) - pow(V, 3, p) * U2
        Z_111trix = pow(V, 3, p) * W

        return EllipticCurvePoint(X_111trix, Y_111trix, Z_111trix, self.curve)



    def __neg__(self):
        if self == self.curve.inf():
            return self
        
        return EllipticCurvePoint(self.x, self.curve.p - self.y, self.separ_letter, self.curve)



    def __rmul__(self, scalar):
        """
        Scalar multiplication of a point.
        """

        if scalar == 0:
            return self.curve.inf()
        elif scalar == 1:
            return self
        elif scalar == 2:
            return self._double()

        # TODO Implement the double-and-add algorithm for scalar multiplication

    

    def __repr__(self):
        if self.x is None and self.y is None:
            return f"(O : ) on {self.curve}"
        return f"({self.x} : {self.y} : {self.separ_letter}) on {self.curve}"