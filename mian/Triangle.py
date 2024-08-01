import math
import numpy as np
from sympy import *
import time
import warnings


class CircumferenceAndArea(object):
    """用于计算三角形周长和面积"""

    def __init__(self, a=0, b=0, c=0, bottom=0, high=0, pax=0, pbx=0, pcx=0, pay=0, pby=0, pcy=0, paz=0, pbz=0, pcz=0,
                 getparms=0, dimensionality=0):
        self.SidesA = a
        self.SidesB = b
        self.SidesC = c
        self.bottom = bottom
        self.high = high
        self.CoordinatePointAX = pax
        self.CoordinatePointBX = pbx
        self.CoordinatePointCX = pcx
        self.CoordinatePointAY = pay
        self.CoordinatePointBY = pby
        self.CoordinatePointCY = pcy
        self.CoordinatePointAZ = paz
        self.CoordinatePointBZ = pbz
        self.CoordinatePointCZ = pcz
        self.Getparms = getparms
        self.Dimensionality = dimensionality

    def circumference(self):
        circumference = self.SidesA + self.SidesB + self.SidesC
        return circumference

    def area_is_bottom_high(self):
        area = self.bottom * self.high * 0.5
        return area

    def area_is_sides(self):
        p = (self.SidesA + self.SidesB + self.SidesC) * 0.5
        area = math.sqrt(p * (p - self.SidesA) * (p - self.SidesB) * (p - self.SidesC))
        return area

    def area_is_planar_vector(self):
        vector_abx = self.CoordinatePointAX - self.CoordinatePointBX
        vector_aby = self.CoordinatePointAY - self.CoordinatePointBY
        vector_acx = self.CoordinatePointAX - self.CoordinatePointCX
        vector_acy = self.CoordinatePointAY - self.CoordinatePointCY
        area = abs(0.5 * (1 * vector_abx * vector_acy) + (-1 * vector_aby * vector_acx))
        return area

    def area_is_spatial_vectors(self):
        vector_ab = np.array([self.CoordinatePointBX - self.CoordinatePointAX,
                              self.CoordinatePointBY - self.CoordinatePointAY,
                              self.CoordinatePointBZ - self.CoordinatePointAZ])
        vector_ac = np.array([self.CoordinatePointCX - self.CoordinatePointAX,
                              self.CoordinatePointCY - self.CoordinatePointAY,
                              self.CoordinatePointCZ - self.CoordinatePointAZ])
        area = 0.5 * np.linalg.norm(np.cross(vector_ab, vector_ac))
        return area


class TriangleCentres(CircumferenceAndArea):
    """用于求三角形四心"""

    def centroid(self):
        if self.Dimensionality == 2:
            planar_centroid_x = (self.CoordinatePointAX + self.CoordinatePointBX + self.CoordinatePointCX) / 3
            planar_centroid_y = (self.CoordinatePointAY + self.CoordinatePointBY + self.CoordinatePointCY) / 3
            planar_centroid_xy = (planar_centroid_x, planar_centroid_y)
            return planar_centroid_xy
        elif self.Dimensionality == 3:
            planar_centroid_x = (self.CoordinatePointAX + self.CoordinatePointBX + self.CoordinatePointCX) / 3
            planar_centroid_y = (self.CoordinatePointAY + self.CoordinatePointBY + self.CoordinatePointCY) / 3
            planar_centroid_z = (self.CoordinatePointAZ + self.CoordinatePointBZ + self.CoordinatePointCZ) / 3
            planar_centroid_xyz = (planar_centroid_x, planar_centroid_y, planar_centroid_z)
            return planar_centroid_xyz
        else:
            warnings.warn("Please enter whether it is 2D or 3D", SyntaxWarning)

    def incentre(self):
        if self.Dimensionality == 2:
            a = math.sqrt((self.CoordinatePointBX - self.CoordinatePointCX) ** 2 + (
                    self.CoordinatePointBY - self.CoordinatePointCY) ** 2)
            b = math.sqrt((self.CoordinatePointCX - self.CoordinatePointAX) ** 2 + (
                    self.CoordinatePointCY - self.CoordinatePointAY) ** 2)
            c = math.sqrt((self.CoordinatePointAX - self.CoordinatePointBX) ** 2 + (
                    self.CoordinatePointAY - self.CoordinatePointBY) ** 2)
            planar_incentre_x = (
                                        a * self.CoordinatePointAX + b * self.CoordinatePointBX + c * self.CoordinatePointCX) / (
                                        a + b + c)
            planar_incentre_y = (
                                        a * self.CoordinatePointAY + b * self.CoordinatePointBY + c * self.CoordinatePointCY) / (
                                        a + b + c)
            planar_incentre_xy = (planar_incentre_x, planar_incentre_y)
            return planar_incentre_xy

        elif self.Dimensionality == 3:
            a = math.sqrt((self.CoordinatePointBX - self.CoordinatePointCX) ** 2 + (
                    self.CoordinatePointBY - self.CoordinatePointCY) ** 2 + (
                                  self.CoordinatePointBZ - self.CoordinatePointCZ) ** 2)
            b = math.sqrt((self.CoordinatePointCX - self.CoordinatePointAX) ** 2 + (
                    self.CoordinatePointCY - self.CoordinatePointAY) ** 2 + (
                                  self.CoordinatePointCZ - self.CoordinatePointAZ) ** 2)
            c = math.sqrt((self.CoordinatePointAX - self.CoordinatePointBX) ** 2 + (
                    self.CoordinatePointAY - self.CoordinatePointBY) ** 2 + (
                                  self.CoordinatePointAZ - self.CoordinatePointBZ) ** 2)

            spatial_incentre_x = (
                                         a * self.CoordinatePointAX + b * self.CoordinatePointBX + c * self.CoordinatePointCX) / (
                                         a + b + c)
            spatial_incentre_y = (
                                         a * self.CoordinatePointAY + b * self.CoordinatePointBY + c * self.CoordinatePointCY) / (
                                         a + b + c)
            spatial_incentre_z = (
                                         a * self.CoordinatePointAZ + b * self.CoordinatePointBZ + c * self.CoordinatePointCZ) / (
                                         a + b + c)
            spatial_incentre_xyz = (spatial_incentre_x, spatial_incentre_y, spatial_incentre_z)
            return spatial_incentre_xyz
        else:
            warnings.warn("Please enter whether it is 2D or 3D", SyntaxWarning)

    def orthocentre(self):
        if self.Dimensionality == 2:
            a = (self.CoordinatePointAX - self.CoordinatePointBX) * (
                    self.CoordinatePointBX - self.CoordinatePointCX) * (
                        self.CoordinatePointCX - self.CoordinatePointAX)
            b = (self.CoordinatePointAY - self.CoordinatePointBY) * (
                    self.CoordinatePointBY - self.CoordinatePointCY) * (
                        self.CoordinatePointCY - self.CoordinatePointAY)
            c = self.CoordinatePointAX * self.CoordinatePointBX * (
                    self.CoordinatePointAY - self.CoordinatePointBY) + self.CoordinatePointBX * self.CoordinatePointCX * (
                        self.CoordinatePointBY - self.CoordinatePointCY) + self.CoordinatePointCX * self.CoordinatePointAX * (
                        self.CoordinatePointCY - self.CoordinatePointAY)
            d = self.CoordinatePointAY * self.CoordinatePointBY * (
                    self.CoordinatePointAX - self.CoordinatePointBX) + self.CoordinatePointBY * self.CoordinatePointCY * (
                        self.CoordinatePointBX - self.CoordinatePointCX) + self.CoordinatePointCY * self.CoordinatePointAY * (
                        self.CoordinatePointCX - self.CoordinatePointAX)
            e = (self.CoordinatePointAX * self.CoordinatePointCY - self.CoordinatePointAX * self.CoordinatePointBY) + (
                    self.CoordinatePointBX * self.CoordinatePointAY - self.CoordinatePointBX * self.CoordinatePointCY) + (
                        self.CoordinatePointCX * self.CoordinatePointBY - self.CoordinatePointCX * self.CoordinatePointAY)
            planar_orthocentre_x = (b + c) / e
            planar_orthocentre_y = (a + d) / e
            planar_orthocentre_xy = (planar_orthocentre_x, planar_orthocentre_y)
            return planar_orthocentre_xy
        elif self.Dimensionality == 3:
            a = (self.CoordinatePointBX - self.CoordinatePointAX) * (
                    self.CoordinatePointCX - self.CoordinatePointAX) + (
                        self.CoordinatePointBY - self.CoordinatePointAY) * (
                        self.CoordinatePointCY - self.CoordinatePointAY) + (
                        self.CoordinatePointBZ - self.CoordinatePointAZ) * (
                        self.CoordinatePointCZ - self.CoordinatePointAZ)
            b = (self.CoordinatePointAX - self.CoordinatePointBX) * (
                    self.CoordinatePointCX - self.CoordinatePointBX) + (
                        self.CoordinatePointAY - self.CoordinatePointBY) * (
                        self.CoordinatePointCY - self.CoordinatePointBY) + (
                        self.CoordinatePointAZ - self.CoordinatePointBZ) * (
                        self.CoordinatePointCZ - self.CoordinatePointBZ)
            c = (self.CoordinatePointAX - self.CoordinatePointCX) * (
                    self.CoordinatePointBX - self.CoordinatePointCX) + (
                        self.CoordinatePointAY - self.CoordinatePointCY) * (
                        self.CoordinatePointBY - self.CoordinatePointCY) + (
                        self.CoordinatePointAZ - self.CoordinatePointCZ) * (
                        self.CoordinatePointBZ - self.CoordinatePointCZ)
            spatial_orthocentre_x = ((
                                             b * c * self.CoordinatePointAX + a * c * self.CoordinatePointBX + a * b * self.CoordinatePointCX) / (
                                             b * c + a * c + a * b))
            spatial_orthocentre_y = ((
                                             b * c * self.CoordinatePointAY + a * c * self.CoordinatePointBY + a * b * self.CoordinatePointCY) / (
                                             b * c + a * c + a * b))
            spatial_orthocentre_z = ((
                                             b * c * self.CoordinatePointAZ + a * c * self.CoordinatePointBZ + a * b * self.CoordinatePointCZ) / (
                                             b * c + a * c + a * b))
            spatial_orthocentre_xyz = (spatial_orthocentre_x, spatial_orthocentre_y, spatial_orthocentre_z)
            return spatial_orthocentre_xyz
        else:
            warnings.warn("Please enter whether it is 2D or 3D", SyntaxWarning)

    def circumcentre(self):
        if self.Dimensionality == 2:
            a = ((self.CoordinatePointBX ** 2) + (self.CoordinatePointBY ** 2)) - (
                    (self.CoordinatePointAX ** 2) + (self.CoordinatePointAY ** 2))
            b = ((self.CoordinatePointCX ** 2) + (self.CoordinatePointCY ** 2)) - (
                    (self.CoordinatePointBX ** 2) + (self.CoordinatePointBY ** 2))
            fm = 2 * ((self.CoordinatePointCY - self.CoordinatePointBY) * (
                    self.CoordinatePointBX - self.CoordinatePointAX) - (
                              self.CoordinatePointBY - self.CoordinatePointAY) * (
                              self.CoordinatePointCX - self.CoordinatePointBX))
            planar_circumcentre_x = ((self.CoordinatePointCY - self.CoordinatePointBY) * a - (
                    self.CoordinatePointBY - self.CoordinatePointAY) * b) / fm
            planar_circumcentre_y = ((self.CoordinatePointBX - self.CoordinatePointAX) * b - (
                    self.CoordinatePointCX - self.CoordinatePointBX) * a) / fm
            planar_circumcentre_xy = (planar_circumcentre_x, planar_circumcentre_y)
            return planar_circumcentre_xy

        elif self.Dimensionality == 3:
            if self.CoordinatePointAZ == self.CoordinatePointBZ == self.CoordinatePointCZ:
                spatial_circumcentre_z = self.CoordinatePointAZ
                a = ((self.CoordinatePointBX ** 2) + (self.CoordinatePointBY ** 2)) - (
                        (self.CoordinatePointAX ** 2) + (self.CoordinatePointAY ** 2))
                b = ((self.CoordinatePointCX ** 2) + (self.CoordinatePointCY ** 2)) - (
                        (self.CoordinatePointBX ** 2) + (self.CoordinatePointBY ** 2))
                fm = 2 * ((self.CoordinatePointCY - self.CoordinatePointBY) * (
                        self.CoordinatePointBX - self.CoordinatePointAX) - (
                                  self.CoordinatePointBY - self.CoordinatePointAY) * (
                                  self.CoordinatePointCX - self.CoordinatePointBX))
                spatial_circumcentre_x = ((self.CoordinatePointCY - self.CoordinatePointBY) * a - (
                        self.CoordinatePointBY - self.CoordinatePointAY) * b) / fm
                spatial_circumcentre_y = ((self.CoordinatePointBX - self.CoordinatePointAX) * b - (
                        self.CoordinatePointCX - self.CoordinatePointBX) * a) / fm
                spatial_circumcentre_xyz = (spatial_circumcentre_x, spatial_circumcentre_y, spatial_circumcentre_z)
                return spatial_circumcentre_xyz
            else:
                a = np.array([self.CoordinatePointAX, self.CoordinatePointAY, self.CoordinatePointAZ])
                b = np.array([self.CoordinatePointBX, self.CoordinatePointBY, self.CoordinatePointBZ])
                c = np.array([self.CoordinatePointCX, self.CoordinatePointCY, self.CoordinatePointCZ])
                a_dot_b = np.dot(a, b)
                b_dot_c = np.dot(b, c)
                c_dot_a = np.dot(c, a)
                a_modulus = np.linalg.norm(a)
                b_modulus = np.linalg.norm(b)
                c_modulus = np.linalg.norm(c)
                cos_a = a_dot_b / a_modulus * b_modulus
                cos_b = c_dot_a / c_modulus * a_modulus
                cos_c = b_dot_c / b_modulus * c_modulus
                sin_2a = 2 * (cos_a - 90)
                sin_2b = 2 * (cos_b - 90)
                sin_2c = 2 * (cos_c - 90)
                spatial_circumcentre_x = (
                                                 sin_2a * self.CoordinatePointAX + sin_2b * self.CoordinatePointBX + sin_2c * self.CoordinatePointCX) / (
                                                 sin_2a + sin_2b + sin_2c)
                spatial_circumcentre_y = (
                                                 sin_2a * self.CoordinatePointAY + sin_2b * self.CoordinatePointBY + sin_2c * self.CoordinatePointCY) / (
                                                 sin_2a + sin_2b + sin_2c)
                spatial_circumcentre_z = (
                                                 sin_2a * self.CoordinatePointAZ + sin_2b * self.CoordinatePointBZ + sin_2c * self.CoordinatePointCZ) / (
                                                 sin_2a + sin_2b + sin_2c)
                spatial_circumcentre_xyz = (spatial_circumcentre_x, spatial_circumcentre_y, spatial_circumcentre_z)
                return spatial_circumcentre_xyz
        else:
            warnings.warn("Please enter whether it is 2D or 3D", SyntaxWarning)


class FermatProblem(CircumferenceAndArea):
    """"费马问题"""

    def fermat_problem(self):
        if self.Dimensionality == 2:
            a = np.array([self.CoordinatePointAX, self.CoordinatePointAY])
            b = np.array([self.CoordinatePointBX, self.CoordinatePointBY])
            c = np.array([self.CoordinatePointCX, self.CoordinatePointCY])
            a_minus_b = a - b
            a_minus_c = a - c
            b_minus_c = b - c
            vector_modulus_a = np.linalg.norm(a_minus_b)
            vector_modulus_b = np.linalg.norm(a_minus_c)
            vector_modulus_c = np.linalg.norm(b_minus_c)
            if vector_modulus_a + vector_modulus_b > vector_modulus_c and vector_modulus_a + vector_modulus_c > vector_modulus_b and vector_modulus_b + vector_modulus_c > vector_modulus_a:
                a_b = np.dot(a_minus_b, a_minus_c)
                a_c = np.dot(a_minus_b, b_minus_c)
                b_c = np.dot(a_minus_c, b_minus_c)
                cos_a = a_b / (vector_modulus_a * vector_modulus_b)
                cos_b = a_c / (vector_modulus_a * vector_modulus_c)
                cos_c = b_c / (vector_modulus_b * vector_modulus_c)
                aa = np.degrees(np.arcsin(cos_a))
                bb = np.degrees(np.arcsin(cos_b))
                cc = np.degrees(np.arccos(cos_c))
                aaa = abs(aa)
                bbb = abs(bb)
                ccc = abs(cc)

                if aaa < 120 and bbb < 120 and ccc < 120:
                    p = (vector_modulus_a + vector_modulus_b + vector_modulus_c) * 0.5
                    area = math.sqrt(p * (p - vector_modulus_a) * (p - vector_modulus_b) * (p - vector_modulus_c))
                    k_2 = (vector_modulus_a ** 2 + vector_modulus_b ** 2 + vector_modulus_c ** 2) * 0.5 + 2 * math.sqrt(
                        3) * area
                    k = k_2 ** 0.5
                    y = (vector_modulus_a ** 2 + vector_modulus_b ** 2 - 2 * vector_modulus_c ** 2 + k_2) / (3 * k)
                    x = (vector_modulus_a ** 2 + vector_modulus_c ** 2 - 2 * vector_modulus_b ** 2 + k_2) / (3 * k)
                    z = (vector_modulus_c ** 2 + vector_modulus_b ** 2 - 2 * vector_modulus_a ** 2 + k_2) / (3 * k)
                    if self.Getparms == 1:
                        l, j = symbols('l j')
                        eq1 = Eq(((self.CoordinatePointCX - l) ** 2 + (self.CoordinatePointCY - j) ** 2) ** 0.5, x)
                        eq2 = Eq(((self.CoordinatePointBX - l) ** 2 + (self.CoordinatePointBY - j) ** 2) ** 0.5, z)
                        solution = solve((eq1, eq2), (l, j))
                        return solution
                    elif self.Getparms == 0:
                        xyz = x + y + z
                        return xyz
                    else:
                        warnings.warn("getparms value is useless", SyntaxWarning)
                if aaa >= 120 or bbb >= 120 or ccc >= 120:
                    max_a = max(aaa, bbb, ccc)
                    if max_a == aaa:
                        return a
                    elif max_a == bbb:
                        return c
                    elif max_a == ccc:
                        return b
            else:
                print("It is not a triangle")

        elif self.Dimensionality == 3:
            a = np.array([self.CoordinatePointAX, self.CoordinatePointAY, self.CoordinatePointAZ])
            b = np.array([self.CoordinatePointBX, self.CoordinatePointBY, self.CoordinatePointBZ])
            c = np.array([self.CoordinatePointCX, self.CoordinatePointCY, self.CoordinatePointCZ])
            a_minus_b = a - b
            a_minus_c = a - c
            b_minus_c = b - c
            vector_modulus_a = np.linalg.norm(a_minus_b)
            vector_modulus_b = np.linalg.norm(a_minus_c)
            vector_modulus_c = np.linalg.norm(b_minus_c)
            if vector_modulus_a + vector_modulus_b > vector_modulus_c and vector_modulus_a + vector_modulus_c > vector_modulus_b and vector_modulus_b + vector_modulus_c > vector_modulus_a:
                a_b = np.dot(a_minus_b, a_minus_c)
                a_c = np.dot(a_minus_b, b_minus_c)
                b_c = np.dot(a_minus_c, b_minus_c)
                cos_a = a_b / (vector_modulus_a * vector_modulus_b)
                cos_b = a_c / (vector_modulus_a * vector_modulus_c)
                cos_c = b_c / (vector_modulus_b * vector_modulus_c)
                aa = np.degrees(np.arcsin(cos_a))
                bb = np.degrees(np.arcsin(cos_b))
                cc = np.degrees(np.arccos(cos_c))
                aaa = abs(aa)
                bbb = abs(bb)
                ccc = abs(cc)

                if aaa < 120 and bbb < 120 and ccc < 120:
                    p = (vector_modulus_a + vector_modulus_b + vector_modulus_c) * 0.5
                    area = math.sqrt(p * (p - vector_modulus_a) * (p - vector_modulus_b) * (p - vector_modulus_c))
                    k_2 = (vector_modulus_a ** 2 + vector_modulus_b ** 2 + vector_modulus_c ** 2) * 0.5 + 2 * math.sqrt(
                        3) * area
                    k = k_2 ** 0.5
                    y = (vector_modulus_a ** 2 + vector_modulus_b ** 2 - 2 * vector_modulus_c ** 2 + k_2) / (3 * k)
                    x = (vector_modulus_a ** 2 + vector_modulus_c ** 2 - 2 * vector_modulus_b ** 2 + k_2) / (3 * k)
                    z = (vector_modulus_c ** 2 + vector_modulus_b ** 2 - 2 * vector_modulus_a ** 2 + k_2) / (3 * k)
                    if self.Getparms == 1:
                        if self.CoordinatePointAZ == self.CoordinatePointBZ == self.CoordinatePointCZ:
                            l, j = symbols('l j')
                            eq1 = Eq(((self.CoordinatePointCX - l) ** 2 + (self.CoordinatePointCY - j) ** 2) ** 0.5, x)
                            eq2 = Eq(((self.CoordinatePointBX - l) ** 2 + (self.CoordinatePointBY - j) ** 2) ** 0.5, z)
                            solution = solve((eq1, eq2), (l, j))
                            return solution
                        else:
                            l, j, f = symbols('l j f')
                            eq1 = Eq(((self.CoordinatePointCX - l) ** 2 + (self.CoordinatePointCY - j) ** 2 + (
                                    self.CoordinatePointCZ - f) ** 2) ** 0.5, x)
                            eq2 = Eq(((self.CoordinatePointBX - l) ** 2 + (self.CoordinatePointBY - j) ** 2 + (
                                    self.CoordinatePointBZ - f) ** 2) ** 0.5, z)
                            eq3 = Eq(((self.CoordinatePointAX - l) ** 2 + (self.CoordinatePointAY - j) ** 2 + (
                                    self.CoordinatePointAZ - f) ** 2) ** 0.5, y)
                            solution = solve((eq1, eq2, eq3), (l, j, f))
                            return solution
                    elif self.Getparms == 0:
                        xyz = x + y + z
                        return xyz
                    else:
                        warnings.warn("getparms value is useless", SyntaxWarning)
                if aaa >= 120 or bbb >= 120 or ccc >= 120:
                    max_a = max(aaa, bbb, ccc)
                    if max_a == aaa:
                        return a
                    elif max_a == bbb:
                        return c
                    elif max_a == ccc:
                        return b
            else:
                print("It is not a triangle")
        else:
            warnings.warn("Please enter whether it is 2D or 3D", SyntaxWarning)


class NapoleonicTriangle(CircumferenceAndArea):
    def napoleonic_triangle(self):
        if self.Dimensionality == 2:
            if (self.CoordinatePointAX == 0 and self.CoordinatePointAY == 0) or (self.CoordinatePointBX == 0 and self.CoordinatePointBY == 0) or (self.CoordinatePointCX == 0 and self.CoordinatePointCY == 0):
                a = np.array([self.CoordinatePointAX, self.CoordinatePointAY])
                u = np.array([self.CoordinatePointBX, self.CoordinatePointBY])
                v = np.array([self.CoordinatePointCX, self.CoordinatePointCY])
                w = u - v
                planar_x = (w[0] + w[0] + w[0]) / 3
                planar_y = (self.CoordinatePointAY + self.CoordinatePointBY + self.CoordinatePointCY) / 3
                planar_xy = (planar_x, planar_y)
                return planar_xy




t = time.time()
tri = TriangleCentres(pax=0, pay=0, pbx=8, pby=0, pcx=4, pcy=2, paz=0, pbz=4, pcz=2, dimensionality=3)
print(tri.circumcentre())
print(f'coast:{time.time() - t:.4f}s')
