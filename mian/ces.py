import numpy as np
from sympy import *


def planar_fermat_problem():
    CoordinatePointAX = 1
    CoordinatePointAY = 1
    CoordinatePointBX = 8
    CoordinatePointBY = 1
    CoordinatePointCX = 4
    CoordinatePointCY = 2
    vector_abx = CoordinatePointAX - CoordinatePointBX
    vector_aby = CoordinatePointAY - CoordinatePointBY
    vector_acx = CoordinatePointAX - CoordinatePointCX
    vector_acy = CoordinatePointAY - CoordinatePointCY
    vector_bcx = CoordinatePointBX - CoordinatePointCX
    vector_bcy = CoordinatePointBY - CoordinatePointCY
    vector_c = np.array([vector_abx, vector_aby])
    vector_a = np.array([vector_bcx, vector_bcy])
    vector_b = np.array([vector_acx, vector_acy])
    c_dot_a = np.dot(vector_c, vector_a)
    c_dot_b = np.dot(vector_c, vector_b)
    b_dot_a = np.dot(vector_b, vector_a)
    vector_modulus_a = np.linalg.norm(vector_a)
    vector_modulus_b = np.linalg.norm(vector_b)
    vector_modulus_c = np.linalg.norm(vector_c)
    list_to_judge = [vector_modulus_a, vector_modulus_b, vector_modulus_c]
    num = 0
    num_occurrences = list_to_judge.count(num)
    print(vector_abx,vector_aby,vector_acx,vector_acy,vector_bcx,vector_bcy,vector_c,vector_a,vector_modulus_b,vector_modulus_c,num_occurrences,vector_b,num_occurrences)

planar_fermat_problem()
