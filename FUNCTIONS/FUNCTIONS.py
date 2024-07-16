import math
import numpy as np

#arccot finds the inverse cotangent
def arccot(x):
    return math.atan(1/x)

#Sign outputs the sign of the value, num, inputted
def sign(num):
    return -1 if num < 0 else 1

# Rotation Matrices

#Rot 1 finds the rotation about the X axis
def Rot1(alpha):
    matrix = np.array([1, 0, 0],
                      [0, np.cos(alpha), np.sin(alpha)],
                      [0, -np.sin(alpha), np.cos(alpha)])
    return matrix

#Rot 2 finds the rotation about the Y axis
def Rot2(alpha):
    matrix = np.array([np.cos(alpha), 0, -np.sin(alpha)],
                                    [0, 1, 0],
                            [np.sin(alpha), 0, np.cos(alpha)])
    return matrix

#Rot 3 finds the rotation about the Z axis
def Rot3(alpha):
    matrix = np.array([np.cos(alpha), np.sin(alpha), 0],
                      [-np.sin(alpha), np.cos(alpha), 0],
                                        [0, 0, 1])
    return matrix