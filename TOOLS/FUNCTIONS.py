import math
import numpy as np
from SOLN_TO_KEPLER.KEPLER import KEPLER

#arccot finds the inverse cotangent
def arccot(x):
    return np.arctan(1/x)

#Sign outputs the sign of the value, num, inputted
def sign(num):
    return -1 if num < 0 else 1

# Rotation Matrices

#Rot 1 finds the rotation about the X axis
def Rot1(alpha):
    matrix = np.array([[1, 0, 0],
                      [0, np.cos(alpha), np.sin(alpha)],
                      [0, -np.sin(alpha), np.cos(alpha)]])
    return matrix

#Rot 2 finds the rotation about the Y axis
def Rot2(alpha):
    matrix = np.array([[np.cos(alpha), 0, -np.sin(alpha)],
                                    [0, 1, 0],
                            [np.sin(alpha), 0, np.cos(alpha)]])
    return matrix

#Rot 3 finds the rotation about the Z axis
def Rot3(alpha):
    matrix = np.array([[np.cos(alpha), np.sin(alpha), 0],
                      [-np.sin(alpha), np.cos(alpha), 0],
                                        [0, 0, 1]])
    return matrix


class DoubleRangeValue:
    def __init__(self, value, range1, range2):
        self.value = value
        self.range1 = range1
        self.range2 = range2

    def is_within_ranges(self):
        return (self.range1[0] <= self.value <= self.range1[1]) or (self.range2[0] <= self.value <= self.range2[1])

    def __repr__(self):
        return f"DoubleRangeValue(value={self.value}, range1={self.range1}, range2={self.range2})"

def withinrange(value, start, end):
    return start <= value <= end








