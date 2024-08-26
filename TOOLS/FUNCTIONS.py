import numpy as np
import scipy as sc
import pandas as pd

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


#INPUTS
# x                     - [1XN] function of
# desired_degree        - [1X1]desired_degree

# Pbar_lm               - Normalized Associated Legendre Function
def LegendreFunction(x, desired_degree, HarmonicValues):

    # Saving x calculation
    xval = 1 - (x**2)

    # Pre-Allocating Arrays for Normalized Plm
    Pbarlm = np.zeros([desired_degree + 2, desired_degree-1])

    #Beginning for Loops to Calculate Recursions
    for degrees in range(2, 1, desired_degree+1):

        # Degree correction to store values in correct column
        store = degrees - 2

        # Calculating Pl,l (degrees - 2 stores values into first column, -2 stores value into second to last column)
        Plm[-2, store] = (sc.special.factorial((2 * degrees) - 1) / (2**(degrees - 1) * sc.special.factorial(degrees - 1))) * (xval**(degrees / 2))

        # Calculating Pl,l-1 (degrees - 2 stores values into first column, -3 stores value into third to last column
        Plm[-3, store] = (x / xval**(1/2)) * Plm[degrees, -1]

        # Normalization
        # Pre-allocating delta variable with 0 and 1 values
        deltam0 = np.zeros([degrees, 1])
        deltam0[1] = 1

        # (degrees - 2 starts for loop at Pl,l-2 since Pl,l and Pl,l-1 have already been calculated)
        for order in range(degrees - 2, -1, -1):
            calc1 = 1 / ((degrees - order) * (degrees + order + 1))
            # Calculating Pl,m
            Plm[order, store] = calc1 * (((2 * (order + 1)) * (x / xval**(1/2)) * Plm[order+1, store]) - Plm[order+2, store])

        orders = np.arange(degrees+1)

        #Normalization
        Pbarlm = ((2 - deltam0) * ((2 * degrees) + 1) * ((sc.special.factorial(degrees - order)) / (sc.special.factorial(degrees + order))))**(1/2) * Plm






