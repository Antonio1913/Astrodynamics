import numpy as np


# arccot finds the inverse cotangent
def arccot(x):
    return np.arctan(1/x)


# Sign outputs the sign of the value, num, inputted
def sign(num):
    return -1 if num < 0 else 1


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


