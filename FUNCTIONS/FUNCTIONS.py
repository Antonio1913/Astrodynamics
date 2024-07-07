import math

#arccot finds the inverse cotangent
def arccot(x):
    return math.atan(1/x)

#Sign outputs the sign of the value, num, inputted
def sign(num):
    return -1 if num < 0 else 1

#Should Input a [N x 3] vector
#Error message will output if formatting is not correct
def norm(vec)
    vec