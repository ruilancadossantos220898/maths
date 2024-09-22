import numpy as np 

print("value of a, b and c for: ax^2 + bx + c")


a=input('a=?\n') 
b=input('b=?\n') 
c=input('c=?\n') 

a=1
b=-2
c=1


a = float(a)
b = float(b)
c = float(c)


if b**2-4*a*c > 0:
    print("Have 2 different real roots")

if b**2-4*a*c == 0:
    print("Have 2 real roots wich are equals")

if b**2-4*a*c < 0:
    print("Have 2 imaginary roots")



x1 = (-b+np.sqrt(b**2-4*a*c))/(2*a)
x2 = (-b-np.sqrt(b**2-4*a*c))/(2*a)

print("(x-" + str(x1) + ")(x-" + str(x2) + ")")


print("x1=" + str(x1))
print("x2=" + str(x2))

