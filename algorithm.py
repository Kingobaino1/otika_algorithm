import numpy as np
import pandas as pd  
from numpy.linalg import norm
import matplotlib.pyplot as plt
# %matplotlib inline
import time
import seaborn as sns

# defining of function parameters
def theta_n(n):
    return (n)/(1.45*(n+1))
def t2heta_n(n):
    return ((n**2)-4)/(7.5*(n**2))
def gamma_n(n):
    return ((3*n+3)**(1/float(n+1)))/32
def newgamma_n(n):
    return ((3 * n)**(1 / float(n))) / 32
def g2amma_n(n):
    return (n+1)/(4*(n+2))
def newg2amma_n(n):
    return (n)/(4*(n+1))
def nom(x):
    return norm(x,2)
def A(x):
    if nom(x) > 3:
      value = x
    else:
      value = (4 - nom(x))*x
    return value
def PC(x):
    if nom(x) > 3:
      value = (3*x)/nom(x)
    else:
      value = x
    return value
x_0 = np.array([[1],[1/2],[1/4],[1/8],[1/16]])
x_1 = np.array([[4/5],[16/25],[64/125],[256/625],[1024/3125]])

st = time.time()
# initialization
count = 1
TOL = []
iterations = []
while (nom(x_1-x_0))**2 > 10**(-7):
    w_n = x_1 + theta_n(count)*(x_1-x_0)
    d_n = w_n -((gamma_n(count)+newgamma_n(count))*A(x_1) - newgamma_n(count)*A(x_0))
    x_0 = x_1
    x_1 = PC(d_n)
    TOL.append(nom(x_1-x_0)**2)
    iterations.append(count)
    count += 1
    print(x_1, "after", count, "iterations")
    #print(x_1, "after", count, "iterations") #Uncomment to see how it iterated
else:
  print("Iteration has ended with x=", x_1, "after", count - 1, "iterations")

et = time.time()
TOTAL_TIME = et - st
print('\nExecution time for iteration 1:', TOTAL_TIME, 'seconds\n')

# initialization
st1 = time.time()
count1 = 1
TOL1 = []
iterations1 = []
y_0 = np.array([[1],[1/2],[1/4],[1/8],[1/16]])
y_1 = np.array([[4/5],[16/25],[64/125],[256/625],[1024/3125]])
while (nom(y_1-y_0))**2 > 10**(-7):
# while (count1 <= 129):
    u_n = y_1 + t2heta_n(count1)*(y_1-y_0)
    v_n = u_n -((g2amma_n(count1)+newg2amma_n(count1))*A(y_1) - newg2amma_n(count1)*A(y_0))
    y_0 = y_1
    y_1 = PC(v_n)
    TOL1.append(nom(y_1-y_0)**2)
    iterations1.append(count1)
    count1 += 1
    print(count1, 'count')
    print(y_1, "after", count1, "iterations")
    #print(x_1, "after", count, "iterations") #Uncomment to see how it iterated
else:
  print("Iteration has ended with x=", y_1, "after", count1 - 1, "iterations")
et1 = time.time()
TOTAL_TIME1 = et1 - st1
print('\nExecution time for iteration 2:', TOTAL_TIME1, 'seconds\n')

# initialization
st2 = time.time()
count2 = 1
TOL2 = []
iterations2 = []
z_0 = np.array([[1],[1/2],[1/4],[1/8],[1/16]])
z_1 = np.array([[4 / 5], [16 / 25], [64 / 125], [256 / 625], [1024 / 3125]])
while (nom(z_1-z_0))**2 > 10**(-7):
    q_n = z_1 + theta_n(count2)*(z_1-z_0)
    h_n = q_n -((g2amma_n(count2)+newg2amma_n(count2))*A(z_1) - newg2amma_n(count2)*A(z_0))
    z_0 = z_1
    z_1 = PC(h_n)
    TOL2.append(nom(z_1-z_0)**2)
    iterations2.append(count2)
    count2 += 1
    # print(z_1, "after", count2, "iterations")
    #print(x_1, "after", count, "iterations") #Uncomment to see how it iterated
else:
  print("Iteration has ended with x=", z_1, "after", count2 - 1, "iterations")
et2 = time.time()
TOTAL_TIME2 = et2 - st2
print('\nExecution time for iteration 3:', TOTAL_TIME2, 'seconds\n')

# initialization
st3 = time.time()
count3 = 1
TOL3 = []
iterations3 = []
za_0 = np.array([[1],[1/2],[1/4],[1/8],[1/16]])
za_1 = np.array([[4/5],[16/25],[64/125],[256/625],[1024/3125]])
while (nom(za_1-za_0))**2 > 10**(-7):
#while (count3 <= 129):
    qa_n = za_1 + t2heta_n(count3)*(za_1-za_0)
    ha_n = qa_n -((gamma_n(count3)+newgamma_n(count3))*A(za_1) - newgamma_n(count3)*A(za_0))
    za_0 = za_1
    za_1 = PC(ha_n)
    TOL3.append(nom(za_1-za_0)**2)
    iterations3.append(count3)
    count3 += 1
    # print(za_1, "after", count3, "iterations")
    #print(x_1, "after", count, "iterations") #Uncomment to see how it iterated
else:
  print("Iteration has ended with x=", za_1, "after", count3 - 1, "iterations")
et3 = time.time()
TOTAL_TIME3 = et3 - st3
print('\nExecution time for iteration 4:', TOTAL_TIME3, 'seconds\n')


plt.figure(figsize=(7, 7), dpi=80)
plt.plot(iterations,TOL, "bx-.", color="r", label="Case 1")
plt.plot(iterations1,TOL1, "bx-.", color="yellow", label="Case 2")
plt.plot(iterations2,TOL2, "bx-.", color="blue", label="Case 3")
plt.plot(iterations3,TOL3, "bx-.", color="green", label="Case 4")
plt.ylim(0, .001)
plt.xlim(1, 80)
plt.xlabel("No of iterations (n)")
plt.ylabel("TOL_n")
plt.legend()

plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)
# plt.legend((arr1, arr2), ('arr1', 'arr2'))
plt.show()# -*- coding: utf-8 -*-