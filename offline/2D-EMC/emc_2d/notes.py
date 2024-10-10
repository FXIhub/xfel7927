import numpy as np
import matplotlib.pyplot as plt


# I would like a way to reduce iterations when trying to solve:
# sum_i a_i / (x + b_i) = c
# since each iteration requires a sum over i and evaluation of a_i and b_i 
# which is expensive in terms of non-local memory and floating point operations
# 
# So fit a taylor series to the target function then use newton on the coefficients
# we then have a single pass algorithm

# this will strugle when the solution is near zero (c ~ sum_i a_i / b_i)
# and in fact the taylor series seems to be divergent there (from tests below)
# 
# perhaps a taylor series expansion could help?
# g(x') = sum_i a_i / (x' + b_i)
# g(x-x') approx. sum_i a_i / (x' + b_i) 
#      + (x-x')   sum_i -  a_i / (x' + b_i)^2 
#      + (x-x')^2 sum_i +2 a_i / (x' + b_i)^3  / 2
#      + (x-x')^3 sum_i -6 a_i / (x' + b_i)^4  / 2*3
#  
#        = sum_n (x-x')^n / n! sum_i a_i (-1)^n n! / (x' + b_i)^(n+1)
#        = sum_n (x-x')^n [ (-1)^n sum_i a_i / (x' + b_i)^(n+1)]
#        = sum_n (x-x')^n c_n 
#
# what about a different function?
# g(x) = a' / (x + b') + sum_n=1 c_n (x-x')^n
#
# b'    = min_i(b) - this ensures the asymtotes of a' / (x + b') and g match
# g(x') = a' / (x' + b') --> a' = g(x') (x' + b')
#
# g(x')  = sum_i   a_i / (x' + b_i)   =  a' / (x' + b')  
# g(x')1 = sum_i  -a_i / (x' + b_i)^2 = -a' / (x + b')^2 + c_1
# g(x')2 = sum_i  2a_i / (x' + b_i)^3 = 2a' / (x + b')^3 + 2 c_2
# g(x')n = sum_i n! (-1)^n a_i / (x' + b_i)^(n+1) = (-1)^n n! a' / (x + b')^(n+1) + n! c_n

# a'  = g(x') (x' + b')
# c_1 = sum_i -a_i / (x' + b_i)^2 - (-a' / (x+b')^2)
# c_n = (-1)^n [ sum_i a_i / (x' + b_i)^(n+1) - a' / (x+b')^(n+1) ]


N = 1000
a = np.random.random(N)
b = 10*np.random.random(N)
c = 100
xmax = np.sum(a) / c

g_calc = lambda x: np.sum(a / (x + b))

x0 = xmax

x = np.linspace(0, xmax, 1000)
g = np.array([g_calc(xx) for xx in x])



# modified Taylor
"""
#b0 = np.min(b)
#a0 = g_calc(x0) * (x0 + b0)
a0 = np.sum(a)
b0 = a0 / np.sum(a/b) 

order = 10
coefs = np.zeros((order,))
for n in range(0, order):
    coefs[n] = (-1)**n * (np.sum( a / (x0 + b)**(n+1) ) - a0 / (x0 + b0)**(n+1))


g_taylor = a0 / (x + b0)
for n in range(order):
    g_taylor += (x-x0)**n * coefs[n]
"""


# Taylor

# fit taylor expansion
#g_taylor = g_calc(x0) - (x-x0) * np.sum( a / (x0 + b)**2 ) + (x-x0)**2 * np.sum( a / (x0 + b)**3 )

order = 10
coefs = np.zeros((order,))
for n in range(order):
    coefs[n] = (-1)**n * np.sum( a / (x0 + b)**(n+1) )

g_taylor = np.zeros((N,))
for n in range(order):
    g_taylor += (x-x0)**n * coefs[n]

fig, ax = plt.subplots()

ax.plot(x, g, c='k')
ax.plot(x, g_taylor, c='b')
ax.hlines(c, 0, xmax, color='k', linestyle='--')
#ax.set_yscale('log')
ax.set_ylim([g.min(), g.max()])
plt.show()
