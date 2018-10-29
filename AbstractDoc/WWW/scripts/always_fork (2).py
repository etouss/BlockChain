import matplotlib.pyplot as plt
import numpy as np
import math

# generating function of Catalan numbers:
def f(x): 
	return (1-np.sqrt(1-4*x))/(2*x)

# Data for plotting

## alpha and beta
a = 0.9999967030
b = 0.9999996156

def af(h):
	psi = (a*b*h)*f(a*b*b*h*(1-h))
	fi = ((a*b*h)/((1-a)*(1-b)))*(f(b*b*h*(1-h))-a*f(a*b*b*h*(1-h)))

	return (1-b)*fi/(1-psi)

def fo(h):
	## fork once strategy (genesis fork)
	k1 = (a*b*h)/((1-a)*(1-b))
	k2 = (a*a*b*h*(a*b+h*b-h*a*b-1))/((1-a)*(1-b)*(1-a*b))

	return (1-b)*(k1*f(b*b*h*(1-h)) + k2*f(a*b*b*h*(1-h)))


def default_comp(h):
	return (1-b)*h*a*b/((1-b)*(1-a*b))

## the range and granularity to consider
h = np.arange(0.001, 0.999, 0.001)

always_fork = af(h)

fork_once = fo(h)

default = default_comp(h)

## plotting
fig, ax = plt.subplots()

# plot the utlities
ax.plot(h,always_fork,color='blue')
ax.plot(h,fork_once,color='red')
ax.plot(h,default,color='pink')

ax.set(xlabel='hash power (h)', ylabel='utility',
       title='Utility for a=%r, b=%r' % (a,b) )
ax.grid()

ax.legend(['always fork','fork once','default'])

fig.savefig("test_AF.png")
plt.show()
