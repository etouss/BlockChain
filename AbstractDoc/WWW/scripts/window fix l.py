import matplotlib.pyplot as plt
import numpy as np
import math

## alpha and beta
a = 0.9999967030
b =  0.9999996156


# generating function of Catalan numbers:
def ff(x): 
	return (1-np.sqrt(1-4*x))/(2*x)


# always fork
def af(h):

	## utility of the always fork strategy as calculated in the paper
	psi = (a*b*h)*ff(a*b*b*h*(1-h))
	fi = ((a*b*h)/((1-a)*(1-b)))*(ff(b*b*h*(1-h))-a*ff(a*b*b*h*(1-h)))

	always_fork = fi/(1-psi)

	return always_fork

# factorial function
def fact(x):
	res=1

	if (x!=0):
		for i in range(1,x+1):
			res*=i;	

	return res

# choose function
def choose(x,y):
	num = fact(x)
	den = fact(y)*fact(x-y)

	res = num/den

	return res	


# Cat_gen:
def Cat_g(x,aa,bb): 
	sum = 0
	for i in range(0,bb+1):
		sum += (((x**i)*(aa+1))/(aa+i+1))*choose(aa+2*i,aa+i)
	return sum

# choose genrealized
def choose_gen(i,j):
	res = (i+1)/(i+j+1) * choose(i+2*j,i+j)

	return res


# functio

def f(aa,r,bb):
	res = 0

	for i in range(0,r-aa+1):
 		res += choose_gen(aa-1,i) * choose_gen(aa+bb-r,r-aa-i)

	return res	

# pentagon generation
def Pent_aux(aa,bb,r):

 	if (r<=aa):
 		res = choose(bb+r,r)	

 	elif (r<=bb):
 		res = choose_gen(bb-r,r)

 		for i in range(1,aa+1):
 			res += f(i,r,bb)

 	else:
 		res = choose_gen(r-bb,bb)

 		for i in range(r-bb+1,aa+1):
 			res += f(i,r,bb)

 	return res

# real Pent
def Pent(aa,bb,r): 	
	return Pent_aux(aa,bb-1,r)

# part a1 of the equation
def a1(k,j,l,h):
	x = b*b*h*(1-h)
	y = a*x

	res = (a*b*h)/((1-a)*(1-b)) * x**j * (Cat_g(x,j-1,l-j) - (a**(j+1) * Cat_g(y,j-1,l-j)))

	return res

# sum a1
def sum_a1(k,l,h):
	sum = 0
	for j in range(1,k+1):
		sum += a1(k,j,l,h)

	return sum

# part a2 of the equation
def a2(k,j,l,h):
	x = b*b*h*(1-h)
	y = a*x

	res = (a*b*h) * y**j * Cat_g(y,j-1,l-j)

	return res

# sum a2
def sum_a2(k,l,h):
	sum = 0
	for j in range(1,k+1):
		sum += a2(k,j,l,h)

	return sum	

# part b1 of the equation
def b1(k,j,l,h):
	x = b*b*h*(1-h)
	y = a*x

	res = (b**(j+k+1) * a**(j-k+1) * h**(k+1) * (1-h)**j)/((1-a)*(1-b)) * (Cat_g(x,k-1,l-k) - (a**(k+1) * Cat_g(y,k-1,l-k)))

	return res

# sum b1
def sum_b1(k,l,h):
	x = b*b*h*(1-h)
	y = a*x

	res = a*a * x**(k+1)/((1-a)*(1-b)*(1-a*b*(1-h))) * (Cat_g(x,k-1,l-k) - (a**(k+1) * Cat_g(y,k-1,l-k)))

	return res		

# part b2 of the equation
def b2(k,j,l,h):
	x = b*b*h*(1-h)
	y = a*x

	res = b**(j+k+1) * a**(j+1) * h**(k+1) * (1-h)**j * Cat_g(y,k-1,l-k)

	return res	


# sum b2
def sum_b2(k,l,h):
	x = b*b*h*(1-h)
	y = a*x

	res = a * y**(k+1)/(1-a*b*(1-h)) * Cat_g(y,k-1,l-k)

	return res		


# E_a,b
def E(aa,bb):
	res = choose(aa+2*bb,aa+bb)
	res *= (aa+1)
	res /= (aa+bb+1)

	return res



# part c of the equation
def c(k,j,l,h):
	x = b*b*h*(1-h)
	y = a*x	

	res = 0

	for r in range(0,l):
		res += (b*h)**r * Pent(j-1,l+1-j,r)
	res *= y * (a*b*(1-h))**l

	
	return res


# sum c
def sum_c(k,l,h):
	sum = 0
	for j in range(1,k+1):
		sum += c(k,j,l,h)

	return sum		


def PP(l,j,x):
	res = 0

	for r in range(0,l):
		res += (x)**r * Pent(j-1,l+1-j,r)

	return res	

# part d of the equation
def d(k,j,l,h):
	x = b*b*h*(1-h)
	y = a*x	

	res = 0

	for r in range(0,l):
		res += (b*h)**r * Pent(k-1,l+1-k,r)

	#res = a**(j+l-k) * b**(l+j-k+1) * h * (1-h)**(j+l-k) * (1-(b*h)**(l-1))/(1-b*h) * E(k-1,l-k)

	res *= y * (a*b*(1-h))**(j+l-k)

	return res


# sum d
def sum_d(k,l,h):
	x = b*b*h*(1-h)
	y = a*x

	#res = b * h * (a*b*(1-h))**(l+1) *  (1-(b*h)**(l-1)) *  E(k-1,l-k) / ((1-a*b*(1-h))*(1-b*h))
	res = y * (a*b*(1-h))**(l+1) / (1 - a*b*(1-h)) * PP(l,k,b*h)

	return res			


def util(k,l,h):
	res = a*b*h/(1-b) + sum_a1(k,l,h) + sum_b1(k,l,h)

	den = 1 - a*b*h - sum_a2(k,l,h) - sum_b2(k,l,h) - sum_c(k,l,h) - sum_d(k,l,h)

	res /= den

	return res*10

# Data for plotting

#print(sum_b2(3,5,0.5))

#print(choose(28,20))

start = 0.0001
end = 0.9999
step = 0.001

h = np.arange(start, end, step)

default = h*a*b/((1-b)*(1-a*b))*10

always_fork = af(h)

window_fork = {}
labels = []

linestyles = [':', '--', '-.', ':']
#linestyles = [None,None,None,None]
colors = [0.1,0.25,0.4,0.55,0.7]

fig, ax = plt.subplots()
#ax.set_ylim([1e12,5e12])

#for j in range(1,2):
#    for i in range(1,5):
#        if (j <= i):
#            labels.append("window: " + str(j) + ", give up: " + str(2*i))
#            ax.plot(h,util(j,2*i,h),linestyle=linestyles[i-1])

window = [1,2,3,3];
gu =     [4,5,6,6];


labels.append("$\mathbf{G}^{k = " + str(window[0]) +
                      "}_{\ell = " + str(gu[0]) + "}$")
w1g2 = util(window[0],gu[0],h)
ax.plot(h,util(window[0],gu[0],h),linestyle=linestyles[0],color=str(colors[0]))

labels.append("$\mathbf{G}^{k = " + str(window[1]) +
                      "}_{\ell = " + str(gu[1]) + "}$")
w1g4 = util(window[1],gu[1],h)
ax.plot(h,util(window[1],gu[1],h),linestyle=linestyles[1],color=str(colors[1]))

labels.append("$\mathbf{G}^{k = " + str(window[2]) +
                      "}_{\ell = " + str(gu[2]) + "}$")
w1g6 = util(window[2],gu[2],h)
ax.plot(h,util(window[2],gu[2],h),linestyle=linestyles[2],color=str(colors[2]))

labels.append("$\mathbf{G}^{k = " + str(window[3]) +
                      "}_{\ell = " + str(gu[3]) + "}$")
w1g8 = util(window[3],gu[3],h)
ax.plot(h,util(window[3],gu[3],h),linestyle=linestyles[3],color=str(colors[3]))


y = 0

for x in np.arange(0.001, 0.999, step):
	intersection = util(window[0],gu[0],x) - util(window[1],gu[1],x)
	if intersection < 0.00000005:
		y = x		
		break

section = np.arange(start,y,step)
#plt.fill_between(section,util(1,2,section),color='Grey').set_hatch('///')

plt.fill_between(section,util(window[0],gu[0],section),facecolor='1',hatch='\\\\\\',edgecolor='0.8')
#plt.fill_between(section,util(1,2,section),facecolor='0.85',edgecolor=str(colors[0]))


y1 = 0

for x in np.arange(0.001, 0.999, step):
	intersection = util(window[1],gu[1],x) - util(window[2],gu[2],x)
	if intersection < 0.00000005:
		y1 = x		
		break

section = np.arange(y-step,y1,step)
plt.fill_between(section,util(window[1],gu[1],section),facecolor='1',hatch='....',edgecolor='0.8')
#plt.fill_between(section,util(1,6,section),facecolor='0.25',edgecolor=str(colors[1]))

y2 = 0

for x in np.arange(0.001, 0.999, step):
	intersection = util(window[2],gu[2],x) - util(window[3],gu[3],x)
	if intersection < 0.00000005:
		y2 = x		
		break

section = np.arange(y1-step,y2,step)
plt.fill_between(section,util(window[2],gu[2],section),facecolor='1',hatch='----',edgecolor='0.8')
#plt.fill_between(section,util(1,6,section),facecolor='0.25',edgecolor=str(colors[1]))


section = np.arange(y2-step,end,step)
plt.fill_between(section,util(window[3],gu[3],section),facecolor='1',hatch='////',edgecolor='0.8')
#plt.fill_between(section,util(1,10,section),facecolor='0.55',edgecolor=str(colors[2]))

#ax.plot(y,util(1,2,y),'ro')

ax.plot(h,default,color='black')

labels.append('default')

#ax.plot(h,always_fork,color='blue')

#labels.append('always fork')


#ax.set(xlabel='hash power (h)', ylabel='utility', title='Utility for a=%r, b=%r' % (a,b) )

ax.set(xlabel='hash power (h)', ylabel='utility')

#ax.grid()
#ax.yaxis.grid() # horizontal lines

major_ticks = np.arange(start, end, 0.1)

#ax.set_xticks([0.5])
ax.set_xticks([0.5])
ax.set_yticks([])

for item in ax.get_xticklabels():
	item.set_fontsize(18)

ax.xaxis.label.set_fontsize(18)
ax.yaxis.label.set_fontsize(18)

ax.xaxis.grid() # vertical lines
ax.legend(labels,fontsize = 'xx-large')

fig.savefig("window_fork.png")
plt.show()



