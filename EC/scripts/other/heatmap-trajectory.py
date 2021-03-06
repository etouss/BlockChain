import matplotlib.pyplot as plt
import numpy as np

# alpha and beta
a = 0.9999967030
b = 0.9999996156


# generating function of Catalan numbers
# Written as 2 / (sqrt(1 - 4*x) + 1) to deal with continuity at 0
# (We don't care about x -> 1/2 because of alpha, beta)
def ff(x):

    return 2 / (np.sqrt(1 - 4 * x) + 1)


# always fork
def af(h):

    # utility of the always fork strategy as calculated in the paper
    psi = (a * b * h) * ff(a * b * b * h * (1 - h))
    fi = ((a * b * h) / ((1 - a) * (1 - b))) * \
        (ff(b * b * h * (1 - h)) - a * ff(a * b * b * h * (1 - h)))

    always_fork = fi / (1 - psi)

    return always_fork


def default(h):
    return h * a * b / ((1 - b) * (1 - a * b))


# factorial function
def fact(x):
    res = 1

    if (x != 0):
        for i in range(1, x + 1):
            res *= i

    return res


# choose function
def choose(x, y):
    num = fact(x)
    den = fact(y) * fact(x - y)

    res = num / den

    return res


# Cat_gen:
def Cat_g(x, aa, bb):
    sum = 0
    for i in range(0, bb + 1):
        sum += (((x ** i) * (aa + 1)) / (aa + i + 1)) * \
            choose(aa + 2 * i, aa + i)
    return sum


# choose genrealized
def choose_gen(i, j):
    if i < 0 or j < 0:
        return 0
    res = (i + 1) / (i + j + 1) * choose(i + 2 * j, i + j)
    return res


# pentagon sum, inner term

def f(aa, r, bb):
    res = 0

    for i in range(0, r - aa + 1):
        res += choose_gen(aa - 1, i) * choose_gen(aa + bb - r, r - aa - i)
    return res


# pentagon generation
def Pent_aux(aa, bb, r):

    if (r <= aa):
        res = choose(bb + r, r)

    elif (r <= bb):
        res = choose_gen(bb - r, r)
        for i in range(1, aa + 1):
            res += f(i, r, bb)
    else:
        res = choose_gen(r - bb, bb)
        for i in range(r - bb + 1, aa + 1):
            res += f(i, r, bb)

    return res


# real Pent
def Pent(aa, bb, r):
    return Pent_aux(aa, bb - 1, r)


# part a1 of the equation
def a1(k, j, l, h):
    x = b * b * h * (1 - h)
    y = a * x
    res = (a * b * h) / ((1 - a) * (1 - b)) * x**j * \
        (Cat_g(x, j - 1, l - j) - (a ** (j + 1) * Cat_g(y, j - 1, l - j)))
    return res


# sum a1
def sum_a1(k, l, h):
    sum = 0
    for j in range(1, k + 1):
        sum += a1(k, j, l, h)

        return sum


# part a2 of the equation
def a2(k, j, l, h):
    x = b * b * h * (1 - h)
    y = a * x

    res = (a * b * h) * y ** j * Cat_g(y, j - 1, l - j)

    return res


# sum a2
def sum_a2(k, l, h):
    sum = 0
    for j in range(1, k + 1):
        sum += a2(k, j, l, h)

        return sum


# part b1 of the equation
def b1(k, j, l, h):
    x = b * b * h * (1 - h)
    y = a * x

    res = (b ** (j + k + 1) * a ** (j - k + 1) *
           h ** (k + 1) * (1 - h) ** j) / ((1 - a) * (1 - b))
    res *= (Cat_g(x, k - 1, l - k) - (a ** (k + 1) * Cat_g(y, k - 1, l - k)))

    return res


# sum b1
def sum_b1(k, l, h):
    x = b * b * h * (1 - h)
    y = a * x

    res = a * a * x ** (k + 1) / ((1 - a) * (1 - b) * (1 - a * b * (1 - h))) \
        * (Cat_g(x, k - 1, l - k) - (a ** (k + 1) * Cat_g(y, k - 1, l - k)))

    return res


# part b2 of the equation
def b2(k, j, l, h):
    x = b * b * h * (1 - h)
    y = a * x

    res = b ** (j + k + 1) * a ** (j + 1) * h ** (k + 1) * (1 - h) ** j \
        * Cat_g(y, k - 1, l - k)

    return res


# sum b2
def sum_b2(k, l, h):
    x = b * b * h * (1 - h)
    y = a * x
    res = a * y ** (k + 1) / (1 - a * b * (1 - h)) * Cat_g(y, k - 1, l - k)

    return res


# E_a,b - Amount of trapezoidal paths
def E(aa, bb):
    res = choose(aa + 2 * bb, aa + bb)
    res *= (aa + 1)
    res /= (aa + bb + 1)

    return res


# part c of the equation
def c(k, j, l, h):
    x = b * b * h * (1 - h)
    y = a * x
    res = 0

    for r in range(0, l):
        res += (b * h) ** r * Pent(j - 1, l + 1 - j, r)

    res *= y * (a * b * (1 - h)) ** l

    return res


# sum c
def sum_c(k, l, h):
    sum = 0
    for j in range(1, k + 1):
        sum += c(k, j, l, h)

    return sum


def PP(l, j, x):
    res = 0
    for r in range(0, l):
        res += x ** r * Pent(j - 1, l + 1 - j, r)

    return res


# part d of the equation
def d(k, j, l, h):
    x = b * b * h * (1 - h)
    y = a * x

    res = 0

    for r in range(0, l):
        res += (b * h) ** r * Pent(k - 1, l + 1 - k, r)

    # res = a**(j+l-k) * b**(l+j-k+1) * h * (1-h)**(j+l-k) * \
    # (1-(b*h)**(l-1))/(1-b*h) * E(k-1,l-k)

    res *= y * (a * b * (1 - h)) ** (j + l - k)
    return res


# sum d
def sum_d(k, l, h):
    x = b * b * h * (1 - h)
    y = a * x

    # res = b * h * (a*b*(1-h))**(l+1) *  (1-(b*h)**(l-1)) * \
    # E(k-1,l-k) / ((1-a*b*(1-h))*(1-b*h))
    res = y * (a * b * (1 - h)) ** (l + 1) / (1 - a * b * (1 - h)) \
        * PP(l, k, b * h)

    return res


def util(k, l, h):
    res = a * b * h / (1 - b) + sum_a1(k, l, h) + sum_b1(k, l, h)

    den = 1 - a * b * h - sum_a2(k, l, h) - sum_b2(k, l, h) \
        - sum_c(k, l, h) - sum_d(k, l, h)

    res /= den

    return res


def max_util(h):

    max_k = 0
    max_l = 0
    max_value = 0
    for i in range(1, 21):
        for j in range(i+1, 21):
            if util(i, j, h) > max_value:
                max_value = util(i, j, h)
                max_k = i
                max_l = j

    return max_k, max_l


for x in np.arange(0, 1, .1):
    print(max_util(x))

h = 0.49

utill = [[util(x, y, h) for x in range(1, 21)] for y in range(1, 21)]

fig, ax = plt.subplots()
im = ax.imshow(utill)
ax.invert_yaxis()

ax.set_xticks(np.arange(1,21,1))
ax.set_yticks(np.arange(1,21,1))

window_fork = {}
labels = []
fig.savefig("window_fork.png")
plt.show()

