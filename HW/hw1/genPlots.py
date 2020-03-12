import matplotlib.pyplot as plt
import numpy as np
from scipy.special import comb
'''
The code below was found and stolen from
https://www.geeksforgeeks.org/program-for-nth-fibonacci-number/
Credit to ChitraNayal
'''

# Helper function that multiplies
# 2 matrices F and M of size 2*2,
# and puts the multiplication
# result back to F[][]

# Helper function that calculates
# F[][] raise to the power n and
# puts the result in F[][]
# Note that this function is
# designed only for fib() and
# won't work as general
# power function
def fib(n):
    F = [[1, 1],
         [1, 0]]
    if (n == 0):
        return 0
    power(F, n - 1)

    return F[0][0]

def multiply(F, M):

    x = (F[0][0] * M[0][0] +
         F[0][1] * M[1][0])
    y = (F[0][0] * M[0][1] +
         F[0][1] * M[1][1])
    z = (F[1][0] * M[0][0] +
         F[1][1] * M[1][0])
    w = (F[1][0] * M[0][1] +
         F[1][1] * M[1][1])

    F[0][0] = x
    F[0][1] = y
    F[1][0] = z
    F[1][1] = w

def power(F, n):

    M = [[1, 1],
         [1, 0]]

    # n - 1 times multiply the
    # matrix to {{1,0},{0,1}}
    for i in range(2, n + 1):
        multiply(F, M)

def calcProb(i,n):
    reg = fib(n+2)
    if(i==1 or i == n):
        return fib(n)/reg
    else:
        LS = fib(i)
        RS = fib(n-i+1)
    return (LS*RS/reg)

# This code is contributed
# by ChitraNayal

def calcPossibleSize(gens,k,sign):
    '''
    elif(gens == 1):
        counter = 1
        tot=0
        while(counter<=k):
            tot+=comb(k,counter)
            counter+=1
        return tot+1
    elif(gens == -1):
        return 0
    '''
    if(gens == 0):
        return 1
    else:
        if(sign==0):
            return (calcPossibleSize(gens-1,k,0)+calcPossibleSize(gens-1,k,1))**k
        else:
            return calcPossibleSize(gens-1,k,0)**k

def CPT(gens,k):
    if(gens==0):
        top = calcPossibleSize(gens,k,1)
        bot = top + calcPossibleSize(gens,k,0)
        return top/bot
    else:
        top = 1
        botp1 = 1
        botp2 = 1/(1-CPT(gens-1,k))
        botp3 = botp2**k
        bot = botp1+botp3
        return top/bot

def CPTalt(gens,k):
    top = calcPossibleSize(gens,k,1)
    bot = calcPossibleSize(gens,k,0)
    return top/(top+bot)


'''
The following code generates plots for the probability according to my formula P(L_x_i) = 1/Z(L_n) * (Z(L_(i-1)) + Z(L_(n-i)))
'''
n= 11
probs = []
if __name__ == "__main__":
    i = 1
    debug=True
    if(debug):
        ls = np.arange(0,51,1)
        ks = [1,2,3,4,10]
        for k in ks:
            vals=[]
            for l in ls:
                vals.append(CPT(l,k))
            plt.clf()
            plt.plot(ls,vals)
            plt.xlabel("l")
            #ax.set_xticklabels(np.arange(1,12,1).astype(str))
            plt.ylabel("Probability of root node = 1")
            plt.title('P(root node = 1) over l generations with '+str(k) + 'children')
            plt.savefig('p2_4_'+str(k))
            print(vals[-1])
    else:
        while(i<=n):
            probs.append(calcProb(i,n))
            i+=1
        fig, ax = plt.subplots()
        plt.plot(np.arange(1,12,1),probs)
        plt.xlabel("possible x_i")
        #ax.set_xticklabels(np.arange(1,12,1).astype(str))
        plt.ylabel("Probability of L_x_i")
        plt.title('Probability over different nodes')
        #plt.show()
        plt.savefig('p2_3')
