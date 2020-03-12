import numpy as np
import random
import matplotlib.pyplot as plt
import pdb
def firstPart():
    epsilon = [0.2,0.4,0.6,0.8,1.0]
    theta=0.05
    l=10
    allas = np.arange(-1,1,0.01)
    for e in epsilon:
        part1 =2*l*l*e*(allas**2)+theta*l*l*allas
        part2 = ((1+allas)/2)*np.log((1+allas)/2)
        part3 = ((1-allas)/2)*np.log((1-allas)/2)
        part4 = part1-l*l*(part2+part3)
        print(allas)
        print(part4)
        plt.plot(allas,part4,label=str(e))
    plt.legend()
    plt.title("free energy Fmv(bv) of various epsilon and a")
    plt.show()
#firstPart()
def dervFHelper(a):
    opa = 1+a
    oma = 1-a
    part1 = 1/opa+np.log(opa/2)+a/opa+1/oma+np.log(oma/2)-a/oma
    return 0.5*part1

def altApproach():
    epsilon = np.arange(0,1,0.05)
    theta=0.05
    l=10
    allas = np.arange(-0.99,1,0.01)
    maxBs=[]
    maxFs=[]
    for e in epsilon:
        part1 =2*l*l*e*(allas**2)+theta*l*l*allas
        part2 = ((1+allas)/2)*np.log((1+allas)/2)
        part3 = ((1-allas)/2)*np.log((1-allas)/2)
        part4 = part1-l*l*(part2+part3)
        #pdb.set_trace()
        maxloc = (np.argmax(np.array(part4)))
        maxA = allas[maxloc]
        maxFs.append(part4[maxloc])
        maxBs.append((maxA+1)/2)
    pdb.set_trace()
    plt.plot(epsilon,maxBs,label=str(e))
    plt.xlabel("theta")
    plt.ylabel("b*")
    plt.title("optimal b* over theta")
    plt.show()
    plt.plot(epsilon,maxFs,label=str(e))
    plt.xlabel("theta")
    plt.ylabel("F(b*)")
    plt.title("optimal F(b*) over theta")
    plt.show()
def secondPart():
    theta=0.05
    l=10
    epsilon= np.arange(0,1,0.01)
    #allas = np.arange(-2,1,0.01)
    allAs = np.arange(-0.99,1,0.001)
    minLoc = []
    bestBs=[]
    for e in epsilon:
        derivs=[]
        #bs=[]
        minDeriv = 1000
        minB=None
        for a in allAs:
            b = (a+1)/2
            part1 = 4*l*l*e*b + theta*l*l
            part2 = l*l/2
            part3=np.log(1-b)+np.log(1+b)+2-np.log(4)
            altRes = part1-l*l*dervFHelper(b)
            #pdb.set_trace()
            '''part1 = 4*l*l*e*b + theta*l*l
            part2 = l*l/2
            part3=np.log(1-b)+np.log(1+b)+2-np.log(4)'''
            result = part1-part2*part3
            if(altRes<minDeriv):
                minDeriv=altRes
                minB = b
            derivs.append(abs(altRes))
        pdb.set_trace()
        bestBs.append(minB)
        minLoc.append((allAs[np.argmin(np.array(derivs))]+1)/2)
    minLoc = np.array(minLoc)
    actualBs = minLoc#(minLoc+1)/2
    pdb.set_trace()
    plt.title("values of bv(1) with derivative 0 over epsilon")
    plt.plot(epsilon,bestBs)
    plt.show()
    return actualBs
def thirdPart(abs):
    epsilon= np.arange(0,1,0.01)
    theta=0.05
    l=10
    b = 0.5
    results=[]
    counter=0
    for e in epsilon:
        b = abs[counter]
        a = 2*b-1
        part1 =2*l*l*e*(a**2)+theta*l*l*a
        part2 = ((1+a)/2)*np.log((1+a)/2)
        part3 = ((1-a)/2)*np.log((1-a)/2)
        part4 = part1-l*l*(part2+part3)
        #   print(part1)
        results.append(part4)
        counter+=1
    plt.plot(epsilon,results)
    plt.title("free energy Fmv(b*) of various epsilon")
    plt.show()
#abs=secondPart()
#thirdPart(abs)
#firstPart()
altApproach()
