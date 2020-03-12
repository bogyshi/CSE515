#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt


stockData = np.genfromtxt('sp500.csv', delimiter=',')

print(stockData)
print(len(stockData))
q=0.9
p_up_good = q
p_up_bad = 1 - p_up_good
upPs = np.array([p_up_good,p_up_bad])
p_down_good=1-q
p_down_bad = 1-p_down_good
downPs = np.array([p_down_good,p_down_bad])
initPs2 = np.array([[p_up_good,p_up_bad],[p_down_good,p_down_bad]])
updownMatrix = np.array([[0.8,0.2],[0.2,0.8]]) # being in the good state is top row, bad state is bottom, 1st col represents good prob, second represent bas prob

pStay = 0.8
pChange = 0.2

def sumProdRecursive(data,index):
    ps = np.array([0,0])
    if(index == len(data)):
        return np.array([0.2,0.8])
    else:
        message = sumProd(data,index-1)



def sumProdIter(data):
    totalProbOfEvents=1
    lenData = len(data)
    i = 0
    measurements = []
    initPs = np.array([0.2,0.8])
    message = initPs
    while(i<lenData):
        pGood = message[0]
        pBad = message[1]
        updown = data[i]
        if(updown == 1):
            diagMatrix = np.diag(initPs2[:,0]) # up on top, down on bottom, good state 0th col, bad state 1th col
            top = p_up_good*pGood
            bottom = p_up_good*pGood + p_up_bad*pBad
            #message[0] = (pGood*pStay*p_up_good + pBad*pChange*p_up_good)
            #message[1] = (pBad*pStay*p_up_bad + pGood*pChange*p_up_bad)
        else:
            diagMatrix = np.diag(initPs2[:,1])
            top = p_down_good*pGood
            bottom = p_down_good*pGood + p_down_bad*pBad
            #message[0] = (pGood*pStay*p_down_good + pBad*pChange*p_down_good)
            #message[1] = (pBad*pStay*p_down_bad + pGood*pChange*p_down_bad)
        ptop = top/bottom
        totalProbOfEvents = totalProbOfEvents*np.sum(message.dot(updownMatrix).dot(diagMatrix))
        message = message.dot(updownMatrix).dot(diagMatrix)/np.sum(message.dot(updownMatrix).dot(diagMatrix))
        print(message)
        print(ptop)

        #message[0] = (pGood*pStay + pBad*pChange)*ptop
        #message[1] = (pBad*pStay + pGood*pChange)*ptop
        measurements.append(ptop)
        i+=1
    return measurements,totalProbOfEvents
res,totalProbOfEvents2= (sumProdIter(stockData))
print(totalProbOfEvents2)

plt.plot(res)
plt.title("P(x='good'|y) over all weeks, q = " + str(q))
plt.xlabel("week")
plt.ylabel("P(x='good'|y)")
plt.show()
