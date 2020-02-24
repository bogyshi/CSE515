import numpy as np
import matplotlib.pyplot as plt

l = 3
r = 4

epsilon1 = 0.05
epsilon2 = 0.1

etest = epsilon2


def fofy(prevw):
    part1 = 1-((1-2*prevw)**(r-1))
    part2 = part1/2
    part3 = part2**(l-1)
    part4 = (1 - part2**(l-1) - (1 - part2)**(l-1))
    part5= part4*etest
    neww = part3+part5
    return neww

ts = np.arange(0,10,1)
somtin = np.arange(0,1,0.01)
ws = []
#ws.append(etest)
counter = 1
for s in somtin:
    ws.append(fofy(s))
    counter+=1

posps = np.arange(0,1,0.01)
plt.title("x vs y and F(y) where epislon = "+str(etest))
plt.plot(posps,ws)
plt.plot(posps,posps)
plt.show()
