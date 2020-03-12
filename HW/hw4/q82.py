import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
theta = 0.45

def setValue(ogpos,rows,cols,grid):
    ogx = ogpos[0]
    ogy = ogpos[1]
    y = 1# grid[ogx][ogy] might need to use the sample that is there? i dont htink so tho
    tosum = 0
    tosum2 = 0
    for r in rows:
        tosum+=theta*grid[r][ogy]*y
        tosum2+=theta*grid[r][ogy]*y*-1

    for c in cols:
        tosum+=theta*grid[ogx][c]*y
        tosum2+=theta*grid[r][ogy]*y*-1

    top = np.exp(tosum)
    pt2 = np.exp(tosum2)
    bot = top+pt2
    prob = top/bot
    sample = random.random()
    if(sample<prob):
        grid[ogx][ogy]=1
    else:
        grid[ogx][ogy]=-1



def updateIndexHelper(ogpos,rows,col,grid):
    if(col>0 and col < grid.shape[1]-1):# we are in the middle
        setValue(ogpos,rows,[col-1,col+1],grid)
    elif(col >0): #we are at the bottom
        setValue(ogpos,rows,[col-1],grid)
    else: # we are at the bottom
        setValue(ogpos,rows,[col+1],grid)



def updateIndex(row,col,grid):
    if(row>0 and row < grid.shape[0]-1):# we are in the middle
        updateIndexHelper((row,col),[row-1,row+1],col,grid)
    elif(row >0): #we are at the bottom
        updateIndexHelper((row,col),[row-1],col,grid)
    elif(row < grid.shape[0]-1): # we are at the bottom
        updateIndexHelper((row,col),[row+1],col,grid)
    else:
        print("wuh")

gridDim = 60
grid = np.zeros((gridDim,gridDim))
numIters = 3600001
##init grid
r = 0
c = 0
while(r<gridDim):
    c=0
    while(c<gridDim):
        if((int(random.random()+0.5))==0):
            grid[r][c]=-1
        else:
            grid[r][c]=1
        c+=1
    r+=1

counter = 0
#uniform sample the node we want to check
print(grid)
posVals = [-1,1]

fig = plt.figure(figsize=(10, 10))
SHOW_EVERY=360000
while(counter<numIters):
    row = int(random.random()*gridDim)
    col = int(random.random()*gridDim)
    updateIndex(row,col,grid)
    if((counter)%360000==0):
        plt.subplot(4, 3, int((counter)/360000)+1)
        im = plt.imshow(grid,cmap='jet')
        #colors = [ im.cmap(im.norm(value)) for value in posVals]
        #patches = [ mpatches.Patch(color=colors[i], label="{l}".format(l=posVals[i]) ) for i in range(len(posVals)) ]
        #plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
    counter+=1
plt.show()
