import matplotlib.pyplot as plt
import numpy as np

epsilon = 0.01

def calcTri(xi,xj):
    if(xi==xj):
        return 0.9
    else:
        return 0.1

def visualizeDiff():
    backgroundFile=('./iter1_0_01.png')
    foregroundFile=('./iter2_0_01.png')
    background = plt.imread(backgroundFile)
    foreground = plt.imread(foregroundFile)

    print(foreground)
    plt.figure()
    plt.imshow(background)
    plt.show()

def calcPhi(y,x,means,covs):
    pMean = means[x]
    pCovs = covs[x]
    front = 1/((2*np.pi)**(3/2) * np.linalg.det(pCovs)**(0.5))
    part1 = -1/2 * (y-pMean)
    part2 = np.linalg.inv(pCovs)
    part3 = (y-pMean)

def beliefProp(image,FGBG,means,covs):
    width = FGBG.shape[1]
    height = FGBG.shape[0]
    numIters = [1,2,3,4]
    betavals = np.linspace(0.2, 3, 15)
    l = 10
    E = 4*l*(l-1)
    Tmaxvals = np.array([15, 25])
    delta = np.zeros((betavals.size, Tmaxvals.size))
    # bhorfixed = np.random.uniform(size=(l, l-1))-0.5
    # bvertfixed = np.random.uniform(size=(l-1, l))-0.5
    # bnodesfixed = np.random.uniform(size=(l, l))-0.5
    '''
    bhorfixed = np.random.uniform(low=-1,high=1,size=(l, l-1))
    bvertfixed = np.random.uniform(low=-1,high=1,size=(l-1, l))
    bnodesfixed = np.random.uniform(low=-1,high=1,size=(l, l))
    '''
    bhorfixed = np.random.uniform(size=(l, l-1))
    bvertfixed = np.random.uniform(size=(l-1, l))
    bnodesfixed = np.random.uniform(size=(l, l))

    for iter1 in numIters:
        beta = betavals[iter1]

        # uniform random initialization of theta's in [0,beta]
        '''
        no need for these thetas anymore
        '''
        #bhor = bhorfixed*beta
        #bvert = bvertfixed*beta
        #bnodes = bnodesfixed*beta

        # The following four matrices contain the messages (in log-likelihood format)
        # such that hhorright[i,j] = (1/2)log(nu_{(i,j)->(i,j+1)}(+1)/nu_{(i,j)->(i,j+1)}(-1))
        hhorright = np.zeros((height, width-1))  # hhorright[i,j] is the message from (i,j) to (i,j+1)
        hhorleft  = np.zeros((height, width-1))  # hhorleft[i,j] is the message from (i,j+1) to (i,j)
        hvertup   = np.zeros((height-1, width))  # hvertup[i,j] is the message from (i+1,j) to (i,j)
        hvertdown = np.zeros((height-1, width))  # hvertdown[i,j] is the message from (i,j) to (i+1,j)

        # The following matrices store the updated messages,
        hhorrightN = np.zeros((height, width-1))
        hhorleftN = np.zeros((height, width-1))
        hvertupN   = np.zeros((height-1, width))
        hvertdownN = np.zeros((height-1, width))


        counter = 0
        while counter < iter1:
            # lets first take care of the corners
            #hhorright[0,0] # this is the message from the topleft corner of the image to its neighbor on the right
            hhorright[0,0] =0
            # looping over edges with 3 inputs
            for i in range(1,l-1):
                for j in range(1, l-1):
                    hhorrightnew[i, j] = bnodes[i, j] + np.sum(np.arctanh(np.array([np.tanh(bhor[i,j-1])*np.tanh(hhorright[i, j-1]), np.tanh(bvert[i-1,j])*np.tanh(hvertdown[i-1, j]), np.tanh(bvert[i, j])*np.tanh(hvertup[i, j])])))
                    hvertdownnew[i, j] = bnodes[i, j] + np.sum(np.arctanh(np.array([np.tanh(bvert[i-1,j])*np.tanh(hvertdown[i-1, j]), np.tanh(bhor[i, j-1])*np.tanh(hhorright[i, j-1]), np.tanh(bhor[i, j])*np.tanh(hhorleft[i,j])])))

            for i in range(1,l-1):
                for j in range(l-2):
                    hhorleftnew[i, j] = bnodes[i, j+1] + np.sum(np.arctanh(np.array([np.tanh(bhor[i, j+1])*np.tanh(hhorleft[i, j+1]), np.tanh(bvert[i-1, j+1])*np.tanh(hvertdown[i-1, j+1]), np.tanh(bvert[i, j+1])*np.tanh(hvertup[i, j+1])])))
                    hvertupnew[j, i]  = bnodes[j+1, i] + np.sum(np.arctanh(np.array([np.tanh(bvert[j+1, i])*np.tanh(hvertup[j+1, i]), np.tanh(bhor[j+1, i-1])*np.tanh(hhorright[j+1, i-1]), np.tanh(bhor[j+1, i])*np.tanh(hhorleft[j+1, i])])))

            # looping over edges with 2 inputs
            for i in range(1,l-1):
                hvertupnew  [-1,  i] = bnodes[-1, i] + np.sum(np.arctanh(np.array([np.tanh(bhor [ -1, i-1])*np.tanh(hhorright[ -1, i-1]), np.tanh(bhor [-1,  i])*np.tanh(hhorleft[-1,  i])])))
                hvertdownnew[ 0,  i] = bnodes[ 0, i] + np.sum(np.arctanh(np.array([np.tanh(bhor [  0, i-1])*np.tanh(hhorright[  0, i-1]), np.tanh(bhor [ 0,  i])*np.tanh(hhorleft[ 0,  i])])))
                hhorrightnew[ i,  0] = bnodes[ i, 0] + np.sum(np.arctanh(np.array([np.tanh(bvert[i-1,   0])*np.tanh(hvertdown[i-1,   0]), np.tanh(bvert[ i,  0])*np.tanh(hvertup [ i,  0])])))
                hhorleftnew [ i, -1] = bnodes[ i,-1] + np.sum(np.arctanh(np.array([np.tanh(bvert[i-1,  -1])*np.tanh(hvertdown[i-1,  -1]), np.tanh(bvert[ i, -1])*np.tanh(hvertup [ i, -1])])))

            for i in range(1,l-1):
                hvertdownnew[ i,  0] = bnodes[ i, 0] + np.sum(np.arctanh(np.array([np.tanh(bvert[i-1,   0])*np.tanh(hvertdown[i-1,   0]), np.tanh(bhor [ i,  0])*np.tanh(hhorleft [  i,  0])])))
                hvertdownnew[ i, -1] = bnodes[ i,-1] + np.sum(np.arctanh(np.array([np.tanh(bvert[i-1,  -1])*np.tanh(hvertdown[i-1,  -1]), np.tanh(bhor [ i, -1])*np.tanh(hhorright[i-1, -1])])))
                hhorrightnew[ 0,  i] = bnodes[ 0, i] + np.sum(np.arctanh(np.array([np.tanh(bhor [  0, i-1])*np.tanh(hhorright[  0, i-1]), np.tanh(bvert[ 0,  i])*np.tanh(hvertup  [  0,  i])])))
                hhorrightnew[-1,  i] = bnodes[-1, i] + np.sum(np.arctanh(np.array([np.tanh(bhor [ -1, i-1])*np.tanh(hhorright[ -1, i-1]), np.tanh(bvert[-1,  i])*np.tanh(hvertdown[ -1,  i])])))

            for i in range(l-2):
                hvertupnew [ i,  0] = bnodes[i+1,  0] + np.sum(np.arctanh(np.array([np.tanh(bvert[i+1,   0])*np.tanh(hvertup [i+1,   0]), np.tanh(bhor [i+1,   0])*np.tanh(hhorleft [i+1,   0])])))
                hvertupnew [ i, -1] = bnodes[i+1, -1] + np.sum(np.arctanh(np.array([np.tanh(bvert[i+1,  -1])*np.tanh(hvertup [i+1,  -1]), np.tanh(bhor [i+1,  -1])*np.tanh(hhorright[i+1,  -1])])))
                hhorleftnew[ 0,  i] = bnodes[ 0, i+1] + np.sum(np.arctanh(np.array([np.tanh(bhor [  0, i+1])*np.tanh(hhorleft[  0, i+1]), np.tanh(bvert[  0, i+1])*np.tanh(hvertup  [  0, i+1])])))
                hhorleftnew[-1,  i] = bnodes[-1, i+1] + np.sum(np.arctanh(np.array([np.tanh(bhor [ -1, i+1])*np.tanh(hhorleft[ -1, i+1]), np.tanh(bvert[ -1, i+1])*np.tanh(hvertdown[ -1, i+1])])))

            # Edges with single input
            hvertdownnew[  0,  0] = bnodes[ 0,  0] + np.arctanh(np.tanh(bhor[0, 0])*np.tanh(hhorleft[0, 0]))
            hvertdownnew[  0, -1] = bnodes[ 0, -1] + np.arctanh(np.tanh(bhor[-2, 0]))
            hvertupnew  [ -1,  0] = bnodes[-1,  0] + np.arctanh(np.tanh(bhor[-1, 0])*np.tanh(hhorleft  [-1,  0]))
            hvertupnew  [ -1, -1] = bnodes[-1, -1] + np.arctanh(np.tanh(bhor[-1, -1])*np.tanh(hhorright[-1, -1]))

            hhorrightnew[ 0,  0] = bnodes[ 0,  0] + np.arctanh(np.tanh(bvert[ 0,  0])*np.tanh(hvertup  [ 0,  0]))
            hhorrightnew[-1,  0] = bnodes[-1,  0] + np.arctanh(np.tanh(bvert[-1,  0])*np.tanh(hvertdown[-1, -1]))
            hhorleftnew[  0, -1] = bnodes[ 0, -1] + np.arctanh(np.tanh(bvert[ 0, -2])*np.tanh(hvertup  [ 0, -2]))
            hhorleftnew[ -1, -1] = bnodes[-1, -1] + np.arctanh(np.tanh(bvert[-1, -1])*np.tanh(hvertdown[-1, -1]))

            # convert from log-likelihoods to probabilities
            nuvertdownnew = np.exp(2*hvertdownnew)/(1+np.exp(2*hvertdownnew))
            nuvertupnew   = np.exp(2*hvertupnew)  /(1+np.exp(2*hvertupnew))
            nuhorrightnew = np.exp(2*hhorrightnew)/(1+np.exp(2*hhorrightnew))
            nuhorleftnew  = np.exp(2*hhorleftnew) /(1+np.exp(2*hhorleftnew))

            nuvertdown = np.exp(2*hvertdown)/(1+np.exp(2*hvertdown))
            nuvertup   = np.exp(2*hvertup)  /(1+np.exp(2*hvertup))
            nuhorright = np.exp(2*hhorright)/(1+np.exp(2*hhorright))
            nuhorleft  = np.exp(2*hhorleft) /(1+np.exp(2*hhorleft))

            # store the new messages as current messages
            hvertdown = hvertdownnew.copy()
            hvertup   = hvertupnew.copy()
            hhorleft  = hhorleftnew.copy()
            hhorright = hhorrightnew.copy()
            #print(bvert)

            # compute the difference delta
            if any(t+1==Tmaxvals):
                iter2 = np.where(t+1==Tmaxvals)
                delta[iter1, iter2] = (np.sum(abs(nuvertdownnew-nuvertdown)) + np.sum(abs(nuvertupnew-nuvertup)) + np.sum(abs(nuhorrightnew-nuhorright)) + np.sum(abs(nuhorleftnew-nuhorleft)))/E




def calcMeans(FGBG,knownImage):
    r=0
    c=0
    filterVals = [0,1] # first background then foreground
    means = []
    for v in filterVals:
        r=0
        c=0
        tempMean=np.zeros(3)
        counter=0
        for y in FGBG: # this grabs a row
            c=0
            for x in y:
                if(x == v):
                    tempMean+=knownImage[r][c]
                    counter+=1
                c+=1
            r+=1
        means.append(tempMean/counter)
    return means

def calcCovariances(FGBG,knownImage,means):
    r=0
    c=0
    filterVals = [0,1] # first background then foreground
    covs=[]
    for v in filterVals:
        r=0
        c=0
        tempCov=np.zeros((3,3))
        counter=0
        for y in FGBG: # this grabs a row
            c=0
            for x in y:
                if(x == v):
                    tempCov += np.dot(np.reshape(knownImage[r][c]-means[v],(3,1)),np.reshape(knownImage[r][c]-means[v],(1,3)))
                    counter+=1
                c+=1
            r+=1
        covs.append(tempCov/(counter-1))
    return covs

imageFile = ('./flower.bmp')
backgroundFile=('./background.bmp')
foregroundFile=('./foreground.bmp')
image = plt.imread(imageFile)/255 # divide by 255 to normalize our data
background = plt.imread(backgroundFile)
foreground = plt.imread(foregroundFile)

'''
each image is of dimension 160x240. Only the flower has RGB values which makes it 160x240x3
background and foreground images are initially from 0 to 255, I have flattend it below to 0 or 1
'''

realizedBG = np.zeros(background.shape) # this will flatten our image into a 0,1 matrix
realizedFG = np.zeros(foreground.shape) # this will also flatten out image int oa 0,1 matrix
FGBGmatrix = np.zeros(background.shape) # this will contain information on foreground and background
knownImage = np.zeros(image.shape)
r=0
c=0
for y in background: # this grabs a row
    c=0
    for x in y:
        if(x>0):
            realizedBG[r][c]=1
            FGBGmatrix[r][c] = 0
        else:
            FGBGmatrix[r][c]=0.5
        c+=1
    r+=1
r=0
c=0
for y in foreground: # this grabs a row
    c=0
    for x in y:
        if(x>0):
            realizedFG[r][c]=1
            FGBGmatrix[r][c]=1
        c+=1
    r+=1
r=0
c=0
for y in FGBGmatrix: # this grabs a row
    c=0
    for x in y:
        if(x!=0.5):
            knownImage[r][c]=image[r][c]
        c+=1
    r+=1

means = (calcMeans(FGBGmatrix,knownImage))
print(means)
covs = (calcCovariances(FGBGmatrix,knownImage,means))
print(covs)
print(np.dot(np.linalg.inv(covs[0]),covs[0]))
'''
plt.figure()
plt.imshow(FGBGmatrix)
plt.figure()
plt.imshow(image)
plt.figure()
plt.imshow(knownImage)
plt.show()'''

visualizeDiff()
