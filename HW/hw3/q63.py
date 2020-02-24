import numpy as np
rho = 0.39
r = rho
jMatrix = np.array(
    [[1,-r,r,r],
    [-r,1,r,0],
    [r,r,1,r],
    [r,0,r,1]]
)
jMatrixcomp1 = np.array(
    [[1,-r,r,r,0,0,0,0],
    [-r,1,0,0,r,0,0,0],
    [r,0,1,0,0,r,r,0],
    [r,0,0,1,0,0,0,r],
    [0,r,0,0,1,0,0,0],
    [0,0,r,0,0,1,0,0],
    [0,0,r,0,0,0,1,0],
    [0,0,0,r,0,0,0,1]]
)


degrees = np.array([3,2,3,2])



def createJFCTree(depth = 3):
    '''
    initial ideas a programatic way to create a matrix for a computation tree
    '''
    startNode = 0
    while(startNode<4):
        counter = 0
        layers=[]
        while(counter < depth):
            pass

def detPosDef(matrix):
    maxS = matrix.shape[0]
    size = 1
    isPosDef = True
    while(size<maxS):
        newMatrix = np.transpose(np.transpose(matrix[0:size])[0:size])
        det = np.linalg.det(newMatrix)
        if(det<=0):
            print(det)
            isPosDef = False
        size+=1
    return isPosDef

def detPosDefv2(matrix):
    maxS = matrix.shape[0]
    size = 1
    isPosDef = True
    vals = np.diag(matrix)
    counter=0
    othersum=0
    while(counter<maxS):
        cc=0
        othersum=0
        for x in matrix[counter]:
            if(cc!=counter): #dont add diags
                othersum+= abs(x)
            cc+=1
        if(othersum>=vals[counter]):
            #print(othersum)
            isPosDef=False
        counter+=1
    return isPosDef

def recursionJs(jMatrix):
    numNodes = jMatrix.shape[0]
    initJMs = np.diag(jMatrix)
    messageMatrix = np.zeros(jMatrix.shape)
    counter=0
    #init our message matrix according to J i -> j = Jii
    for x in initJMs:
        messageMatrix[counter]=np.repeat(x,numNodes)
        counter+=1
    startNode = 0
    numIters=0
    iterMax=1001#shouldnt change as we go past 2 as we have seen everyones message
    iter = 0
    calcDiffs = np.zeros(jMatrix.shape)
    #print(jMatrix)
    while(iter<iterMax):
        #print("alihgkjsgkjsdgkj")
        copyMM = messageMatrix.copy()
        startNode = 0
        while(startNode < numNodes):
            neighbor = 0
            while(neighbor<numNodes):
                if(neighbor==startNode):
                    pass
                else:
                    sum=0
                    passNode = 0
                    while(passNode<numNodes):
                        if(passNode == startNode or passNode == neighbor):
                            pass
                        else:
                            #print('start')
                            #print(startNode)
                            #print(neighbor)
                            #print(passNode)
                            #print(jMatrix[startNode][passNode]*(1/messageMatrix[passNode][startNode])*jMatrix[passNode][startNode])
                            #print('end')
                            sum+= jMatrix[startNode][passNode]*(1/messageMatrix[passNode][startNode])*jMatrix[passNode][startNode]
                        passNode+=1
                    copyMM[startNode][neighbor] = initJMs[startNode]-sum
                    #copyMM[neighbor][startNo] = initJMs[startNode]-sum

                neighbor+=1
            startNode+=1
        messageMatrix = copyMM # has to be done to prevent out of line updates
        if(iter == iterMax-3):
            calcDiffs = messageMatrix
        if(iter == iterMax-2):
            calcDiffs = calcDiffs-messageMatrix
        iter+=1
    print(calcDiffs)
    return messageMatrix

def recursionJs2(jMatrix):
    gamma = 1 # not sure if this is okay
    numNodes = jMatrix.shape[0]
    initJMs = np.diag(jMatrix)
    messageMatrix = np.zeros(jMatrix.shape)
    counter=0
    #init our message matrix according to slide 7-32
    rows = 0
    while(rows < numNodes):
        cols = 0
        while(cols<numNodes):
            if(rows == cols):
                initm = jMatrix[rows][rows]
                #initm = 1+gamma*degrees[rows]
            else:
                initm = jMatrix[rows][rows]
            messageMatrix[rows][cols]=initm
            cols+=1
        rows+=1
    startNode = 0
    numIters=0
    iterMax=3#shouldnt change as we go past 2 as we have seen everyones message
    iter = 0
    while(iter<iterMax):
        copyMM = messageMatrix.copy()
        while(startNode < numNodes):
            neighbor = 0
            while(neighbor<numNodes):
                if(neighbor==startNode):
                    pass
                else:
                    sum=0
                    passNode = 0
                    while(passNode<numNodes):
                        if(passNode == startNode or passNode == neighbor):
                            pass
                        else:
                            numer = gamma**2
                            denom = copyMM[passNode][startNode]
                            sum+= numer/denom
                        passNode+=1
                    copyMM[startNode][neighbor] = 1+gamma*degrees[startNode]-sum
                neighbor+=1
            startNode+=1
        messageMatrix = copyMM
        iter+=1
    return messageMatrix

def getMarginals(jMatrix, messageMatrix):
    numNodes = jMatrix.shape[0]
    marginals=[]
    startNode = 0
    numIters=0
    while(startNode < numNodes):
        neighbor = 0
        sum=0
        while(neighbor<numNodes):
            if(neighbor==startNode):
                pass
            else:
                passNode = 0
                sum += jMatrix[startNode][neighbor]*(1/messageMatrix[neighbor][startNode])*jMatrix[neighbor][startNode]
            neighbor+=1
        marginals.append(jMatrix[startNode][startNode]-sum)
        startNode+=1
    return marginals

print(detPosDef(jMatrix))
print(detPosDefv2(jMatrix))
print(detPosDef(jMatrixcomp1))
print(np.diag(np.linalg.inv(jMatrix)))
MM=(recursionJs(jMatrix))
print(MM) # our main objective
#MM=(recursionJs2(jMatrix))
#print(MM)
print(getMarginals(jMatrix,MM))
