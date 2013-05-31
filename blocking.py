#read data from file
from numpy import *
maxBlockSizeTreshold=1e4
directoryName="/home/jonathan/projectsCP/project1/data/diberyllium/"
#directoryName="/uio/arkimedes/s79/jonatdh/"
fileName="data"
extension=".dat"
nNodes=2
deltaBlockSize = 10
outFilename = "blockingbe2energy.mat"

def importdata():
        filename=directoryName+fileName+"0"+extension
        f = open(filename)
        lines = [line.strip() for line in open(filename)]
        temp = lines
        f.close()
        for i in range(1,nNodes):
                filename=directoryName+fileName+str(i)+extension
                f = open(filename)
                lines = [line.strip() for line in open(filename)]
                temp = temp + lines
                f.close()
        return double(temp)

def computeBlock(data, blocksize):
        samples = len(data)
        blocks = int(samples / blocksize)
        averageEnergyBlock = [0] * blocks
        for i in range(0,blocks):
                energySum=0
                for j in range(i*blocksize, (i+1)*blocksize):
                        energySum = energySum+data[j]
                averageEnergyBlock[i]=energySum/blocksize
        
        #calculating the mean and variance of all blocks

        E=mean(averageEnergyBlock)
        squaredEnergies=[x*x for x in averageEnergyBlock]
        E2=mean(squaredEnergies)
        sigma=sqrt((E2-E*E)/blocks)
        
        results = [blocksize, E, sigma]
        return results

def doBlocking():
        data=importdata()
        N = len(data)
        maxBlockSize= int(N/10)
        if(maxBlockSize>maxBlockSizeTreshold):
                maxBlockSize=maxBlockSizeTreshold
        
        results=zeros((maxBlockSize/deltaBlockSize,3))
        outputFileName= directoryName+outFilename

        i=1
        while(i*deltaBlockSize<=maxBlockSize):
                blocksize=i*deltaBlockSize
                temp = computeBlock(data, blocksize)
                for j in range (0,3):
                        results[i-1,j]=temp[j]
                i=i+1
        savetxt(outputFileName, results)
doBlocking()






