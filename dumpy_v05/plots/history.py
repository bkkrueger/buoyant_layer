from pylab import *
import numpy as np

class HistDisk:

    def __init__(self,filename='history.txt'):
        
        #read basic file data
        (self.time,self.dt,self.totmass,self.totbx2,self.totby2,self.totbz2)=\
            np.loadtxt(filename,unpack='true',usecols=(0,1,2,3,4,5))
        (self.maxwell,self.reynolds,self.alpha,self.divb,self.meanbx,self.meanby,self.meanbz)=\
            np.loadtxt(filename,unpack='true',usecols=(6,7,8,9,10,11,12))


class HistBox:

    def __init__(self,filename='history.txt'):
        
        #read basic file data
        (self.time,self.dt,self.totmass,self.maxwell,self.reynolds,self.alpha)=\
            np.loadtxt(filename,unpack='true',usecols=(0,1,2,3,4,5))

def alpha(data,Torb=2.*pi*3**1.5):

    plot(data.time/Torb,data.alpha/1.e-6)
    plot(data.time/Torb,data.maxwell/1.e-6)
    plot(data.time/Torb,data.reynolds/1.e-6)
    xlabel('Time (in orbits)')
    ylabel('Alpha')

    Tmin=20 ; istart=where(data.time/Torb>Tmin)[0][0]
    print 'Reynolds=',np.mean(data.reynolds[istart:]/1.e-6)
    print 'Maxwell =',np.mean(data.maxwell[istart:]/1.e-6)
    print 'alpha   =',np.mean(data.alpha[istart:]/1.e-6)
    print 'ratio   =',np.mean(data.maxwell[istart:]/data.reynolds[istart:])

def maxwell(data,Torb=2.*pi*3**1.5):

    plot(data.time/Torb,data.maxwell/1.e-6)
    xlabel('Time (in orbits)')
    ylabel('Maxwell stress')

def divb(data,Torb=2.*pi*3**1.5):

    plot(data.time/Torb,data.divb)
    xlabel('Time (in orbits)')
    ylabel('div(B)')

def mass(data,Torb=2.*pi*3**1.5,thick=1.0,norm=1.):

    plot(data.time/Torb,data.totmass/data.totmass[0]*norm,linewidth=thick)
    xlabel('Time (in orbits)')
    ylabel('M/M0')

def getMassLossRate(data,Lz=10.):
    "Calculate mass loss rate from history file"

    sizeOfData=shape(data.time)[0]
    lossRate=zeros(sizeOfData)
    for i in range(sizeOfData-1):
        lossRate[i]=data.totmass[i+1]-data.totmass[i]
        lossRate[i]=-lossRate[i]/(data.time[i+1]-data.time[i])
        if (data.time[i+1]-data.time[i]==0):
            print 'fuck'
    lossRate[sizeOfData-1]=lossRate[sizeOfData-2]

    return lossRate*Lz #Lz accounts for a bug in the history data

def getRunningMassLossRate(massLossRate,istart=1):
    "Calculate running time averaged mass loss rate from history file"

    c0=1.e-3

    sizeOfData=shape(massLossRate)[0]
    runningMassLossRate=zeros(sizeOfData-istart)
    runningMassLossRate[0]=massLossRate[istart-1]
    for i in range(1,sizeOfData-istart):
        runningMassLossRate[i]=runningMassLossRate[i-1]*i+massLossRate[i+istart-1]
        runningMassLossRate[i]=runningMassLossRate[i]/(i+1)

    finalValue=runningMassLossRate[sizeOfData-istart-1]
    print 'Running Time Average Mass Loss Rate: ',finalValue/c0

    return runningMassLossRate

def massLossRate(data,Torb=2.*pi*3**1.5,Lz=10.,thick=1.0,norm=1.):
    "Plot mass loss rate for a run"
    
    c0=1.e-3

    massLossRate=getMassLossRate(data,Lz=Lz)

    plot(data.time/Torb,massLossRate/c0*norm,linewidth=thick)
    xlabel('Time (in orbits)')
    ylabel('Mass loss rate')

def runningMassLossRate(data,Torb=2.*pi*3**1.5,Tmin=10.,Lz=10.,thick=1.0,norm=1.):
    "Plot running time average mass loss rate for a run"
    
    c0=1.e-3

    sizeOfData=shape(data.time)[0]
    istart=where(data.time/Torb>Tmin)[0][0]
    massLossRate=getMassLossRate(data,Lz=Lz)
    runningMassLossRate=getRunningMassLossRate(massLossRate,istart=istart)

    plot(data.time[istart:]/Torb,runningMassLossRate/c0*norm,linewidth=thick)
    xlabel('Time (in orbits)')
    ylabel('Mass loss rate')
    
