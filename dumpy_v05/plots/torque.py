from pylab import *
import numpy as np

import pickle

class Torque:

    def __init__(self,filename='torque.txt'):
        
        #read basic file data
        (self.time,self.mi,self.mo,self.t1i,self.t1o,self.t2i,self.t2o,self.t3i,self.t3o)=\
            np.loadtxt(filename,unpack='true')


    def GetMean(self):

        return self.t1i+self.t1o#+self.t2i+self.t3i+self.t2o+self.t3o

    def rMean(self,save=0):

        torque_i=self.t1i#+self.t2i+self.t3i
        torque_o=self.t1o#+self.t2o+self.t3o
        tot_torque=torque_o+torque_i
        n=shape(tot_torque)[0]

        for i in range(n-1):
            tot_torque[i+1]=tot_torque[i]+tot_torque[i+1]
        tot_torque=tot_torque/(arange(n)+1)

        print "Final running mean average: ",tot_torque[n-1]

        if (save==1):
            file=open('alpha.txt','w')
            for i in range(shape(x)[0]):
                chaine=str(x[i])+' '+str(alpha[i])+'\n'
                file.write(chaine)
            file.close()
            

        return tot_torque

def GetMeanTorque(filename='2d_avg_snapshot.pic'):

    f=file(filename,'rb')
    x    =pickle.load(f)
    y    =pickle.load(f)
    rho  =pickle.load(f)
    f.close()

    nx=shape(x)[0] ; ny=shape(y)[0]
    rpl=3. ; mass=1. ; phipl=pi/2. ; mpl=3.e-4 ; eps=0.06 ; dz=0.6
    Roche_l = rpl * (mpl / (3. * mass) )**(1./3.)

    torque1i=np.zeros(nx) ; torque1o=np.zeros(nx)
    torque2i=np.zeros(nx) ; torque2o=np.zeros(nx)
    torque3i=np.zeros(nx) ; torque3o=np.zeros(nx)

    x_pl=rpl*cos(phipl)
    y_pl=rpl*sin(phipl)
    for j in range(ny):
        for i in range(nx):
           xloc=x[i]*cos(y[j]) ; yloc=x[i]*sin(y[j])
           ypl =x[i]*sin(y[j]-phipl)
           dist=sqrt((xloc-x_pl)**2+(yloc-y_pl)**2)
           r=x[i] ; dr=x[1]-x[0] ; dphi=y[1]-y[0]
           dtorque=rho[i,j]*ypl*rpl*r*dr*dphi*dz/(dist*dist+eps*eps)**(1.5)
           if (x[i]<rpl):
              if (dist>Roche_l):
                  torque1i[i] = torque1i[i] + dtorque
              if ((dist>0.5*Roche_l) and (dist<Roche_l)):
                  torque2i[i] = torque2i[i] + dtorque
              if (dist<=0.5*Roche_l):
                  torque3i[i] = torque3i[i] + dtorque
           else:
              if (dist>Roche_l):
                  torque1o[i] = torque1o[i] + dtorque
              if ((dist>0.5*Roche_l) and (dist<Roche_l)):
                  torque2o[i] = torque2o[i] + dtorque
              if (dist<=0.5*Roche_l):
                  torque3o[i] = torque3o[i] + dtorque

    print np.sum(torque1i+torque1o),np.sum(torque1i+torque1o+torque2i+torque2o+torque3i+torque3o)

    figure(1)
    plot(x,torque1i+torque1o)
    #plot(x,torque1i+torque1o+torque2i+torque2o+torque3i+torque3o)
    xlim((2,4))

    figure(2)
    plot(x,np.cumsum(torque1i+torque1o))
             
