import numpy as np


def pressure(data,gamma=1.4):

    datB=np.zeros(data.B.shape)
    i=0
    (a,)=data.B[:,0,0,0].shape
    (b,)=data.B[0,:,0,0].shape
    (c,)=data.B[0,0,:,0].shape
    while (i < a):
        j=0
        while (j < b):
            k=0
            while (k < c):
                if (j == (b-1)):
                    datB[i,j,k,1]=0.5*(data.B[i,j,k,1]+data.B[i,0,k,1])
                else:    
                    datB[i,j,k,1]=0.5*(data.B[i,j,k,1]+data.B[i,j+1,k,1])
                if (k == (c-1)):
                    datB[i,j,k,2]=0.5*(data.B[i,j,k,2]+data.B[i,j,0,2])
                else:
                    datB[i,j,k,2]=0.5*(data.B[i,j,k,2]+data.B[i,j,k+1,2])
                if (i == ((a-1))) :
                    if (j == (b-1)):
                        dbydy=(data.B[i,0,k,1]-data.B[i,j,k,1])/(data.y[j]-data.y[j-1])
                    else :
                        dbydy=(data.B[i,j+1,k,1]-data.B[i,j,k,1])/(data.y[j+1]-data.y[j])
                    if (k == (c-1)):
                        dbzdz=(data.B[i,j,0,2]-data.B[i,j,k,2])/(data.z[k]-data.z[k-1])
                    else :    
                        dbzdz=(data.B[i,j,k+1,2]-data.B[i,j,k,2])/(data.z[k+1]-data.z[k])
                    xiphalf=(data.x[i]-data.x[i-1])/2+data.x[i]
                    ximhalf=data.x[i]-(data.x[i]-data.x[i-1])/2
                    datB[i,j,k,0]=data.B[i-1,j,k,0]*ximhalf/xiphalf-(xiphalf**2-ximhalf**2)/2*(1/xiphalf*(dbydy/data.x[i]+dbzdz))
                else:
                    datB[i,j,k,0]=0.5*(data.B[i,j,k,0]+data.B[i+1,j,k,0])
                k += 1
            j += 1

        i += 1    
    #     ciso=np.zeros(nx)
    #    ciso2=c0**2/data.x**q
    #    pressure=data.rho[:,:,:]*ciso2
    pressure=0.5*(2*data.E[:,:,:]-datB[:,:,:,0]**2-datB[:,:,:,1]**2-datB[:,:,:,2]**2-(data.rhou[:,:,:,0]**2+data.rhou[:,:,:,1]**2+data.rhou[:,:,:,2]**2)/data.rho[:,:,:])*(gamma-1)

    return pressure
