import numpy as np
from pylab import *
from matplotlib.pyplot import subplots_adjust

from scipy.special import lpn
from math import sqrt

from dumpy_v05.data import rd_dumses

def PlotHist(filename='history.txt',save=0,thick=1.0,n=10,beta=100,Tend=2.5):
    pi=2.*np.arcsin(1.0) ; Omega=1.e-3 ; Torb=2.*pi/Omega
    #(time,stress,B0,rms_br,rms_bphi)=np.loadtxt(filename,usecols=(0,5,8,10,12),unpack='true')
    (time,stress,B0,rms_br,rms_bphi)=np.loadtxt(filename,usecols=(0,3,3,3,3),unpack='true')
    time=time*Torb ; stress=stress*1.e-6

    P0=1.e-6
    alpha=stress/P0
    rms_b=[sqrt(item) for item in rms_br+rms_bphi]

    time_ind=where(time/Torb<Tend)
    size_time_ind=time_ind[0].shape[0]
    
    s=GrowthRate(n,Omega,beta=beta)
    Maxwell=5.e-13*exp(2*s*time[time_ind])/1.e-6
    figure(0)
    xlabel('Time (in orbits)') ; ylabel('Maxwell stress/pressure')
    semilogy(time/Torb,stress/P0,'-',linewidth=thick)
    semilogy(time[time_ind]/Torb,Maxwell,linewidth=max(1.0,thick-1.0))
    ylim(1.e-8,1.e1)

    figure(1)
    Maxwell=5.e-7*exp(2*s*time)
    GrowthRate_num=np.zeros(size_time_ind)
    GrowthRate_th =np.zeros(size_time_ind)
    for ipos in range(size_time_ind):
        dMax   =alpha[ipos+1]-alpha[ipos]
        dMax_th=Maxwell[ipos+1]-Maxwell[ipos]
        dt  =(time[ipos+1]-time[ipos])
        GrowthRate_num[ipos]=0.5*dMax   /(alpha[ipos]+alpha[ipos+1])/0.5/dt
        GrowthRate_th [ipos]=0.5*dMax_th/(Maxwell[ipos]+Maxwell[ipos+1])/0.5/dt
    plot(time[time_ind]/Torb,    GrowthRate_num/Omega)
    plot(time[time_ind]/Torb,    GrowthRate_th /Omega)
    plot(time[time_ind]/Torb,0.8*GrowthRate_th /Omega)
    plot(time[time_ind]/Torb,0.9*GrowthRate_th /Omega)
    plot(time[time_ind]/Torb,1.1*GrowthRate_th /Omega)
    ylim(0.,1.)

    fit_end_ind   =where(alpha>1.e-4)
    Omega=1.e-3
    print 'sigma_th =',s/Omega
    print 'sigma_num=',mean(GrowthRate_num[0:fit_end_ind[0][0]])/Omega,'+/-',std(GrowthRate_num[0:fit_end_ind[0][0]])/Omega
    linear_end_ind=where(GrowthRate_num<0.9*GrowthRate_th)
    print '    dB/B0=',rms_b[linear_end_ind[0][0]]/B0[0],alpha[linear_end_ind[0][0]],time[linear_end_ind[0][0]]/Torb

    if save==1:
        savefig('hist.ps')

def PlotMaxwell(filename='history.txt',save=0,thick=1.0,n=10,beta=100,Tend=2.5):
    pi=2.*np.arcsin(1.0) ; Omega=1.e-3 ; Torb=2.*pi/Omega
    (time,stress,B0,rms_br,rms_bphi)=np.loadtxt(filename,usecols=(0,5,8,10,11),unpack='true')
    xlabel('Time (in orbits)') ; ylabel('Maxwell stress/pressure')
    semilogy(time/Torb,stress/1.e-6,'-',linewidth=thick)

    P0=1.e-6
    alpha=stress/P0
    rms_b=[sqrt(item) for item in rms_br+rms_bphi]

    time_ind=where(time/Torb<Tend)
    size_time_ind=time_ind[0].shape[0]
    
    s=GrowthRate(n,Omega,beta=beta)
    Maxwell=2.e-12*exp(2*s*time[time_ind])/1.e-6
    #semilogy(time[time_ind]/Torb,Maxwell,linewidth=max(1.0,thick-1.0))
    semilogy(time/Torb,rms_bphi,linewidth=max(1.0,thick-1.0))
    #ylim(1.e-8,1.e1)

    if save==1:
        savefig('hist.ps')

def GrowthRate(n,Omega,beta=100):
    d0=1. ; ciso=Omega ; H=ciso/Omega
    sigma=[0,1.1584483822,2.07955511702,2.98292656105,3.87983788353,4.7735461198,\
           5.66539225379,6.55604318513,7.44587302068,8.33510992911,9.22390203811,\
           10.112350311,11.0005262181,11.8884819764,12.7762567094,13.6638804027,\
           14.5513765898,15.4387639714,16.3260575284,17.2132694602,18.1004098526,\
           18.9874871359,19.8745084414,20.7614798382,21.6484064802,22.5352927846,\
           23.4221425833,24.308959224,25.1957456517,26.0825044223,26.9692377404]
    B0=sqrt(2.*d0*ciso**2/beta)
    #kn=1./H*sqrt(n+n**2)
    if n==68:
        kn=60.6561711644
    else:
        kn=sigma[n]
    va=B0/sqrt(d0)
    xi=va/H/Omega
    return Omega*sqrt((-(1.+2.*xi**2*kn**2*H**2)+sqrt(1.+16.*xi**2*kn**2*H**2))/2.)

def GetWeight(idump=1,nvar=8,nmax=5,type='vx'):
    data=rd_dumses.DumsesData(idump,nvar)
    vx=data.get_1d_y(type)
    nz=data.y.size
    Pn=np.zeros((nz,nmax+1))
    In=np.zeros(nmax+1)
    ciso=1.e-3 ; Omega0=1.e-3 ; H=ciso/Omega0 ; Lz=10.*H
    dz=Lz/nz

    if type in ('vx','vz'):
        for n in range(1,nmax+1):
            for k in range(nz):
                Pn[k,n]=lpn(n,np.tanh(data.y[k]))[0][n]
                In[n]=In[n]+vx[k]*Pn[k,n]*(1.-np.tanh(data.y[k])**2)*dz
            In[n]=(2*n+1)/2.*In[n]
    else:
        for n in range(1,nmax+1):
            for k in range(nz):
                Pn[k,n]=lpn(n,np.tanh(data.y[k]))[1][n]
                In[n]=In[n]+vx[k]*Pn[k,n]*(1.-np.tanh(data.y[k])**2)*dz
            In[n]=(2*n+1)/2./n/(n+1)*In[n]

    return In

def PlotWeight(In,save=0):

    xlabel('Mode number') ; ylabel('Mode amplitude')
    nmax=In.size-1
    maxval=max(log10(abs(In)))
    plot(range(nmax+1),log10(abs(In)),'r-')
    plot(range(nmax+1),log10(abs(In)),'ro')
    axis([0,nmax,maxval-2.5,maxval+0.5])
    if save==1:
        savefig('mode_weight.eps')


def ModesConstruct(idump=0,nrange=(0,10),nmode=0,nvar=8,type='vx',save=0):
    
    if nmode==0:
        if size(nrange)==1:
            nmodes=nrange
        else:
            nmodes=nrange[0]+np.array(range(nrange[1]-nrange[0]+1))
    if nmode==1:
        nmodes=np.array(nrange)

    #Get Data
    data=rd_dumses.DumsesData(idump,nvar)
    vx=data.get_1d_y(type)
    rho=data.get_1d_y('rho')
    if (type=='bx'):
        by=data.get_1d_y('by')
    else:
        by=1
    nz=data.y.size

    #Get Normal modes amplitudes
    if size(nmodes)==1:
        In=GetWeight(idump=idump,nmax=nmodes,type=type)
    else:
        In=GetWeight(idump=idump,nmax=max(nmodes),type=type)
    Mode=np.zeros(nz)

    #Reconstruct modes
    if size(nmodes)==1:
        mode_nb=nmodes
        if type in ('vx','vz'):
            for k in range(nz):
                z=np.tanh(data.y[k])
                Mode[k]=In[mode_nb]*lpn(mode_nb,z)[0][mode_nb]
        else:
            for k in range(nz):
                z=np.tanh(data.y[k])
                Mode[k]=In[mode_nb]*lpn(mode_nb,z)[1][mode_nb]*(1-z**2)
    else:        
        if type in ('vx','vz'):
            for n in range(size(nmodes)):
                mode_nb=nmodes[n]
                for k in range(nz):
                    z=np.tanh(data.y[k])
                    Mode[k]=Mode[k]+In[mode_nb]*lpn(mode_nb,z)[0][mode_nb]
        else:
            for n in range(size(nmodes)):
                mode_nb=nmodes[n]
                for k in range(nz):
                    z=np.tanh(data.y[k])
                    Mode[k]=Mode[k]+In[mode_nb]*lpn(mode_nb,z)[1][mode_nb]*(1-z**2)
    #Conv2beta
    #vx=conv2beta(vx,exp(-data.y**2/2.))
    #Mode=conv2beta(Mode,exp(-data.y**2/2.))

    #Plot results
    axes=Axes(figure(1),[0.25,0.25,0.5,0.5])
    xlabel('z/H') ; ylabel('Bx/Bz')
    plot(data.y,vx/by,linewidth=2.5,linestyle='-',color='b')
    plot(data.y,Mode/by,linewidth=2,linestyle='--',color='r')
    xlim((min(data.y),max(data.y)))
    subplots_adjust(left=0.15)
    #show()
    #ylim((0,30))
    if save==1:
        savefig('mode_'+type+'.eps')

def conv2beta(B,rho,ciso=1.e-3):

    P=rho*ciso**2
    Pmag=B**2/2
    return P/Pmag
