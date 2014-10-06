from pylab import *
import numpy as np
import os as os
import pickle

from dumpy_v05.data import rd_dumses
from dumpy_v05.plots import snapshots
from dumpy_v05.plots import streamlines

def xprofile(data,nvar=8,type='rho',thick=1.0,save=0,scale=1):

    #Get radial cut
    #xcut=data.get_slice_x(type)
    xcut=data.get_yz_mean()

    #Plot data
    rcParams['image.origin']='lower'
    xlabel('Radius')
    ylabel(type)
    plot(data.x,xcut)

    if save==1:
	savefig('rprofile.eps')

    return

def alpha_profile(data,save=0):

    #Get alpha radial profile
    x,alpha=GetAlpha(data)

    plot(x,alpha)
    xlabel('Radius')
    ylabel('alpha')

    if save==1:
	savefig('alpha.eps')

    return

def GetAlpha(data,q=0,c0=0.057735027):
    """ Compute alpha profile from Dumses output"""
    
    #Get dimensions
    nx,ny,nz,dim=shape(data.B)
    
    ciso=np.zeros(nx)
    ciso2=c0**2/data.x**q

    alpha_max=-yz_mean(data.B[:,:,:,0]*data.B[:,:,:,1])/yz_mean(data.rho[:,:,:])/ciso2
    vymean=yz_mean(data.rhou[:,:,:,1])/yz_mean(data.rho)
    dvy=data.rhou[:,:,:,1]/data.rho
    for i in range(nx):
        dvy[i,:,:]=dvy[i,:,:]-vymean[i]
    alpha_rey=yz_mean(data.rhou[:,:,:,0]*dvy)/yz_mean(data.rho[:,:,:])/ciso2

    alpha=alpha_max+alpha_rey

    return data.x,alpha_rey

def MeanAlpha(istart,istop,save=0):
    "Compute alpha radial profile averaged between output istart and istop"

    for i in range(istart,istop+1):
        print i
        data=rd_dumses.DumsesData(i)
        x,rd_alpha=GetAlpha(data)
        if (i==istart):
            alpha=rd_alpha
        alpha=(rd_alpha+alpha*(i-istart))/(i-istart+1)

    plot(x,alpha)
    xlabel('Radius')
    ylabel('alpha')

    if save==1:
	savefig('mean_alpha.eps')
    if save==2:
        file=open('alpha.txt','w')
        for i in range(shape(x)[0]):
            chaine=str(x[i])+' '+str(alpha[i])+'\n'
            file.write(chaine)
        file.close()

    return

def GetMeanField1d(filename="2d_avg_snapshot.pic",p=0):

    f=file(filename,'rb')
    x      =pickle.load(f)
    y      =pickle.load(f)
    field  =pickle.load(f)
    f.close()

    field1D=np.mean(field,1)

    plot(x,field1D*x**p)
    xlabel('Radius')
    ylabel('Surface Density')

def MeanDensity(istart,istop,save=1,p=0.):
    "Compute mean gas surface density averaged between istart and istop"

    for i in range(istart,istop+1):
        print i
        data=rd_dumses.DumsesData(i)
        x=data.x ; rd_rho=yz_mean(data.rho)
        if (i==istart):
            rho=rd_rho
        rho=(rd_rho+rho*(i-istart))/(i-istart+1)

    plot(x,rho*x**p)
    xlabel('Radius')
    ylabel('Surface Density')

#    x100=1.+7.*arange(100)/100.
#    plot(x100,1./x100**0.5)
#    ylim(0.,1.5)

    if save==1:
	savefig('mean_density.eps')
    if save==2:
        file=open('surf_dens.txt','w')
        for i in range(shape(x)[0]):
            chaine=str(x[i])+' '+str(rho[i])+'\n'
            file.write(chaine)
        file.close()

    return

def MeanDensity2D(istart,istop,save=1):
    "Compute mean gas surface density averaged between istart and istop"

    for i in range(istart,istop+1):
        print i
        data=rd_dumses.DumsesData(i)
        x=data.x ; y=data.y ; rd_rho=z_mean(data.rho)
        if (i==istart):
            rho=rd_rho
        rho=(rd_rho+rho*(i-istart))/(i-istart+1)

    xmin=1. ; xmax=8.
    ymin=0. ; ymax=pi

    rcParams['image.origin']='lower'
    hot()
    xlabel('Phi') ; ylabel('R')
    imshow(transpose(rho),aspect='auto',extent=(xmin,xmax,ymin,ymax))
    colorbar()

    if save==1:
	savefig('mean_density.eps')
    if save==2:
	savefig('mean_density.eps')
        f=file('2d_avg_snapshot.pic','wb')
        pickle.dump(x,f)
        pickle.dump(y,f)
        pickle.dump(rho,f)
        f.close()

    return

def MeanBy2D(istart,istop,save=1):
    "Compute 2D mean By averaged between istart and istop"

    xmin=1. ; xmax=8.
    ymin=0. ; ymax=pi

    rcParams['image.origin']='lower'
    hot()
    xlabel('Phi') ; ylabel('R')

    for i in range(istart,istop+1):
        print i
        data=rd_dumses.DumsesData(i)
        x=data.x ; y=data.y
#        rd_rho=z_mean(data.B[:,:,:,1])
#        rd_rho=z_mean(data.B[:,:,:,1]*data.B[:,:,:,1])**0.5
#        bymean=z_mean(data.B[:,:,:,1])
#        bymean_3d=data.B[:,:,:,1]
#        for k in range(shape(data.B)[2]):
#            bymean_3d[:,:,k]=bymean
#        rd_rho=(z_mean((data.B[:,:,:,1]-bymean_3d)*(data.B[:,:,:,1]-bymean_3d)))**0.5
        rd_rho=(z_mean(data.B[:,:,:,1]*data.B[:,:,:,1])-(z_mean(data.B[:,:,:,1]))**2)**0.5
        if (i==istart):
            rho=rd_rho
        rho=(rd_rho+rho*(i-istart))/(i-istart+1)
        imshow(transpose(rho),aspect='auto',extent=(xmin,xmax,ymin,ymax))
        draw()

    imshow(transpose(rho),aspect='auto',extent=(xmin,xmax,ymin,ymax))
    colorbar()

    if save==1:
	savefig('mean_by.eps')
    if save==2:
	savefig('mean_by2.eps')
        f=file('2d_avg_snapshot_by2.pic','wb')
        pickle.dump(x,f)
        pickle.dump(y,f)
        pickle.dump(rho,f)
        f.close()

    return

def GetMeanField2d(filename="2d_avg_snapshot.pic",zoom=0):

    f=file(filename,'rb')
    x      =pickle.load(f)
    y      =pickle.load(f)
    field  =pickle.load(f)
    f.close()

    xmin=min(x) ; xmax=max(x)
    ymin=min(y) ; ymax=max(y) 
    rcParams['image.origin']='lower'
    hot()
    ylabel('Phi') ; xlabel('R')
    if zoom==0:
        imshow(transpose(field),aspect='auto',extent=(xmin,xmax,ymin,ymax))
    if zoom==1:
        xmin=x[60] ; xmax=x[130]
        ymin=y[210] ; ymax=y[270] 
        imshow(transpose(field[60:130,210:270]),aspect='auto',extent=(xmin,xmax,ymin,ymax))
    colorbar()


def GetMeanDensity2d(filename="2d_avg_snapshot.txt",idump=20,p=0.):

    f=file('2d_avg_snapshot.pic','rb')
    x    =pickle.load(f)
    y    =pickle.load(f)
    rho  =pickle.load(f)
    f.close()

    figure(1)
    xmin=x[60] ; xmax=x[130]
    ymin=y[210] ; ymax=y[270] 
    rcParams['image.origin']='lower'
    hot()
    ylabel('Phi') ; xlabel('R')
    imshow(transpose(rho[60:130,210:270]),aspect='auto',extent=(xmin,xmax,ymin,ymax))
    colorbar()

    os.chdir('../../hydro_320_480_2_pi_q0_p1.5/nu0_newsurf/')
    data=rd_dumses.DumsesData(idump)

    figure(2)
    ylabel('Phi') ; xlabel('R')
    imshow(transpose(data.rho[60:130,210:270,0]),aspect='auto',extent=(xmin,xmax,ymin,ymax))
    colorbar()

    figure(3)
    ylabel('Phi') ; xlabel('R')
    imshow(transpose(rho-data.rho[:,:,0]),aspect='auto',extent=(xmin,xmax,ymin,ymax))
    colorbar()
    
    os.chdir('../../beta400_320_480_40_pi_q0_p1.5/wplanet/')
    #os.chdir('../../beta400_320_480_40_pi_q1/wplanet/')

    figure(4)
    plot(data.x,np.mean(data.rho[:,:,0],1)*data.x**p,color='b')
    plot(data.x,np.mean(rho,1)*data.x**p,color='r')

    figure(5) ; p=0.5
    os.chdir('../../beta400_320_480_40_pi_q1/wplanet/')
    f=file('2d_avg_snapshot.pic','rb')
    x    =pickle.load(f)
    y    =pickle.load(f)
    rho  =pickle.load(f)
    f.close()
    os.chdir('../../hydro_320_480_4_wpl_pi_q1/nu0_newsurf/')
    data=rd_dumses.DumsesData(50)
    plot(data.x,np.mean(data.rho[:,:,0],1)*data.x**p,color='b')
    plot(data.x,np.mean(rho,1)*data.x**p,color='r')

    figure(6)
    ylabel('Phi') ; xlabel('R')
    imshow(transpose(rho[60:130,210:270]),aspect='auto',extent=(xmin,xmax,ymin,ymax))
    colorbar()

    figure(7)
    ylabel('Phi') ; xlabel('R')
    imshow(transpose(data.rho[60:130,210:270,0]),aspect='auto',extent=(xmin,xmax,ymin,ymax))
    colorbar()

    os.chdir('../../beta400_320_480_40_pi_q0_p1.5/wplanet/')

    return

def MeanBy(istart,istop,save=1,color='r'):
    "Compute mean By averaged between istart and istop"

    for i in range(istart,istop+1):
        print i
        data=rd_dumses.DumsesData(i)
        x=data.x ; rd_by=yz_mean(data.B[:,:,:,1])
        if (i==istart):
            by=rd_by
        by=(rd_by+by*(i-istart))/(i-istart+1)

    plot(x,by,color)
    xlabel('Radius')
    ylabel('By')

    if save==1:
	savefig('mean_density.eps')

    return

def ubAngle(istart,istop,Lz=10,binWidth=0.5):
    "Compute angle between mean velocity and magnetic field"

    Nbin=int(Lz/binWidth)
    meanAbsTheta=zeros(Nbin)
    zBin=binWidth/2.+arange(-Nbin/2,Nbin/2)*binWidth ; left=arange(-Nbin/2,Nbin/2)*binWidth

    for ifile in range(istart,istop+1):
        print "Reading file number ",ifile," ..."
        data=rd_dumses.DumsesData(ifile)
        print "  File number ",ifile," read."

        bx=xy_mean(data.B[:,:,:,0])
        by=xy_mean(data.B[:,:,:,1])
        vx=xy_mean(data.rhou[:,:,:,0]/data.rho[:,:,:])
        vy=xy_mean(data.rhou[:,:,:,1]/data.rho[:,:,:])
        
        v=sqrt(vx*vx+vy*vy)
        b=sqrt(bx*bx+by*by)

        cosTheta=(bx*vx+by*vy)/v/b
        sinTheta=(vx*by-vy*bx)#/v/b

        #figure(1)
        #plot(data.z,bx)
        #plot(data.z,by)
        #figure(2)
        #plot(data.z,vx)
        #plot(data.z,vy)
        #plot(data.z,data.B[0,0,:,0])
        #stop

        theta=np.arccos(cosTheta)

        nz=shape(data.z)[0]
        for k in range(nz):
            if sinTheta[k]<0:
                theta[k]=-theta[k]

        for ibin in range(Nbin):
            kmin=where(data.z>(-Nbin/2.+ibin)*binWidth)[0][0]
            if ibin<Nbin-1:
                kmax=where(data.z>(-Nbin/2.+ibin+1)*binWidth)[0][0]
            else:
                kmax=nz
            meanAbsTheta[ibin]=meanAbsTheta[ibin]*(ifile-istart)+abs(np.mean(theta[kmin:kmax]))
            meanAbsTheta[ibin]=meanAbsTheta[ibin]/(ifile-istart+1)

        if ifile==istart:
            meanTheta=theta

        meanTheta=(theta + meanTheta*(ifile-istart))/(ifile-istart+1)

    figure(1)
    plot(data.z,meanTheta)
    figure(2)
    #plot(zBin,meanAbsTheta,'ro')
    bar(left,meanAbsTheta,width=binWidth,color='b')
    ylim(0,2.) ; xlim(-5,5)
    xlabel('Z/H') ; ylabel('Angle')
    #plot(data.z,cosTheta)
    #figure(4)
    #plot(data.z,v*b)
    #plot(data.z,sinTheta)

def byVertStruct(istart,istop,save=0,readPickle=False,color='r'):
    "Compute horizontally and time averaged mean By and fluctuations"

    matplotlib.rc('text', usetex=True)
    matplotlib.rcParams['font.size']=16
    rho0=1. ; c0sq=1.e-6 ; Pth=rho0*c0sq

    for i in range(istart,istop+1):
        print i
        if (readPickle==True):
            f2=file('betaProf'+str(i)+'.pickle','rb')
            z=pickle.load(f2) ; nz=size(z)
            mBy2=pickle.load(f2)
            dBy2=pickle.load(f2)
            mBeta=pickle.load(f2)
            f2.close()
            meanBy2=mBy2 ; deltaBy2=dBy2 ; meanBeta=mBeta 
        else:
            data=rd_dumses.DumsesData(i)
            mBy  =xy_mean(data.B[:,:,:,1]) ; mBy2=mBy**2
            dBy2 =xy_mean(data.B[:,:,:,1]**2)
            mBeta=2.*xy_mean(data.rho)*c0sq/xy_mean(data.B[:,:,:,0]**2+data.B[:,:,:,1]**2+data.B[:,:,:,2]**2)
            if (i==istart):
                z=data.z ; nz=size(z)
                meanBy2 =mBy**2
                deltaBy2=dBy2-mBy**2
                meanBeta=mBeta

        zmid=np.int(shape(meanBeta)[0]/2.)

        #dump beta profile in betaProf.pickle file
        if (readPickle==False):
            f2=file('betaProf'+str(i)+'.pickle','wb')
            pickle.dump(data.z,f2)
            pickle.dump(mBy**2,f2)
            pickle.dump((dBy2-mBy**2),f2)
            pickle.dump(mBeta,f2)
            f2.close()

        #Compute time average
        if (readPickle==False):
            meanBy2 =( mBy2      + meanBy2*(i-istart))/(i-istart+1)
            deltaBy2=((dBy2-mBy2)+deltaBy2*(i-istart))/(i-istart+1)
            meanBeta=(mBeta      +meanBeta*(i-istart))/(i-istart+1)
        else:
            meanBy2 =( mBy2      + meanBy2*(i-istart))/(i-istart+1)
            deltaBy2=(dBy2       +deltaBy2*(i-istart))/(i-istart+1)
            meanBeta=(mBeta      +meanBeta*(i-istart))/(i-istart+1)

    figure(1)
#    semilogy(z, meanBy2/Pth,'r',linewidth=2.)
#    semilogy(z,deltaBy2/Pth,'b',linewidth=2.)
    plot(z, meanBy2/Pth,'r',linewidth=2.)
    plot(z,deltaBy2/Pth,'b',linewidth=2.)
    xlabel('Z/H')
    ylabel('$<$By$>$$^2$,dBy$^2$')

    figure(2)
    plot(z, sqrt(meanBy2/deltaBy2),linewidth=2.,color=color)
    xlabel('Z/H')
    ylabel('$<$By$>$/dBy')

    figure(3)
    semilogy(z,meanBeta,linewidth=2.,color=color)
    xlabel('Z/H')
    ylabel(r"$\beta$")
    zbot,ztop=getBetaOne(meanBeta,z,nz/2)
    print zbot,ztop
    plotArrow(zbot,0.12,0.2,color=color)
    plotArrow(ztop,0.12,0.2,color=color)

    if save==1:
	savefig('byStruct.eps')

    return

def fieldlines(istart,istop,readPickle=False):
    
    matplotlib.rc('text', usetex=True)
    matplotlib.rcParams['font.size']=16

    for i in range(istart,istop+1):
        print i
        if readPickle:
            f=open('meany'+str(i)+'.pickle','rb')
            x=pickle.load(f)
            z=pickle.load(f)
            rho=pickle.load(f)
            bx=pickle.load(f)
            by=pickle.load(f)
            bz=pickle.load(f)
            vx=pickle.load(f)
            vz=pickle.load(f)
            f.close()
        else:
            data=rd_dumses.DumsesData(i)
            x=data.x ; z=data.z
            rho=np.mean(data.rho,1)
            bx=np.mean(data.B[:,:,:,0],1)
            by=np.mean(data.B[:,:,:,1],1)
            bz=np.mean(data.B[:,:,:,2],1)
            vx=np.mean(data.rhou[:,:,:,0]/data.rho,1)
            vz=np.mean(data.rhou[:,:,:,2]/data.rho,1)
        if (i==istart):
            mbx=bx ; mby=by ; mbz=bz
            mvx=vx ; mvz=vz
        else:
            mbx=mbx+bx ; mby=mby+by ; mbz=mbz+bz
            mvx=mvx+vx ; mvz=mvz+vz

        #dump stress profile in stressProf.pickle file
        if not(readPickle):
            f=file('meany'+str(i)+'.pickle','wb')
            pickle.dump(x,f)
            pickle.dump(z,f)
            pickle.dump(rho,f)
            pickle.dump(bx,f)
            pickle.dump(by,f)
            pickle.dump(bz,f)
            pickle.dump(vx,f)
            pickle.dump(vz,f)
            f.close()

    mbx=mbx/(istop-istart+1)
    mby=mby/(istop-istart+1)
    mbz=mbz/(istop-istart+1)
    mvx=mvx/(istop-istart+1)
    mvz=mvz/(istop-istart+1)

    nx,nz=shape(mbx)
    dx=x[1]-x[0] ; dz=z[1]-z[0]
    Ay=np.zeros((nx,nz))
    for i in range(1,nx):
        Ay[i,0]=Ay[i-1,0]+mbz[i,0]*dx
        for k in range(1,nz):
            Ay[i,k]=Ay[i,k-1]-mbx[i,k]*dz
        
    xmin=min(x) ; xmax=max(x)
    zmin=min(z) ; zmax=max(z)

    figure(1,figsize=(4,4*nz/nx))
    contour(x,z,transpose(Ay),10,colors='white')
    contour(x,z,transpose(-Ay),10,colors='white')
    contourf(x,z,transpose(mby),256,figsize=(xmin,xmax,zmin,zmax))
    xlabel('x/H') ; ylabel('z/H')
    xlim(xmin,xmax) ; ylim(zmin,zmax)
    gcf().subplots_adjust(left=0.15)
    gcf().subplots_adjust(bottom=0.05)
    gcf().subplots_adjust(top=0.95)

    figure(2,figsize=(4,4*nz/nx))
    streams=streamlines.Streamlines(x,z,transpose(mvx),transpose(mvz),res=1, spacing=10)
    streams.plot()
    contourf(x,z,transpose(rho),256,figsize=(xmin,xmax,zmin,zmax))
    xlabel('x/H') ; ylabel('z/H')
    xlim(xmin,xmax) ; ylim(zmin,zmax)
    gcf().subplots_adjust(left=0.15)
    gcf().subplots_adjust(bottom=0.05)
    gcf().subplots_adjust(top=0.95)

def computeBetaOne(istart,istop):

    rho0=1. ; c0sq=1.e-6 ; Pth=rho0*c0sq

    filename='betaOne.pickle'
    if os.path.isfile(filename):
        f=file(filename,'ab')
    else:
        f=file(filename,'wb')

    for i in range(istart,istop+1):
        print i
        data=rd_dumses.DumsesData(i)
        mBeta=2.*xy_mean(data.rho)*c0sq/xy_mean(data.B[:,:,:,0]**2+data.B[:,:,:,1]**2+data.B[:,:,:,2]**2)
        z=data.z
        zmid=np.int(shape(mBeta)[0]/2.)

        #Compute location where beta=1
        zbottom,ztop=getBetaOne(mBeta,z,zmid)
        pickle.dump((data.time,zbottom,ztop),f)

    f.close()

def getBetaOne(array,z,zstart):
    "Find location in disk where beta=1"

    zTop=zstart ; zBottom=zstart ; zMax=zstart*2
    while ((array[zTop]>1.) and (zTop<zMax-1)):
        zTop=zTop+1
    while ((array[zBottom]>1.) and (zBottom>1.)):
        zBottom=zBottom-1

    return z[zBottom],z[zTop]

def readBetaOne(filename):
    "Read pickle file containing beta=1 locations"

    time=[] ; zBottom=[] ; zTop=[]

    f=open(filename,'r')
    eOf=False
    while not(eOf):
        try:
            toto=pickle.load(f)
        except:
            print 'End of File reached.'
            f.close()
            eOf=True
        else:
            time.append(toto[0])
            zBottom.append(toto[1])
            zTop.append(toto[2])

    return time,zBottom,zTop

def get_transport(istart,istop,readPickle=False,zcut=2.5,verbose=False):
    "Compute vertical profiles BxBy and BxBz"

    matplotlib.rc('text', usetex=True)
    matplotlib.rcParams['font.size']=16

    for i in range(istart,istop+1):
        if (verbose):
            print i
        if readPickle:
            f=open('stressProf'+str(i)+'.pickle','rb')
            z=pickle.load(f) ; nz=size(z)
            alphax=pickle.load(f)
            alphaz=pickle.load(f)
            windx=pickle.load(f)
            windz=pickle.load(f)
            fwind=pickle.load(f)
            fturb=pickle.load(f)
            reyx=pickle.load(f)
            reyz=pickle.load(f)
            f.close()
        else:
            data=rd_dumses.DumsesData(i)
            bx=xy_mean(data.B[:,:,:,0])
            by=xy_mean(data.B[:,:,:,1])
            bz=xy_mean(data.B[:,:,:,2])
            windx=bx*by
            windz=bz*by
            reyx=xy_mean(data.rhou[:,:,:,0]*data.rhou[:,:,:,1]/data.rho[:,:,:])
            reyz=xy_mean(data.rhou[:,:,:,0]*data.rhou[:,:,:,2]/data.rho[:,:,:])
            alphax=xy_mean((data.B[:,:,:,0]-bx)*(data.B[:,:,:,1]-by))
            alphaz=xy_mean((data.B[:,:,:,2]-bz)*(data.B[:,:,:,1]-by))
            fwind=sqrt(windx*windx+windz*windz)
            fturb=sqrt(alphax*alphax+alphaz*alphaz)
            z=data.z ; nz=size(z)
        if (i==istart):
            mwindx=windx
            mwindz=windz
            mreyx=reyx
            mreyz=reyz
            malphax=alphax
            malphaz=alphaz
            mfwind=fwind
            mfturb=fturb
        else:
            mwindx=mwindx+windx
            mwindz=mwindz+windz
            mreyx=mreyx+reyx
            mreyz=mreyz+reyz
            malphax=malphax+alphax
            malphaz=malphaz+alphaz
            mfwind=mfwind+fwind
            mfturb=mfturb+fturb

        #dump stress profile in stressProf.pickle file
        if not(readPickle):
            f=file('stressProf'+str(i)+'.pickle','wb')
            pickle.dump(data.z,f)
            pickle.dump(alphax,f)
            pickle.dump(alphaz,f)
            pickle.dump(windx,f)
            pickle.dump(windz,f)
            pickle.dump(fwind,f)
            pickle.dump(fturb,f)
            pickle.dump(reyx,f)
            pickle.dump(reyz,f)
            f.close()

    mwindx=mwindx/(istop-istart+1)
    mwindz=mwindz/(istop-istart+1)
    malphax=malphax/(istop-istart+1)
    malphaz=malphaz/(istop-istart+1)
    mfwind=mfwind/(istop-istart+1)
    mfturb=mfturb/(istop-istart+1)
    mreyx=mreyx/(istop-istart+1)
    mreyz=mreyz/(istop-istart+1)

    P0=1.e-6

    figure(1)
    plot(z,-windx/P0,color='r',linewidth=1.5)
    plot(z,-malphax/P0,color='b',linewidth=1.5)
    xlabel('Z/H') ; ylabel(r"$|Fx|$,$|\delta Fx|$")
    xlim(-5,5)

    figure(2)
    plot(z,mreyx/P0-(malphax+mwindx)/P0,color='b',linewidth=1.5)
    plot(z,mreyz/P0-(malphaz+mwindz)/P0,color='r',linewidth=1.5)
    plot([-5,5],[0.,0.],linestyle='--',color='k')
    xlabel('Z/H') ; ylabel(r"$\alpha_R$,$\alpha_Z$")
    plotArrow(-2.3,-0.0045,-0.003,length=7.5e-4,width=0.15,color='r')
    plotArrow(2.6,-0.0045,-0.003,length=7.5e-4,width=0.15,color='r')
    xlim(-5,5)

    figure(3)
    plot(z,sqrt(windx*windx+windz*windz)/P0,color='r',linewidth=1.5)
    plot(z,sqrt(malphax*malphax+malphaz*malphaz)/P0,color='b',linewidth=1.5)
    xlabel('Z/H') ; ylabel(r"$|F|$,$|\delta F|$")
    xlim(-5,5)

    figure(4)
    plot(z,mfwind/P0,color='r',linewidth=1.5)
    plot(z,mfturb/P0,color='b',linewidth=1.5)
    xlabel('Z/H') ; ylabel(r"$F_w, \, F_t$")
    plotArrow(-2.3,0.001,0.003,length=7.5e-4,width=0.15,color='r')
    plotArrow(2.6,0.001,0.003,length=7.5e-4,width=0.15,color='r')
    xlim(-5,5)

    wind_d=0. ; wind_s=0. 
    turb_d=0. ; turb_s=0.
    dz=z[1]-z[0]
    for k in range(nz):
        if (abs(z[k])<zcut):
            turb_d=turb_d+mreyx[k]*dz-(malphax[k]+mwindx[k])*dz
        else:
            turb_s=turb_s+mreyx[k]*dz-(malphax[k]+mwindx[k])*dz
    for k in range(nz-1):
        if ((z[k]<-zcut) and (z[k+1]>-zcut)):
            wind_s=mreyz[k]-malphaz[k]-mwindz[k]
        if ((z[k]<zcut) and (z[k+1]>zcut)):
            wind_s=wind_s-mreyz[k]+malphaz[k]+mwindz[k]
#    print 'Trphi=',Trphi,' (alpha=',Trphi/sqrt(2.*pi)/1.e-6,')'
#    print 'Trphi_w=',Trphiw
    print 'turb=',turb_d,turb_s
    print 'wind=',wind_s
    print 'ratio=',turb_d/wind_s

def plotArrow(xval,ymin,ymax,length=0.1,width=0.2,color='b'):

    arrow(xval,ymin,0.,ymax-ymin,linewidth=1.5,head_length=length,head_width=width,color=color)

def xy_mean(array3d):
    "Compute averages along 1st and 2nd dimension of a 3D array"

    mean_data=np.mean(np.mean(array3d,0),0)
    return mean_data

def yz_mean(array3d):
    "Compute averages along 2nd and 3rd dimension of a 3D array"

    mean_data=np.mean(np.mean(array3d,1),1)
    return mean_data

def z_mean(array3d):
    "Compute averages along 3rd dimension of a 3D array"

    mean_data=np.mean(array3d,2)
    return mean_data

def plt_mean(data0,data55,data41):

    by0 =yz_mean(data0.B[:,:,:,1])
    by55=yz_mean(data55.B[:,:,:,1])
    by41=yz_mean(data41.B[:,:,:,1])

    figure(1)

    plot(data0.x,by0,'r')
    plot(data55.x,by55,'b')
    plot(data41.x,by41,'g')

    xlabel('R') ; ylabel('Mean Bphi')

    plot([2.5,3],[0.045,0.045],'r')
    plot([2.5,3],[0.04,0.04],'b')
    plot([2.5,3],[0.035,0.035],'g')

    text(3.25,0.044,'Initial Bphi')
    text(3.25,0.039,'t=55, w/o resistivity, no flux redistrib')
    text(3.25,0.034,'t=41, resistivity, flux redistrib')

    savefig('mean_bphi.eps')
    
    figure(2)

    plot(data55.x,np.cumsum(by55)/np.max(np.cumsum(by55)),'b')
    plot(data41.x,np.cumsum(by41)/np.max(np.cumsum(by41)),'g')
    xlabel('R') ; ylabel('Cumulative Mean Bphi')

    plot([3,3.5],[0.1,0.1],'b')
    plot([3,3.5],[0.2,0.2],'g')

    text(3.75,0.08,'t=55, w/o resistivity, no flux redistrib')
    text(3.75,0.18,'t=41, resistivity, flux redistrib')

    savefig('cum_mean_bphi.eps')

    return

def plt_mean_II():

    MeanBy(0,0,color='r')
    MeanBy(38,41,color='r')

    plot([4,4.5],[0.01,0.01],'r')
    plot([4,4.5],[0.011,0.011],'b')

    text(4.75,0.0098,'beta=50, t=0, 38<t<41')
    text(4.75,0.0108,'beta=400, t=0, 57<t<60')

    os.chdir('../beta400_320_240_40_piover2/')

    MeanBy(0,0,color='b')
    MeanBy(57,60,color='b')

    os.chdir('../beta50_320_240_40_piover2/')
