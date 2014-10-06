from pylab import *
import numpy as np

from dumpy_v05.utils.utils import *
from dumpy_v05.data import rd_zprof
from dumpy_v05.plots import cut

def plot1D(idump=1,type='rho',thick=1.0,save=0,log=0):
    filename='zprof/zprof_'+str_suffix(idump)+'.h5'
    data=rd_zprof.zData(filename)
    field=data.get_array(type)
    xlabel('z/H') ; ylabel(type)
    if (log==1):
        semilogy(data.z,field,linewidth=thick)
    else:
        plot(data.z,field,linewidth=thick)
    if save==1:
        savefig(type+'.eps')

    return

def spacetime(nframe=10,save=0,type='rho',norm=1.,withBetaOne=False):

    rc('text', usetex=True)
    rcParams['image.origin']='lower'
    c0=1.e-3 ; c0sq=c0*c0
    time=np.zeros(nframe)
    Torb=4.*np.arcsin(1.e0)/1.e-3

    for iframe in range(nframe-1):
        filename='zprof/zprof_'+str_suffix(iframe+1)+'.h5'
        data=rd_zprof.zData(filename)
        if (type=='by'):
            field=data.get_array('by')/data.get_array('bz')
        elif (type=='bx'):
            field=data.get_array('bx')/data.get_array('bz')
        elif (type=='vx'):
            field=data.get_array('vx')#/data.get_array('vz')
        elif (type=='massFlux'):
            field=data.get_array('vz')*data.get_array('rho')
        elif (type=='pmagNorm'):
            Pmag=(data.get_array('bx')**2+data.get_array('by')**2)
#            field=np.log10(Pmag/data.get_array('rho')/c0sq)
#            field=np.log10(Pmag/data.get_array('rho')/c0sq)
            field=np.log10(Pmag/data.get_array('bz')**2)
#            field=Pmag/data.get_array('bz')**2
        else:
            field=data.get_array(type)
        if iframe==0: 
            nz=shape(field)[0]
            image=np.zeros((nz,nframe))
            z=np.zeros(nz) ; z=data.z
        if (type=='rho'):
            image[:,iframe]=np.log10(field)
        else:
            image[:,iframe]=field
        time[iframe]=data.time
    time=time/Torb

    xlabel('Time (in orbits)') ; ylabel('z/H')
    imshow(image*norm,aspect='auto',extent=(min(time),max(time),min(z),max(z)))

    if withBetaOne:
        timeDump,zBottom,zTop=cut.readBetaOne('betaOne.pickle')
        timeDumpNorm=array(timeDump)/Torb
        plot_curve(timeDumpNorm,zBottom,thick=1.5)
        plot_curve(timeDumpNorm,zTop   ,xtitle='Time (in orbits)',ytitle='z/H',thick=1.5)
        xlim(min(time),max(time))
        ylim(min(z),max(z))

    if save==1:
        savefig('spacetime.eps')

    return

def Mean1d(istart,istop,type='rho',dir='zprof',title='rho',save=0,norm=1,log=False,thick=1.0,style='-',color='black',fold=0):

    rc('text', usetex=True)
    c0=1.e-3 ; c0sq=c0*c0

    for iframe in range(istart,istop+1):
        filename=dir+'/zprof_'+str_suffix(iframe+1)+'.h5'
        data=rd_zprof.zData(filename)
        if (type=='bx'):
            field=data.get_array('bx')#/data.get_array('bz')
        if (type=='bxoby'):
            field=np.abs(data.get_array('bx')/data.get_array('by'))
        elif (type=='vax'):
            field=data.get_array('bx')/data.get_array('rho')
        elif (type=='vay'):
            field=data.get_array('by')/data.get_array('rho')
        elif (type=='alpha'):
            field=data.get_array('maxwell')+data.get_array('reynolds')
            norm=1.e6
        elif (type=='beta'):
            Pmag=0.5*(data.get_array('bx')**2+data.get_array('by')**2+data.get_array('bz')**2)
            field=data.get_array('rho')*c0sq/Pmag
        elif (type=='pmagNorm'):
            Pmag=0.5*(data.get_array('bx')**2+data.get_array('by')**2)
            field=Pmag/data.get_array('rho')/c0sq
        elif (type=='angle'):
#            field=np.arctan(abs(data.get_array('bx')/data.get_array('bz')))*180./pi
            field=np.arctan(data.get_array('bx')/data.get_array('bz'))*180./pi
        elif (type=='angleVel'):
#            field=np.arctan(abs(data.get_array('vx')/data.get_array('vz')))*180./pi
            field=np.arctan(data.get_array('vx')/data.get_array('vz'))*180./pi
        else:
            field=data.get_array(type)
        if iframe==istart: 
            nz=shape(field)[0]
            image=np.zeros(nz)
            z=np.zeros(nz) ; z=data.z
            image=field
        else:
            image=image+field
    image=image/(istop-istart+1)

    if (fold<>0):
        halfZ=z[nz/2:]
        halfImage=foldArray(image,fold)
        plot_curve(halfZ,halfImage*norm,xtitle=r'z/H',ytitle=title,log=log,thick=thick,style=style,color=color)
        print halfImage[nz/2-1]*norm
    else:
        plot_curve(z,image*norm,xtitle=r'z/H',ytitle=title,log=log,thick=thick,style=style,color=color)

    if save==1:
        savefig('mean_zprofile.eps')

    rho=data.get_array('rho')
    print sum(image*norm*rho)/sum(rho)

    return image

def invariant(istart,istop,zcut=3.):

    dir='zprof'
    rc('text', usetex=True)
    matplotlib.rcParams['font.size']=14
    c0=1.e-3 ; c0sq=c0*c0 ; q=1.5 ; Omega=c0

    for iframe in range(istart,istop+1):
        filename=dir+'/zprof_'+str_suffix(iframe+1)+'.h5'
        data=rd_zprof.zData(filename)
        if iframe==istart:
            rho=data.get_array('rho')
            bx=data.get_array('bx')
            by=data.get_array('by')
            bz=data.get_array('bz')
            vx=data.get_array('vx')
            vy=data.get_array('vy')
            vz=data.get_array('vz')
            rhovz=data.get_array('rhovz')
            rhovy=data.get_array('rhovy')
            rhovx=data.get_array('rhovx')
            nz=shape(bx)[0]
            image=np.zeros(nz)
            z=np.zeros(nz) ; z=data.z ; dz=z[1]-z[0]
        else:
            rho=rho+data.get_array('rho')
            bx=bx+data.get_array('bx')
            by=by+data.get_array('by')
            vx=vx+data.get_array('vx')
            vy=vy+data.get_array('vy')
            vz=vz+data.get_array('vz')
            rhovx=rhovx+data.get_array('rhovx')
            rhovy=rhovy+data.get_array('rhovy')
            rhovz=rhovz+data.get_array('rhovz')
    rho=rho/(istop-istart+1)
    bx=bx/(istop-istart+1)
    by=by/(istop-istart+1)
    vx=vx/(istop-istart+1)
    vy=vy/(istop-istart+1)
    vz=vz/(istop-istart+1)
    rhovz=rhovz/(istop-istart+1)
    rhovx=rhovx/(istop-istart+1)
    rhovy=rhovy/(istop-istart+1)
    #vx=rhovx/rho ; vz=rhovz/rho ; vy=rhovy/rho

    tanTheta=bx/bz
    alpha=rhovz/bz

    L0=np.zeros(nz)
    Ek=np.zeros(nz) ; Epsi=np.zeros(nz) ; Et=np.zeros(nz) ; Eb=np.zeros(nz)
    vystar=np.zeros(nz) ; xstar=np.zeros(nz)
    istart=np.int(64*(5-zcut))
    for i in range(istart,-1,-1):
        alpha[i]=-0.12
        L0[i]=vy[i]+(2-q)*Omega*np.sum(-tanTheta[i:istart])*dz
        Ek[i]=0.5*(vx[i]**2+(vy[i]-q*Omega*np.sum(-tanTheta[i:istart])*dz)**2+vz[i]**2)
        Epsi[i]=0.5*Omega**2*(z[i]**2-z[istart]**2)-q*Omega**2*(np.sum(-tanTheta[i:istart])*dz)**2
        vystar[i]=vy[i]-q*Omega*np.sum(-tanTheta[i:istart])*dz-alpha[i]*by[i]/rho[i]
#        Eb[i]=-by[i]/alpha[i]*((vy[i]-q*Omega*np.sum(-tanTheta[i:istart])*dz)-q*Omega*np.sum(-tanTheta[istart:i:-1])*dz-alpha[i]/rho[i]*by[i])
        Eb[i]=-by[i]/alpha[i]*vystar[i]
        xstar[i]=np.sum(-tanTheta[i:istart])*dz
        Et[i]=c0**2*np.log(rho[i])
        print i,xstar[i],0.5*Omega**2*(z[i]**2-3.5**2-3.*xstar[i]),Epsi[i]
    istart=np.int((5+zcut)*64)
    for i in range(istart,nz):
        alpha[i]=0.05
        L0[i]=vy[i]+(2-q)*Omega*np.sum(tanTheta[istart:i])*dz
        Ek[i]=0.5*(vx[i]**2+(vy[i]-q*Omega*np.sum(tanTheta[istart:i])*dz)**2+vz[i]**2)
        Epsi[i]=0.5*Omega**2*(z[i]**2-z[istart]**2)-q*Omega**2*(np.sum(tanTheta[istart:i])*dz)**2
        vystar[i]=vy[i]-q*Omega*np.sum(tanTheta[istart:i])*dz-alpha[i]*by[i]/rho[i]
        Eb[i]=-by[i]/alpha[i]*vystar[i]
#        Eb[i]=-by[i]/alpha[i]*((vy[i]-q*Omega*np.sum(tanTheta[istart:i])*dz)-q*Omega*np.sum(tanTheta[istart:i])*dz-alpha[i]/rho[i]*by[i])
#        vystar[i]=vy[i]-alpha[i]*by[i]/rho[i]
        xstar[i]=np.sum(tanTheta[istart:i])*dz
        Et[i]=c0**2*np.log(rho[i])
#        print z[i],xstar[i]

    #plot(z,alpha)
    #stop

    #L0[(z>-2)*(z<2)]=0.
    #by[(z>-2)*(z<2)]=0.


    #plot_curve(z,L0,xtitle='Z/H',color='b')    
    #plot_curve(z,-by/alpha,xtitle='Z/H',color='r')
    #plot_curve(z,L0-by/alpha,color='g')

    figure(1)
    
    plot_curve(z[z<-zcut],L0[z<-zcut],xtitle='z/H',color='b',thick=2)    
    plot_curve(z[z<-zcut],-by[z<-zcut]/alpha[z<-zcut],xtitle='z/H',color='r',thick=2)
    plot_curve(z[z<-zcut],L0[z<-zcut]-by[z<-zcut]/alpha[z<-zcut],color='g',thick=2)

    plot_curve(z[z>+zcut],L0[z>+zcut],xtitle='z/H',color='b',thick=2)    
    plot_curve(z[z>+zcut],-by[z>+zcut]/alpha[z>+zcut],xtitle='z/H',color='r',thick=2)
    plot_curve(z[z>+zcut],L0[z>+zcut]-by[z>+zcut]/alpha[z>+zcut],color='g',thick=2)

    #plot([-zcut,-zcut],[-5.e-3,+5.e-3],color='k',linewidth=1.5)
    #plot([+zcut,+zcut],[-5.e-3,+5.e-3],color='k',linewidth=1.5)
    #plot([-zcut,+zcut],[-5.e-3,+5.e-3],color='k',linewidth=1.5)
    #plot([-zcut,+zcut],[+5.e-3,-5.e-3],color='k',linewidth=1.5)

    xlabel('z/H') ; xlim(-5,5)
    ylabel('') ; ylim(-5.e-3,5.e-3)

    figure(2)
    
    plot_curve(z[z<-zcut],Ek[z<-zcut],color='b',xtitle='z/H',thick=2)
    plot_curve(z[z<-zcut],Et[z<-zcut],color='g',thick=2)
    plot_curve(z[z<-zcut],Epsi[z<-zcut],color='y',thick=2)
    plot_curve(z[z<-zcut],Eb[z<-zcut],color='r',thick=2)
    plot_curve(z[z<-zcut],(Ek[z<-zcut]+Et[z<-zcut]+Epsi[z<-zcut]+Eb[z<-zcut]),xtitle='z/H',color='k',thick=2)
    #plot_curve(z[z<-zcut],(Ek[z<-zcut]+Et[z<-zcut]+Epsi[z<-zcut]),xtitle='Z/H',color='k',thick=2,style='--')

    plot_curve(z[z>+zcut],Ek[z>+zcut],color='b',xtitle='z/H',thick=2)
    plot_curve(z[z>+zcut],Et[z>+zcut],color='g',thick=2)
    plot_curve(z[z>+zcut],Epsi[z>+zcut],color='y',thick=2)
    plot_curve(z[z>+zcut],Eb[z>+zcut],color='r',thick=2)
    plot_curve(z[z>+zcut],(Ek[z>+zcut]+Et[z>+zcut]+Epsi[z>+zcut]+Eb[z>+zcut]),xtitle='z/H',color='k',thick=2)
    #plot_curve(z[z>+zcut],(Ek[z>+zcut]+Et[z>+zcut]+Epsi[z>+zcut]),xtitle='Z/H',color='k',thick=2,style='--')

    ylimit=1.5e-5
    #plot([-zcut,-zcut],[-ylimit,+ylimit],color='k',linewidth=1.5)
    #plot([+zcut,+zcut],[-ylimit,+ylimit],color='k',linewidth=1.5)
    #plot([-zcut,+zcut],[-ylimit,+ylimit],color='k',linewidth=1.5)
    #plot([-zcut,+zcut],[+ylimit,-ylimit],color='k',linewidth=1.5)

    ylim(-ylimit,ylimit) ; xlim(-5.,5.)

    figure(3)
#    plot(z,vx-vz*bx/bz,color='b')
    plot(z,vystar,color='r')
    plot(z,-alpha*by/rho,color='b')
    #plot(z,-q*Omega*xstar,color='g')
    plot(z,vy-q*Omega*xstar,color='g')
#    plot(z,vy**2)#-vz*bx/bz)
#    plot(z,vz**2)#-vz*bx/bz)

    figure(4)
    plot(z,rhovz*1.e3,color='k',linewidth=1.5)
    xlabel('z/H') ; ylabel(r"$\rho v_z/(\rho_0 c_0)$")
    gcf().subplots_adjust(left=0.15)
    xlim(min(z),max(z))

#    plot(z,rho*vx/bx)
#    plot(z,rho*vz/bz,color='b')
#    plot(z,rhovx/bx,color='g')

    figure(5)
    plot(z,alpha)

    figure(6)
    plot(z,z**2)
    plot(z,-3.*xstar**3)
    plot(z,(z**2-3.*xstar**2))


def alfven(istart,istop,title='rho',dir='zprof',save=0,norm=1,log=False,thick=1.0,style='-',color='green'):

    for iframe in range(istart,istop+1):
        filename=dir+'/zprof_'+str_suffix(iframe+1)+'.h5'
        data=rd_zprof.zData(filename)
        #Alfven speed
        c2=1.e-6
        field=sqrt(data.get_array('bz')**2/data.get_array('rho'))/c2**0.5
        if iframe==istart: 
            nz=shape(field)[0]
            image=np.zeros(nz)
            z=np.zeros(nz) ; z=data.z
            image=field
        else:
            image=image+field
    image=image/(istop-istart+1)

    halfImage=foldArray(image,1.)
    plot_curve(z[nz/2:],halfImage*norm,xtitle='Z/H',ytitle=title,log=log,thick=thick,style=style,color=color)

    return

def mhdwave(istart,istop,title='rho',dir='zprof',save=0,norm=1,log=False,thick=1.0,style='-',type=1,color='black'):
    "Compute mhd wave speed in the vertical direction: type=1 is for fast, type=2 for slow"

    for iframe in range(istart,istop+1):
        filename=dir+'/zprof_'+str_suffix(iframe+1)+'.h5'
        data=rd_zprof.zData(filename)
        #Fast magnetosonic speed
        va2=data.get_array('bx')**2+data.get_array('by')**2+data.get_array('bz')**2
        va2=va2/data.get_array('rho')
        vaz2=data.get_array('bz')**2/data.get_array('rho')
        c2=1.e-6
        if type==1:
            field=sqrt((va2+c2)/2+sqrt((va2+c2)**2/4.-c2*vaz2))/c2**0.5
        else:
            field=sqrt((va2+c2)/2-sqrt((va2+c2)**2/4.-c2*vaz2))/c2**0.5
        if iframe==istart: 
            nz=shape(field)[0]
            image=np.zeros(nz)
            z=np.zeros(nz) ; z=data.z
            image=field
        else:
            image=image+field
    image=image/(istop-istart+1)

    halfImage=foldArray(image,1.)
    plot_curve(z[nz/2:],halfImage*norm,xtitle='Z/H',ytitle=title,log=log,thick=thick,style=style,color=color)

    return

def sonicPoints(istart,istop,zmax=5.):
    "Plot vertical profile of vertical gas & waves velocities"

    matplotlib.rcParams['font.size']=16
    Mean1d(istart,istop,type='vz',thick=2,log=True,norm=1.e3,fold=-1,color='b')
    alfven(istart,istop,thick=2,log=True,style='--',color='g')
    mhdwave(istart,istop,type=1,thick=2,log=True,style='--',color='r')
    mhdwave(istart,istop,type=2,thick=2,log=True,style='--',color='c')
    xlim(0.,zmax) ; ylim(1.e-3,1.e1)
    xlabel('z/H') ; ylabel('Vertical velocities')

def foldArray(array,sign):
    "Fold a 1d array in two"

    N=shape(array)[0]
    outArray=zeros(N/2)
    for i in range(N/2):
        outArray[i]=array[N/2+i]+sign*array[N/2-1-i]
    outArray=outArray/2.

    return outArray

def plot_curve(x,y,xtitle='x',ytitle='y',log=False,thick=1.0,style='-',color='black'):
    
    matplotlib.rcParams['font.size']=16
    xlabel(xtitle)
    ylabel(ytitle)
    if (log==False):
        plot(x,y,linewidth=thick,linestyle=style,color=color)
    else:
        semilogy(x,y,linewidth=thick,linestyle=style,color=color)

    return

