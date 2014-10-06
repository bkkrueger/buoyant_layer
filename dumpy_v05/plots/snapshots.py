import time
from matplotlib.ticker import MaxNLocator
from pylab import *

from dumpy_v05.data import rd_dumses
from dumpy_v05.data import rd_movie

def plot1d(idump=1,nvar=8,type='rho',linestyle='-',color='g',thick=1.0,save=0,scale=1,bin=False):
    data =rd_dumses.DumsesData(idump,nvar,bin=bin)
    data0=rd_dumses.DumsesData(0    ,nvar,bin=bin)
    data_y =data.get_1d_y(type)
    data_y0=data0.get_1d_y(type)
    if (type=='bx'):
        data_by=data.get_1d_y('by')
    else:
        data_by=1
    print 'Scaling factor since t=0: ',max(abs(data_y))/max(abs(data_y0))
    yname=type
    if (type=='bx'):
        yname='$B_x/B_z$'
    params = {'axes.labelsize': 20,
              'xtick.labelsize': 16,
              'ytick.labelsize': 16,
              'text.usetex': True}
    rcParams.update(params)
    xlabel('$z/H$') ; ylabel(yname)
    #if (type=='rho'):
    #    plot(data.y,data_y/data_y0-1,linewidth=thick,linestyle=linestyle)
    #    xlim((-2.,2.)) ; ylim((-1.,1.))
    #else:
    plot(data.y,scale*data_y/data_by,linewidth=thick,linestyle=linestyle,color=color)
    #xlim((min(data.y),max(data.y)))
    xlim((-2.,2.))
    if save==1:
        savefig(type+'.eps')

def plot2d(data,nvar=8,type='rho',slice='xy',thick=1.0,save=0):
    "Extract 2D slice of the data object before plotting"

    fig = figure()
    ax = fig.add_subplot(111)
    ax.xaxis.set_major_locator(MaxNLocator(4))

    if slice=='xy':
        if type=='beta':
            data_xy =data.get_2d_xy('bz')
            data_xy2=data.rho[:,:,2]#data.get_2d_xy('rho')
            data_xy=data_xy**2/(data_xy2*1.e-3**2*2.)
        else:
            data_xy=data.get_2d_xy(type)
    elif slice=='xz':
        data_xy=data.get_2d_xz(type)
    elif slice=='yz':
        if type=='alfven':
            data_yz =data.get_2d_yz('by')**2/2./data.get_2d_yz('rho')
            data_yz[:,:]=log(data_yz[:,:]/1.e-6)
            indices=where(data_yz<-1.)
            data_yz[indices]=-1.
        else:
            data_yz=data.get_2d_yz(type)
    xmin=min(data.x) ; xmax=max(data.x)
    ymin=min(data.y) ; ymax=max(data.y)
    zmin=min(data.z) ; zmax=max(data.z)

#    data_xy=data.rho[:,:,2]
#    print 'fuck',data_xy.min()

    #Plot data
    rcParams['image.origin']='lower'
    matplotlib.rcParams['font.size']=14
    hot()
    #xlabel('phi') ; ylabel('R')
    #imshow(data_xy,aspect='auto',extent=(ymin,ymax,xmin,xmax))
    if (slice=='xy'):
        xlabel('x/H') ; ylabel('y/H')
        imshow(transpose(data_xy),aspect='auto',extent=(xmin,xmax,ymin,ymax))
    if (slice=='xz'):
        xlabel('x/H') ; ylabel('z/H')
        imshow(transpose(data_xy),extent=(xmin,xmax,zmin,zmax))
    if (slice=='yz'):
        xlabel('y/H') ; ylabel('z/H')
        print data_yz.max()
        imshow(transpose(data_yz),extent=(ymin,ymax,zmin,zmax),vmax=3.43)
        xlim(int(round(ymin)),int(round(ymax)))
    colorbar()

    
    if save==1:
        savefig(type+'.eps')

def vortensity(data):

    r=data.x ; phi=data.y
    dr=r[1]-r[0] ; dphi=phi[1]-phi[0]
    vr  =data.rhou[:,:,0,0]/data.rho[:,:,0]
    vphi=data.rhou[:,:,0,1]/data.rho[:,:,0]

    xmin=min(r)   ; xmax=max(r)
    ymin=min(phi) ; ymax=max(phi)

    (nx,ny)=vr.shape

    vortensity=np.zeros((nx,ny))

    for j in range(1,ny-1):
        for i in range(1,nx-2):
            rc_i=r[i] ; rc_ip1=r[i+1] ; rc_im1=r[i-1]
            vortensity[i,j]=                (rc_ip1*vphi[i+1,j]-rc_im1*vphi[i-1,j])/dr/rc_i/2.
            vortensity[i,j]=vortensity[i,j]-(vr[i,j+1]-vr[i,j-1])/dphi/rc_i/2.
            #vortensity[i,j]=vortensity[i,j]/data.rho[i,j,0]

    #imshow(vortensity[1:nx-1,1:ny-1],aspect='auto',extent=(ymin,ymax,xmin,xmax))
    #ylabel('radius')
    #xlabel('phi')
    #colorbar()

    init_vort=init_vorticity(data.x)
    for j in range(1,ny-1):
        for i in range(1,nx-2):
            vortensity[i,j]=(vortensity[i,j]-init_vort[i])/init_vort[i]
            vortensity[i,j]=max(min(vortensity[i,j],2.),-2.)

    return vortensity

def movie_vortensity(start,stop):

    #hot()
    spectral()

    for idump in range(start,stop+1):
        print idump
        data=rd_dumses.DumsesData(idump)
        vdata=vortensity(data)
        imshow(vdata,aspect='auto',animated=True)
        colorbar()
        draw()
        savefig('png/snap'+str(idump)+'.png')
        clf()

def movie(start,stop):

    hot()
    #spectral()

    ifile=0
    for idump in range(start,stop+1):
        print idump
        data=rd_movie.DumsesMovie(idump)
        rho=data.get_2d_xy()
        imshow(rho,aspect='auto',animated=True)
        colorbar()
        draw()
        for iloop in range(2):
            savefig('png/rho'+str(ifile)+'.png')
            ifile=ifile+1
        clf()

def init_vorticity(r):

    return 0.5*r**(-0.5)


def GetAz(data):

    bx=data.get_2d_xy('bx')
    by=data.get_2d_xy('by')

    (nx,ny)=bx.shape
    dx=data.x[1]-data.x[0]
    dy=data.y[1]-data.y[0]
    xmin=min(data.x)-dx/2. ; xmax=max(data.x)+dx/2.
    ymin=min(data.y)-dy/2. ; ymax=max(data.y)+dy/2.

    Az=np.zeros((nx,ny))
    for j in range(ny):
        for i in range(1,nx):
            Az[i,j]=Az[i-1,j]-by[i-1,j]*dx
        if (j<ny-1):
            Az[0,j+1]=Az[0,j]+bx[0,j]*dy
    azmin=Az.min()
    Az[0:nx,0:ny]=Az[0:nx,0:ny]-azmin

    N=21 ; limit=0.#5e-7
    V=np.zeros(N) ; eps=(Az.max()-Az.min())/N/10.
    #V=Az.min()+arange(N)*(Az.max()-Az.min())/N
    V=arange(N)*Az.max()/N+limit
    #print V,Az.min()
    #print Az

    return Az,V,xmin,xmax,ymin,ymax

def ContAnimate(start,stop,save=0):

    rcParams['image.origin']='lower'

    #clf()
    #figure(figsize=(3,8))
    for ifile in range(start,stop):
        print ifile
        data=rd_dumses.DumsesData(ifile)
        Az,V,xmin,xmax,ymin,ymax=GetAz(data)

        if (ifile==start):
            figure(figsize=(xmax-xmin,ymax-ymin))
        contour(data.x,data.y,transpose(Az),V[1:],colors='black',figsize=(xmax-xmin,ymax-ymin),aspect='auto')
#        imshow(transpose(Az),extent=(xmin,xmax,zmin,zmax))#,aspect='auto')
        #clabel(cs)
        draw()
        
        if (save==1):
            savefig('png/snap'+str(ifile)+'.png')
        
        time.sleep(0.0)
        if (ifile<stop-1):
            clf()


