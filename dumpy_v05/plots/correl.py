from scipy.fftpack import ifft,fftn
from pylab import *
import numpy as np
import pickle

from dumpy_v05.data import rd_dumses
from dumpy_v05.utils import utils

def correl2d(istart,iend,type='rho',periodic=False):

    # Read data
    print "Reading data..."
    data=rd_dumses.DumsesData(istart)

    # Get dims
    nx,ny,nz=shape(data.rho)

    # y-array
    y=data.y ; Ly=4.

    # Compute correlation function
    print "Computing correlation function..."
    xi=GetXi2d(data.B,data.x,data.time,type=type,periodic=periodic)
    Xitofile(xi,data.x,data.y,data.time,'toto.bin')
    #if (iend>istart):
    #    ifile=range(istart+1,iend+1)
    #    for i in range(len(ifile)):
    #        array,x,time=GetData(ifile[i],dims,type)
    #        xi_i=GetXi(array,x,time,type=type,periodic=periodic)
    #        filename='../correl/'+type+'/'+type+utils.str_suffix(ifile[i],length=2)+'.bin'
    #        Xitofile(xi_i,x,y,time,filename)
    #        xi=xi+xi_i
    #    xi=xi/(len(ifile)+1)

    # Average for 1.5H<|Z|<2.5H
    zmin=1.5 ; kmin=where(data.z>zmin)[0][0]
    zmax=2.5 ; kmax=where(data.z>zmax)[0][0]
    xiCorona=np.zeros((nx,ny))
    for k in range(kmin,kmax+1):
        xiCorona[:,:]=xiCorona[:,:]+xi[:,:,k]
    ktot=kmax-kmin+1
    zmin=-2.5 ; kmin=where(data.z>zmin)[0][0]
    zmax=-1.5 ; kmax=where(data.z>zmax)[0][0]
    for k in range(kmin,kmax+1):
        xiCorona[:,:]=xiCorona[:,:]+xi[:,:,k]
    ktot=ktot+kmax-kmin+1
    xiCorona=xiCorona/ktot
    # Average for |Z|<0.5H
    zmin=-0.5 ; kmin=where(data.z>zmin)[0][0]
    zmax= 0.5 ; kmax=where(data.z>zmax)[0][0]
    xiMidplane=np.zeros((nx,ny))
    for k in range(kmin,kmax+1):
        xiMidplane[:,:]=xiMidplane[:,:]+xi[:,:,k]
    xiMidplane=xiMidplane/(kmax-kmin+1)

    maxXiCorona=np.max(xiCorona)
    maxXiMidplane=np.max(xiMidplane)

    #Get angle
    pi=2.*np.arcsin(1.)
    R=0.10 ; i1=get_angle(xiCorona[:,:],data.x,data.y,R); print 'i(R=0.10)=',i1*180./pi
    R=0.25 ; i2=get_angle(xiCorona[:,:],data.x,data.y,R); print 'i(R=0.25)=',i2*180./pi
    R=0.45 ; i3=get_angle(xiCorona[:,:],data.x,data.y,R); print 'i(R=0.45)=',i3*180./pi
    i=(i1+i2+i3)/3.
    print '   -> i=',i*180./pi

    #Get xi along major and minor axis
    xMajCorona,xiMajCorona=GetAxis(i,data.x,data.y,xiCorona[:,:],major=True)
    xMinCorona,xiMinCorona=GetAxis(i,data.x,data.y,xiCorona[:,:],major=False)
    xMajMidplane,xiMajMidplane=GetAxis(i,data.x,data.y,xiMidplane[:,:],major=True)
    xMinMidplane,xiMinMidplane=GetAxis(i,data.x,data.y,xiMidplane[:,:],major=False)

    #Plot data
    print "Plotting results..."

    vmin=min(xiCorona.min(),xiMidplane.min())
    vmax=1.

    rcParams['image.origin']='lower'
    xmin=min(data.x) ; xmax=max(data.x)
    ymin=min(data.y)-Ly/2. ; ymax=max(data.y)-Ly/2.
    xlabel('Delta_x/H') ; ylabel('Delta_y/H')
    imshow(transpose(xiCorona[:,:])/maxXiCorona,extent=(xmin,xmax,ymin,ymax),vmin=vmin,vmax=1.)
    colorbar()

    figure(2)
    xlabel('Delta_x/H') ; ylabel('Delta_y/H')
    imshow(transpose(xiMidplane[:,:])/maxXiMidplane,extent=(xmin,xmax,ymin,ymax),vmin=vmin,vmax=1.)
    colorbar()

    figure(3)
    xlabel('Delta_x/H') ; ylabel('xi')
    #semilogy(x[nx/2:nx],xi[nx/2-1,ny/2-1,nz/2:nz]/max(xi[nx/2-1,ny/2-1,nz/2:nz]))
    semilogy(xMajCorona,xiMajCorona/max(xiMajCorona))
    semilogy(xMinCorona,xiMinCorona/max(xiMinCorona))
    semilogy(xMajMidplane,xiMajMidplane/max(xiMajMidplane),linestyle='--')
    semilogy(xMinMidplane,xiMinMidplane/max(xiMinMidplane),linestyle='--')
    #semilogy(x[nx/2:nx],exp(-x[nx/2:nx]/0.05),linestyle='--')
    #semilogy(x[nx/2:nx],exp(-x[nx/2:nx]/0.45),linestyle='--')
    xlim((0.,1.))
    ylim((1.e-2,2.))

    return

def DataCorrel(istart,iend,type='v',save=0):
    # Get correlation function from previously calculated data...

    ifile=range(istart,iend+1)
    for i in range(len(ifile)):
        filename=type+utils.str_suffix(ifile[i],length=2)+'.bin'
        print 'Reading data from file',filename,'...'
        xi_loc,x,y,time=Xifromfile(filename)
        if ifile[i]==istart:
            xi=xi_loc
        else:
            xi=xi+xi_loc
    xi=xi/len(ifile)
    maxxi=np.max(xi)
    nx,ny,nz=shape(xi)

    #Get angle
    pi=2.*np.arcsin(1.)
    R=0.10 ; i1=get_angle(xi[:,:,nz/2-1],x,y,R); print 'i(R=0.10)=',i1*180./pi
    R=0.25 ; i2=get_angle(xi[:,:,nz/2-1],x,y,R); print 'i(R=0.25)=',i2*180./pi
    R=0.45 ; i3=get_angle(xi[:,:,nz/2-1],x,y,R); print 'i(R=0.45)=',i3*180./pi
    i=(i1+i2+i3)/3.
    print '   -> i=',i*180./pi

    #Get xi along major and minor axis
    d_maj,xi_maj=GetAxis(i,x,y,xi[:,:,nz/2-1],major=True)
    d_min,xi_min=GetAxis(i,x,y,xi[:,:,nz/2-1],major=False)

    #Plot data
    figure(1,figsize=(5,10))
    rcParams['image.origin']='lower'
    xmin=min(x) ; xmax=max(x) ; Ly=2.*np.arcsin(1.0)
    ymin=min(y)-Ly/2. ; ymax=max(y)-Ly/2.
    xlabel('$\Delta x/H$') ; ylabel('$\Delta y/H$')
    imshow(transpose(xi[:,:,nz/2-1])/maxxi,extent=(xmin,xmax,ymin,ymax))
    colorbar()
    if (save==1):
        savefig(type+'_corr.eps')    

    figure(2)
    xlabel('$\Delta x/H$') ; ylabel('$xi$')
    semilogy(x[nx/2:nx],xi[nx/2-1,ny/2-1,nz/2:nz]/max(xi[nx/2-1,ny/2-1,nz/2:nz]))
    semilogy(d_maj,xi_maj/max(xi_maj))
    semilogy(d_min,xi_min/max(xi_min))
    semilogy(x[nx/2:nx],1.2*exp(-x[nx/2:nx]/0.08),linestyle='--')
    semilogy(x[nx/2:nx],1.2*exp(-x[nx/2:nx]/0.45),linestyle='--')
    xlim((0.,1.))
    ylim((1.e-2,2.))
        

def GetData(i,dims=(128,192,128),type='rho'):
    # Get data from datafile...

    nx,ny,nz=dims

    filename='save_'+utils.str_suffix(i,length=2)+'_f'
    print 'Reading file:',filename,'...'
    cube,x,time=rd_data.get_cube(filename,nx,ny,nz)
    print 'File',filename,'read.'
    print

    if (type=='rho'):
        array=cube[:,:,:,0]
    if (type=='v'):
        array=cube[:,:,:,1:4]
    if (type=='B'):
        array=cube[:,:,:,4:7]

    return array,x,time

def GetXi2d(array,x,time,type,periodic):
    #Return correlation function

    ciso2=1.e-6

    # Get dims
    if rank(array)==3:
        nx,ny,nz=shape(array)
    if rank(array)==4:
        nx,ny,nz,ndim=shape(array)

    if ((type=='rho') or (type=='B')):
        rm_shear=0
    else:
        rm_shear=1

    # y-array
    Ly=2.*np.arcsin(1.0) ; dy=Ly/ny
    y=dy/2.+arange(ny)*dy

    # Perform FFTs
    if periodic:
        if rank(array)==3:
            fft_array=fftn(array)
        if rank(array)==4:
            fft_array=array
            fft_array[:,:,:,0]=fftn(array[:,:,:,0])
            fft_array[:,:,:,1]=fftn(array[:,:,:,1])
            fft_array[:,:,:,2]=fftn(array[:,:,:,2])
    else:
        fft_array=shear_fft2d(array,x,time,type=type,rm_shear=rm_shear)

    # Compute correlation function
    if (type=='rho'):
        xi=fft_array
        for k in range(nz):
            xi[:,:,k]=abs(ifftn(fft_array[:,:,k]*conjugate(fft_array[:,:,k])))/nx/ny
    else:
        xi=fft_array[:,:,:,0]
        for k in range(nz):
            xi[:,:,k]=          abs(ifftn(fft_array[:,:,k,1]*conjugate(fft_array[:,:,k,1])))/nx/ny
#            xi[:,:,k]=          abs(ifftn(fft_array[:,:,k,0]*conjugate(fft_array[:,:,k,0])))/nx/ny
#            xi[:,:,k]=xi[:,:,k]+abs(ifftn(fft_array[:,:,k,1]*conjugate(fft_array[:,:,k,1])))/nx/ny
#            xi[:,:,k]=xi[:,:,k]+abs(ifftn(fft_array[:,:,k,2]*conjugate(fft_array[:,:,k,2])))/nx/ny
            xi[:,:,k]=xi[:,:,k]/ciso2
    xi=np.roll(np.roll(xi,nx/2,0),ny/2,1)
    if not(periodic):
        xi=unshear(xi,x,time,type='rho',direct=-1.e0)

    return xi

def Xitofile(array,x,y,time,filename):

    f=file(filename,'wb')
    pickle.dump(time,f)
    pickle.dump(array,f)
    pickle.dump(x,f)
    pickle.dump(y,f)
    f.close()

    return

def Xifromfile(filename):

    f=file(filename,'rb')
    time =pickle.load(f)
    array=pickle.load(f)
    x    =pickle.load(f)
    y    =pickle.load(f)
    f.close()

    return array,x,y,time

def shear_fft2d(array,x,time,type='rho',rm_shear=0):
    #array of size (nx,ny,nz,3) if type='v'
    #array of size (nx,ny,nz) if type='rho'

    q=1.5e0 ; Omega=1.e-3
    Lx=1. ; Ly=2.*np.arcsin(1.e0) ; Lz=1.
    twopi=4.*np.arcsin(1.e0)

    if rank(array)==3:
        nx,ny,nz=shape(array)
    if rank(array)==4:
        nx,ny,nz,ndim=shape(array)

    # Remove background velocity shear if needed...
    if (rm_shear==1):
        array[:,:,:,1]=remove_shear(array[:,:,:,1],x)

    # Unshear data in real space...
    array=unshear(array,x,time,type=type)

    # Compute FFT...
    fft_array=array
    if (type=='rho'):
        for k in range(nz):
            fft_array[:,:,k]=fftn(array[:,:,k]-np.mean(array[:,:,k]))
    else:
        for k in range(nz):
            fft_array[:,:,k,0]=fftn(array[:,:,k,0])
            fft_array[:,:,k,1]=fftn(array[:,:,k,1])
            fft_array[:,:,k,2]=fftn(array[:,:,k,2])

    return fft_array
    

def unshear(v,x,time,type='v',L=(1.,3.1415,1.),direct=1.e0):

    if rank(v)==3:
        xdim,ydim,zdim=shape(v)
    if rank(v)==4:
        xdim,ydim,zdim,ndim=shape(v)

    q=1.5e0 ; Omega=1.e-3
    dx=L[0]/xdim ; dy=L[1]/ydim
    
    tn=round(time/(L[1]/(q*Omega*L[0])))*L[1]/(q*Omega*L[0])

    if ((type=='v') or (type=='B')):
        for i in range(xdim):
            xx=x[i]-0.5*dx
            jshiftreal=direct*q*Omega*xx*(time-tn)/dy
            jshift    =int(floor(jshiftreal))
            jshiftp1  =int(ceil (jshiftreal))
            vshift  =np.roll(v[i,:,:,0],jshift  ,0)
            vshiftp1=np.roll(v[i,:,:,0],jshiftp1,0)
            v[i,:,:,0]= vshift + (jshiftreal-jshift)*(vshiftp1-vshift)

            xx=x[i]
            jshiftreal=direct*q*Omega*xx*(time-tn)/dy
            jshift    =int(floor(jshiftreal))
            jshiftp1  =int(ceil (jshiftreal))
            vshift  =np.roll(v[i,:,:,1:2],jshift  ,0)
            vshiftp1=np.roll(v[i,:,:,1:2],jshiftp1,0)
            v[i,:,:,1:2]= vshift + (jshiftreal-jshift)*(vshiftp1-vshift)

    if (type=='rho'):
        for i in range(xdim):
            xx=x[i]
            jshiftreal=direct*q*Omega*xx*(time-tn)/dy
            jshift    =int(floor(jshiftreal))
            jshiftp1  =int(ceil (jshiftreal))
            vshift  =np.roll(v[i,:,:],jshift  ,0)
            vshiftp1=np.roll(v[i,:,:],jshiftp1,0)
            v[i,:,:]= vshift + (jshiftreal-jshift)*(vshiftp1-vshift)

    return v

def remove_shear(vy,x):

    xdim,ydim,zdim=shape(vy)
    dx=x[1]-x[0]
    q=1.5e0 ; Omega=1.e-3
    
    for i in range(xdim):
        xx=x[i]
        vy[i,:,:]=vy[i,:,:]+q*Omega*xx
    
    return vy

def get_angle(array,x,y,radius):
    # Return tilt angle of the correlatio function in radians

    N=1300 ; dx=radius/N
    xcirc=-radius+arange(N)*dx
    ymid=0.5*(min(y)+max(y)) ; ycirc=ymid+(radius**2-xcirc**2)**0.5
    f=np.zeros(N)

    for ix in range(N):
        f[ix]=GetInterpolate(xcirc[ix],ycirc[ix],x,y,array)

    pos=[i for i,x in enumerate(f) if x==max(f)][0]

    angle=np.arctan(abs(xcirc[pos])/(ycirc[pos]-ymid))

    return angle

def GetInterpolate(x0,y0,x,y,array,order=2):
    #Return 2d interpolation of array at position x0,y0

    xmin=x[0] ; dx=x[1]-x[0]
    ymin=y[0] ; dy=y[1]-y[0]
    nx=shape(x)[0] ; ny=shape(y)[0]

    if order==1:
        i=int(round((x0-xmin)/dx)) ; j=int(round((y0-ymin)/dy))
        interpolate_array=array[i,j]

    if order==2:
        i0=floor((x0-xmin)/dx) ; i1=min(ceil((x0-xmin)/dx)  ,nx-1)
        j0=floor((y0-ymin)/dy) ; j1=min(ceil((y0-ymin)/dy),ny-1)
        if (i0==-1):
            t=(x0-xmin+dx)/dx ; i0=0
        else:
            t=(x0-x[i0])/dx
        if (j0==-1):
            u=(y0-ymin+dy)/dy ; j0=0
        else:
            u=(y0-y[j0])/dy

        interpolate_array = (1-t)*(1-u)*array[i0,j0] + \
                                t*(1-u)*array[i1,j0] + \
                            (1-t)*u    *array[i0,j1] + \
                                t*u    *array[i1,j1]

    return interpolate_array

def GetAxis(i,x,y,array,major=True):

    N=200
    x0=0.5*(min(x)+max(x)) ; y0=0.5*(min(y)+max(y))

    x_axe=np.zeros(N) ; y_axe=np.zeros(N)
    d_axe=np.zeros(N) ; array_axe=np.zeros(N)

    if major:
        dx=(x0-min(x))/N
        x_axe=x0-arange(N)*dx
        y_axe=y0+(x0-x_axe)/np.tan(i)
    else:
        dx=(max(x)-x0)/N
        x_axe=x0+arange(N)*dx
        y_axe=y0+(x_axe-x0)*np.tan(i)

    d_axe=((x_axe-x0)**2+(y_axe-y0)**2)**0.5
    for pos in range(N):
        if (x_axe[pos]>min(x) and x_axe[pos]<max(x) and y_axe[pos]<max(y)):
            array_axe[pos]=GetInterpolate(x_axe[pos],y_axe[pos],x,y,array)

    return d_axe,array_axe

def double_size(array,time):

    xdim,ydim,zdim=shape(array)
    xdim=xdim-2 ; ydim=ydim-2 ; zdim=zdim-2

    q=1.5e0 ; Omega=1.e-3 ; Lx=1. ; Ly=2.*np.arcsin(1.e0)
    dy=Ly/ydim ; dx=Lx/xdim ; xmin=-0.5 ; xmax=+0.5

    dbl_array=np.zeros((2*xdim,2*ydim,2*zdim))

    dbl_array[:,:,:]=1.e0

    dbl_array[xdim/2:3*xdim/2,ydim/2:3*ydim/2,zdim/2:3*zdim/2]=\
        array[1:xdim+1,1:ydim+1,1:zdim+1]

    #y-bound
    dbl_array[xdim/2:3*xdim/2,3*ydim/2:2*ydim,zdim/2:3*zdim/2]=\
        dbl_array[xdim/2:3*xdim/2,ydim/2:ydim,zdim/2:3*zdim/2]
    dbl_array[xdim/2:3*xdim/2,0:ydim/2,zdim/2:3*zdim/2]=\
        dbl_array[xdim/2:3*xdim/2,ydim:3*ydim/2,zdim/2:3*zdim/2]

    #z-bound
    dbl_array[xdim/2:3*xdim/2,:,3*zdim/2:2*zdim]=\
        dbl_array[xdim/2:3*xdim/2,:,zdim/2:zdim]
    dbl_array[xdim/2:3*xdim/2,:,0:zdim/2]=\
        dbl_array[xdim/2:3*xdim/2,:,zdim:3*zdim/2]

    #x-bound
    tn=round(time/(Ly/(q*Omega*Lx)))*Ly/(q*Omega*Lx)
    for i in range(xdim/2):
        xx=Lx
        jshiftreal=q*Omega*xx*(time-tn)/dy
        jshift    =int(floor(jshiftreal))
        jshiftp1  =int(ceil (jshiftreal))
        vshift  =np.roll(dbl_array[xdim+i,:,:],jshift  ,0)
        vshiftp1=np.roll(dbl_array[xdim+i,:,:],jshiftp1,0)
        dbl_array[i,:,:]= vshift + (jshiftreal-jshift)*(vshiftp1-vshift)
        jshift    =int(floor(-jshiftreal))
        jshiftp1  =int(ceil (-jshiftreal))
        vshift  =np.roll(dbl_array[xdim/2+i,:,:],jshift  ,0)
        vshiftp1=np.roll(dbl_array[xdim/2+i,:,:],jshiftp1,0)
        dbl_array[3*xdim/2+i,:,:]= vshift + (-jshiftreal-jshift)*(vshiftp1-vshift)

    return dbl_array

#    #Compute wavenumbers...
#    wave_nb_x=np.zeros(nx)
#    wave_nb_y=np.zeros(ny)
#    wave_nb_z=np.zeros(nz)
#
#    for i in range(nx/2+1):
#        wave_nb_x[i]=i
#        if (i<nx/2-1):
#            wave_nb_x[i+nx/2+1]=-nx/2+1+i
#    wave_nb_x=wave_nb_x*twopi/Lx
#    for j in range(ny/2+1):
#        wave_nb_y[j]=j
#        if (j<ny/2-1):
#            wave_nb_y[j+ny/2+1]=-ny/2+1+j
#    wave_nb_y=wave_nb_y*twopi/Ly
#    for k in range(nz/2+1):
#        wave_nb_y[k]=k
#        if (k<nz/2-1):
##            wave_nb_z[k+nz/2+1]=-nz/2+1+k
#    wave_nb_z=wave_nb_z*twopi/Lz
#
#    # Unshear FFTs in k-space...
#    array_unshear=fft_array
#    Tshift=2*Ly/(q*Omega*Lx)
#    time_mod_Tshift=time-int(time/Tshift)*Tshift
#    dwave=wave_nb_x[1]-wave_nb_x[0]
#    for j in range(ny):
#        kshift_x=q*Omega*time_mod_Tshift*wave_nb_y[j] 
#        for i in range(nx):
#            k_lag_x=wave_nb_x[i]-kshift_x
#            if (k_lag_x>max(wave_nb_x)):
#                k_lag_x=k_lag_x-twopi/Lx*nx
#            if (k_lag_x<min(wave_nb_x)):
#                k_lag_x=k_lag_x+twopi/Lx*nx
#            
#            if (k_lag_x/twopi*Lx>=0):
#                i_lag_x=int(k_lag_x/twopi*Lx)
#            else:
#                i_lag_x=int(k_lag_x/twopi*Lx)-1
#
#            if (i_lag_x<0):
#                i_lag_x=i_lag_x+nx
#            i_lag_x_p1=i_lag_x+1
#            if (i_lag_x_p1==nx):
#                i_lag_x_p1=0
#
#            for k in range(nz):
#                darray=fft_array[i_lag_x_p1,j,k]-fft_array[i_lag_x,j,k]
#                dk    =k_lag_x-wave_nb_x[i_lag_x]
#                #array_unshear[i,j,k]=fft_array[i_lag_x,j,k]+darray*dk/dwave


