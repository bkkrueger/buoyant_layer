from pylab import *
import numpy as np

from dumpy_v05.utils.utils import *
from dumpy_v05.data import rd_yprof

def plot1D(idump=1,type='rho',thick=1.0,save=0,log=0):
    filename='zprof/zprof_'+str_suffix(idump)+'.h5'
    data=rd_yprof.zData(filename)
    field=data.get_array(type)
    xlabel('z/H') ; ylabel(type)
    if (log==1):
        semilogy(data.y,field,linewidth=thick)
    else:
        plot(data.y,field,linewidth=thick)
    if save==1:
        savefig(type+'.eps')

    return

def spacetime(nframe=10,save=0,type='rho'):

    time=np.zeros(nframe)
    Torb=4.*np.arcsin(1.e0)/1.e-3

    for iframe in range(nframe):
        filename='zprof/zprof_'+str_suffix(iframe+1)+'.h5'
        data=rd_yprof.zData(filename)
        field=data.get_array(type)
        if iframe==0: 
            nz=shape(field)[0]
            image=np.zeros((nz,nframe))
            z=np.zeros(nz) ; z=data.y
        if (type=='rho'):
            image[:,iframe]=np.log10(field)
        else:
            image[:,iframe]=field
        time[iframe]=data.time
    time=time/Torb

    xlabel('Time (in orbits)') ; ylabel('z/H')
    imshow(image,aspect='auto',extent=(min(time),max(time),min(z),max(z)))

    if save==1:
        savefig('spacetime.eps')

    return

def Mean1d(istart,istop,type='rho',save=0):

    for iframe in range(istart,istop+1):
        filename='zprof/zprof_'+str_suffix(iframe+1)+'.h5'
        data=rd_yprof.zData(filename)
        field=(data.get_array(type)**2/2.)/(data.get_array('rho')*1.e-6)#**0.5#/1.e-3
        field=(data.get_array(type)**2/2.)#/(data.get_array('rho')*1.e-6)#**0.5#/1.e-3
        field=data.get_array(type)#/(data.get_array('rho')*1.e-6)#**0.5#/1.e-3
        field=data.get_array(type)/data.get_array('by')
#        field=data.get_array('bz')/data.get_array('rho')**0.5#/1.e-3
#        field=2.*data.get_array('rho')*1.e-6/(data.get_array('bx')**2+data.get_array('by')**2+data.get_array('bz')**2)
#        field=2.*data.get_array('rho')*1.e-6/(data.get_array('by')**2)
#        field=data.get_array('rho')*data.get_array('vy')/1.e-3
        if iframe==istart: 
            nz=shape(field)[0]
            image=np.zeros(nz)
            z=np.zeros(nz) ; z=data.y
            image=field
        image=image+field
    image=image/(istop-istart+1)

    print nz

    xlabel('z')
    #ylabel('field quantity')
    ylabel('vz/c0')
    if (type=='rho'):
        semilogy(z,image)
        semilogy(z,exp(-z**2/2.),linestyle='--')
        semilogy(z,4.e-2*exp(-z/2.),linestyle='--')
        ylim(1.e-4,1.)
    else:
        plot(z,image)
#        semilogy(z,image)

    if save==1:
        savefig('mean_zprofile.eps')

    print -image[0]+image[nz-1]

    return
