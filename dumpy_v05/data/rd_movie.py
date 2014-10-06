import tables
import numpy as np

from dumpy_v05.utils.utils import *

class DumsesMovie:

    def __init__(self,idump=1,nbuf=3):

        #Dump directory name
        dir='movie_'+str_suffix(idump)+'/'

        #Read first output to get basic info
        filename=dir+'rho.000000'
        data=DumsesMovieFile_h5(filename)
        [nxslice,nyslice,nzslice]=data.slice
        [nxglob,nyglob,nzglob]   =data.dim_glob
        ndim=1
        if (nyglob>1):
            ndim=2
            if (nzglob>1):
                ndim=3
        nx=nxglob/nxslice ; ny=nyglob/nyslice ; nz=nzglob/nzslice
        npes=nxslice*nyslice*nzslice
        
        self.time=data.time

        #Define output arrays
        self.x   =np.zeros(nxglob)
        self.y   =np.zeros(nyglob)
        self.z   =np.zeros(nzglob)
        self.rho =np.zeros((nxglob,nyglob,nzglob  ))

        #Loop on cpus
        for mype in range(npes):

            filename=dir+'rho.'+str_suffix(mype)
            
            xposition = mype % nxslice
            yposition = mype/(nxslice*nzslice) % nyslice
            zposition = mype/(nxslice) % nzslice
            
            i0=xposition*nx ; i1=(xposition+1)*nx
            j0=yposition*ny ; j1=(yposition+1)*ny
            k0=zposition*nz ; k1=(zposition+1)*nz

            data=DumsesMovieFile_h5(filename)

            self.x[i0:i1]=data.x[nbuf:nx+nbuf]
            if nyglob>1:
                self.y[j0:j1]=data.y[nbuf:ny+nbuf]
            else:
                self.y[j0:j1]=data.y[0:1]
            if nzglob>1:
                self.z[k0:k1]=data.z[nbuf:nz+nbuf]
            else:
                self.z[k0:k1]=data.z[0:1]

            self.rho [i0:i1,j0:j1,k0:k1  ]=GetSubScalar(data.xyz,nx,ny,nz,nbuf,ndim)

    def get_1d_y(self,type='rho'):
        if type=='rho':
            slice=self.rho[0,:,0]
        elif type=='vx':
            slice=self.rhou[0,:,0,0]/self.rho[0,:,0]
        elif type=='vy':
            slice=self.rhou[0,:,0,1]/self.rho[0,:,0]
        elif type=='vz':
            slice=self.rhou[0,:,0,2]/self.rho[0,:,0]
        elif type=='bx':
            slice=self.B[0,:,0,0]
        elif type=='by':
            slice=self.B[0,:,0,1]
        elif type=='bz':
            slice=self.B[0,:,0,2]#/self.B[0,:,0,1]
        return slice
                        
    def get_2d_xy(self,type='rho'):
        if type=='rho':
            slice=self.rho[:,:,0]
        elif type=='vx':
            slice=self.rhou[:,:,0,0]/self.rho[:,:,0]
        elif type=='vy':
            slice=self.rhou[:,:,0,1]/self.rho[:,:,0]
        elif type=='vz':
            slice=self.rhou[:,:,0,2]/self.rho[:,:,0]
        elif type=='bx':
            slice=self.B[:,:,0,0]
        elif type=='by':
            slice=self.B[:,:,0,1]
        elif type=='bz':
            slice=self.B[:,:,0,2]
        return slice
                        
    def get_2d_xz(self,type='rho'):
        if type=='rho':
            slice=self.rho[:,0,:]
        elif type=='vx':
            slice=self.rhou[:,0,:,0]/self.rho[:,0,:]
        elif type=='vy':
            slice=self.rhou[:,0,:,1]/self.rho[:,0,:]
        elif type=='vz':
            slice=self.rhou[:,0,:,2]/self.rho[:,0,:]
        elif type=='bx':
            slice=self.B[:,0,:,0]
        elif type=='by':
            slice=self.B[:,0,:,1]
        elif type=='bz':
            slice=self.B[:,0,:,2]
        return slice
                        
    def get_yz_mean(self):

        mean_data=np.mean(np.mean(self.rho,1),1)
        return mean_data

class DumsesMovieFile_h5:

    def __init__(self,filename='rho.000000',nbuf=3):
        """ Read HDF5 data output from Dumses."""

        #Open file
        f=tables.openFile(filename)

        #Dataset "para_real"
        self.time=f.root.para_real[0]
        self.dt  =f.root.para_real[1]
        self.dx  =f.root.para_real[2]
        self.dy  =f.root.para_real[3]
        self.dz  =f.root.para_real[4]

        #Dataset "para_int"
        self.nspec  =f.root.para_int[0]
        self.nx     =f.root.para_int[1]
        self.ny     =f.root.para_int[2]
        self.nz     =f.root.para_int[3]
        self.nxslice=f.root.para_int[4]
        self.nyslice=f.root.para_int[5]
        self.nzslice=f.root.para_int[6]

        self.dim    =f.root.para_int[1:4]
        self.slice  =f.root.para_int[4:7]
        self.dim_glob=self.dim*self.slice

        #Dataset "x", "y" and "z
        self.x=f.root.x[:]
        self.y=f.root.y[:]
        self.z=f.root.z[:]

        #Dataset "rho"
        self.xyz=f.root.rho[:,:,:]

        #Dataset "para_mpi"
        if (self.nxslice*self.nyslice*self.nzslice>1):
            self.xleft =f.root.para_mpi[0]
            self.xright=f.root.para_mpi[1]
            self.yleft =f.root.para_mpi[2]
            self.yright=f.root.para_mpi[3]
            self.zleft =f.root.para_mpi[4]
            self.zright=f.root.para_mpi[5]
            self.xposition=f.root.para_mpi[6]
            self.yposition=f.root.para_mpi[7]
            self.zposition=f.root.para_mpi[8]

        #Close file
        f.close()

def get_array(fid,nb_count,type):
    bits_nb='i4'
    pad=np.fromfile(fid,count=1,dtype=bits_nb)
    array=np.fromfile(fid,count=nb_count,dtype=type)
    pad=np.fromfile(fid,count=1,dtype=bits_nb)
    return array

def GetSubArray3d(array,ivar,nx,ny,nz,nbuf,ndim):
    if ndim==1:
        SubArray3d=np.transpose(array[ivar,0:nz,0:ny,nbuf:nx+nbuf],(2,1,0))  
    if ndim==2:
        SubArray3d=np.transpose(array[ivar,0:nz,nbuf:ny+nbuf,nbuf:nx+nbuf],(2,1,0))  
    if ndim==3:
        SubArray3d=np.transpose(array[ivar,nbuf:nz+nbuf,nbuf:ny+nbuf,nbuf:nx+nbuf],(2,1,0))  
    
def GetSubScalar(array,nx,ny,nz,nbuf,ndim):
    if ndim==1:
        SubScalar=np.transpose(array[0:nz,0:ny,nbuf:nx+nbuf],(2,1,0))  
    if ndim==2:
        SubScalar=np.transpose(array[0:nz,nbuf:ny+nbuf,nbuf:nx+nbuf],(2,1,0))  
    if ndim==3:
        SubScalar=np.transpose(array[nbuf:nz+nbuf,nbuf:ny+nbuf,nbuf:nx+nbuf],(2,1,0))  
    
    return SubScalar
