import tables
import numpy as np
import netCDF4 as ncdf

from dumpy_v05.utils.utils import *

class DumsesData:

    def __init__(self,idump=1,filedir='./',nvar=8,nbuf=3,bin=False,collective=False,ghost=False,pnetcdf=False,inline=False,gather=False,field=None):

        if pnetcdf: 
            collective = True
            ghost      = True
        if inline:
            if not(gather):
                ghost  = False
        if gather:
            collective = True
        #Read data from file
        if collective:
            if pnetcdf:
                data=DumsesDataCollective_nc(idump=idump,filedir=filedir,nvar=nvar,nbuf=nbuf,ghost=ghost,inline=inline,field=field)
            else:
                data=DumsesDataCollective_h5(idump=idump,filedir=filedir,nvar=nvar,nbuf=nbuf,ghost=ghost,inline=inline,gather=gather,field=field)
        else:
            data=DumsesDataNoCollective(idump=idump,filedir=filedir,nvar=nvar,nbuf=nbuf,bin=bin,field=field)

        #Define output data
        self.time=data.time

        #Define output arrays
        self.x = data.x
        self.y = data.y
        self.z = data.z
        if not field:
            self.rho  = data.rho
            self.rhou = data.rhou
            self.E    = data.E
            self.B    = data.B
        else:
            if field == 'rho':   self.rho   = data.rho
            if field == 'E':     self.E     = data.E
            if field == 'rhoux': self.rhoux = data.rhoux
            if field == 'rhouy': self.rhouy = data.rhouy
            if field == 'rhouz': self.rhouz = data.rhouz
            if field == 'Bx':    self.Bx    = data.Bx
            if field == 'By':    self.By    = data.By
            if field == 'Bz':    self.Bz    = data.Bz

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

    def get_xy_mean(self):

        mean_data=np.mean(np.mean(self.rho,0),0)
        return mean_data

    def get_rho_rhd(self,gamma=5./3.):
        (nx,ny,nz)=self.rho.shape
        rho=np.zeros((nx,ny,nz))
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    M=self.rhou[i,j,k,0]**2+self.rhou[i,j,k,1]**2+self.rhou[i,j,k,2]**2
                    M=np.sqrt(M)
                    if M==0:
                        rho[i,j,k]=self.rho[i,j,k]
                    else:
                        D=self.rho[i,j,k] ; E=self.E[i,j,k]
                        #Method taken from Mignone et al,2007
                        # Make a guess for the value of R. Solve 3R**2-4ER+M**2
                        Delta = 16.0*E**2-12.0*M**2
               		rplus = (4.0*E+np.sqrt(Delta))/6.
                            
                        #Use qplus to find Q. Substract density
                        rplus=rplus-D
                        E=E-D
                            
			R=GetRoot(D,M,E,gamma,rplus)
                        
		        #Go back to Q instead of Qprime
                        R=R+D

                        #Get the lorentz factor
                        u2 = M**2/(R**2-M**2)
                        lor=(1.0+u2)**(1./2.)
                        rho[i,j,k]=self.rho[i,j,k]/lor
         
        return rho

    def get_v_rhd(self,gamma=5./3.):
        (nx,ny,nz)=self.rho.shape
        v_rhd=np.zeros((nx,ny,nz,3))
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    M=self.rhou[i,j,k,0]**2+self.rhou[i,j,k,1]**2+self.rhou[i,j,k,2]**2
                    M=np.sqrt(M)
                    if M==0:
                        v_rhd[i,j,k,:]=0
                    else:
                        D=self.rho[i,j,k] ; E=self.E[i,j,k]
                        Delta = 16.0*E**2-12.0*M**2
               		rplus = (4.0*E+np.sqrt(Delta))/6.
               		
                        rplus=rplus-D
                        E=E-D
                        
                        #Use qplus to find the root
			R=GetRoot(D,M,E,gamma,rplus)
                        
                        #Go back to Q instead of Qprime
                        R=R+D
                        
                        v_rhd[i,j,k,0] = self.rhou[i,j,k,0]/R
                	v_rhd[i,j,k,1] = self.rhou[i,j,k,1]/R		
			v_rhd[i,j,k,2] = self.rhou[i,j,k,2]/R	
		
        return v_rhd

    def get_p_rhd(self,gamma=5./3.):
        (nx,ny,nz)=self.rho.shape
        p=np.zeros((nx,ny,nz))
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    M=self.rhou[i,j,k,0]**2+self.rhou[i,j,k,1]**2+self.rhou[i,j,k,2]**2
                    M=np.sqrt(M)
                    if M==0:
                        p[i,j,k]=(gamma-1)*(self.E[i,j,k]-self.rho[i,j,k])
                    else:
                        D=self.rho[i,j,k] ; E=self.E[i,j,k]
                        Delta = 16.0*E**2-12.0*M**2
               		rplus = (4.0*E+np.sqrt(Delta))/6.
               		
                        rplus=rplus-D
                        E=E-D
                                
                        #Use qplus to find the root
			R=GetRoot(D,M,E,gamma,rplus)	
			#Go back to Q instead of Qprime
                        R=R+D

                        #Get the lorentz factor
                        u2 = M**2/(R**2-M**2)
                        lor=(1.0+u2)**(1./2.)
                        p[i,j,k]=((R-self.rho[i,j,k])/lor**2-D*u2/((1.+lor)*lor**2))*(gamma-1.)/gamma
        return p


class DumsesDataNoCollective:
    def __init__(self,idump=1,filedir='./',nvar=8,nbuf=3,bin=False, field=None):
        #Dump directory name
        dir=filedir+'output_'+str_suffix(idump)+'/'

        #Read first output to get basic info
        filename=dir+'slices.000000'
        if bin:
            data=DumsesFile(filename)
        else:
            data=DumsesFile_h5(filename)
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
        if not field:
            self.rho =np.zeros((nxglob,nyglob,nzglob  ))
            self.rhou=np.zeros((nxglob,nyglob,nzglob,3))
            self.E   =np.zeros((nxglob,nyglob,nzglob  ))
            self.B   =np.zeros((nxglob,nyglob,nzglob,3))
        else:
            if field == 'rho'  : self.rho   = np.zeros((nxglob,nyglob,nzglob))
            if field == 'E'    : self.E     = np.zeros((nxglob,nyglob,nzglob))
            if field == 'rhoux': self.rhoux = np.zeros((nxglob,nyglob,nzglob))
            if field == 'rhouy': self.rhouy = np.zeros((nxglob,nyglob,nzglob))
            if field == 'rhouz': self.rhouz = np.zeros((nxglob,nyglob,nzglob))
            if field == 'Bx'   : self.Bx    = np.zeros((nxglob,nyglob,nzglob))
            if field == 'By'   : self.By    = np.zeros((nxglob,nyglob,nzglob))
            if field == 'Bz'   : self.Bz    = np.zeros((nxglob,nyglob,nzglob))

        #Loop on cpus
        for mype in range(npes):

            filename=dir+'slices.'+str_suffix(mype)
            
            xposition = mype % nxslice
            yposition = mype/(nxslice*nzslice) % nyslice
            zposition = mype/(nxslice) % nzslice
            
            i0=xposition*nx ; i1=(xposition+1)*nx
            j0=yposition*ny ; j1=(yposition+1)*ny
            k0=zposition*nz ; k1=(zposition+1)*nz

            if bin:
                data=DumsesFile(filename,field=field)
            else:
                data=DumsesFile_h5(filename,field=field)

            self.x[i0:i1]=data.x[nbuf:nx+nbuf]
            if nyglob>1:
                self.y[j0:j1]=data.y[nbuf:ny+nbuf]
            else:
                self.y[j0:j1]=data.y[0:1]
            if nzglob>1:
                self.z[k0:k1]=data.z[nbuf:nz+nbuf]
            else:
                self.z[k0:k1]=data.z[0:1]

            if not field:
                self.rho [i0:i1,j0:j1,k0:k1  ]=GetSubArray3d(data.xyz,0,nx,ny,nz,nbuf,ndim)
                self.rhou[i0:i1,j0:j1,k0:k1,0]=GetSubArray3d(data.xyz,1,nx,ny,nz,nbuf,ndim)
                self.rhou[i0:i1,j0:j1,k0:k1,1]=GetSubArray3d(data.xyz,2,nx,ny,nz,nbuf,ndim)
                self.rhou[i0:i1,j0:j1,k0:k1,2]=GetSubArray3d(data.xyz,3,nx,ny,nz,nbuf,ndim)
                self.E   [i0:i1,j0:j1,k0:k1  ]=GetSubArray3d(data.xyz,4,nx,ny,nz,nbuf,ndim)
                self.B   [i0:i1,j0:j1,k0:k1,0]=GetSubArray3d(data.xyz,5,nx,ny,nz,nbuf,ndim)
                self.B   [i0:i1,j0:j1,k0:k1,1]=GetSubArray3d(data.xyz,6,nx,ny,nz,nbuf,ndim)
                self.B   [i0:i1,j0:j1,k0:k1,2]=GetSubArray3d(data.xyz,7,nx,ny,nz,nbuf,ndim)
            else:
                if field == 'rho'  : 
                    self.rho[i0:i1,j0:j1,k0:k1]   = GetArray3d(data.xyz,nx,ny,nz,nbuf,ndim)
                if field == 'E'    :
                    self.E[i0:i1,j0:j1,k0:k1]     = GetArray3d(data.xyz,nx,ny,nz,nbuf,ndim)
                if field == 'rhoux':
                    self.rhoux[i0:i1,j0:j1,k0:k1] = GetArray3d(data.xyz,nx,ny,nz,nbuf,ndim)
                if field == 'rhouy':
                    self.rhouy[i0:i1,j0:j1,k0:k1] = GetArray3d(data.xyz,nx,ny,nz,nbuf,ndim)
                if field == 'rhouz':
                    self.rhouz[i0:i1,j0:j1,k0:k1] = GetArray3d(data.xyz,nx,ny,nz,nbuf,ndim)
                if field == 'Bx'   :
                    self.Bx[i0:i1,j0:j1,k0:k1]    = GetArray3d(data.xyz,nx,ny,nz,nbuf,ndim)
                if field == 'By'   :
                    self.By[i0:i1,j0:j1,k0:k1]    = GetArray3d(data.xyz,nx,ny,nz,nbuf,ndim)
                if field == 'Bz'   :
                    self.Bz[i0:i1,j0:j1,k0:k1]    = GetArray3d(data.xyz,nx,ny,nz,nbuf,ndim)


class DumsesDataCollective_h5:
    def __init__(self,idump=1,filedir='./',nvar=8,nbuf=3,ghost=False,inline=False,gather=False,field=None):
        #Dump directory name
        dir=filedir+'output_'+str_suffix(idump)+'/'

        #Read first output to get basic info
        if gather:
            filename=dir+'gather.000000'
        else:
            filename=dir+'slices.000000'

        #Open file
        f=tables.openFile(filename)

        #Dataset "para_real"
        self.time=f.root.para_real[0]
        self.dt  =f.root.para_real[1]
        self.dx  =f.root.para_real[2]
        self.dy  =f.root.para_real[3]
        self.dz  =f.root.para_real[4]

        #Dataset "para_int"
        self.ndump  =f.root.para_int[0]
        self.nhist  =f.root.para_int[1]
        self.nspec  =f.root.para_int[2]
        self.nx     =f.root.para_int[3]
        self.ny     =f.root.para_int[4]
        self.nz     =f.root.para_int[5]
        self.nxslice=f.root.para_int[6]
        self.nyslice=f.root.para_int[7]
        self.nzslice=f.root.para_int[8]
        npes = self.nxslice*self.nyslice*self.nzslice
        nx = self.nx; nxslice = self.nxslice
        ny = self.ny; nyslice = self.nyslice
        nz = self.nz; nzslice = self.nzslice

        self.dim    =f.root.para_int[3:6]
        self.slice  =f.root.para_int[6:9]
        self.dim_glob=self.dim*self.slice

        #Dataset "boxSize"
        self.xmin=f.root.boxSize[0]
        self.xmax=f.root.boxSize[1]
        self.ymin=f.root.boxSize[2]
        self.ymax=f.root.boxSize[3]
        self.zmin=f.root.boxSize[4]
        self.zmax=f.root.boxSize[5]

        #Dataset "x", "y" and "z"
        xroot  = f.root.x[:]
        x      = xroot[:,nbuf:-nbuf].flatten()
        x.sort()
        self.x = x[::nyslice*nzslice]

        yroot  = f.root.y[:]
        y      = yroot[:,nbuf:-nbuf].flatten()
        y.sort()
        self.y = y[::nxslice*nzslice]

        zroot  = f.root.z[:]
        z      = zroot[:,nbuf:-nbuf].flatten()
        z.sort()
        self.z = z[::nxslice*nyslice]

        #Get global sizes
        nxglob=self.dim_glob[0]
        nyglob=self.dim_glob[1]
        nzglob=self.dim_glob[2]
        ndim=1
        if (nyglob>1):
            ndim=ndim+1
            if (nzglob>1):
                ndim=ndim+1

        #Define arrays
        self.rhou=np.zeros((nxglob,nyglob,nzglob,ndim))
        self.B=np.zeros((nxglob,nyglob,nzglob,ndim))
        if inline:
            if not field:
                rho   = np.zeros((self.dim_glob))
                rhoux = np.zeros((self.dim_glob))
                rhouy = np.zeros((self.dim_glob))
                rhouz = np.zeros((self.dim_glob))
                E     = np.zeros((self.dim_glob))
                Bx    = np.zeros((self.dim_glob))
                By    = np.zeros((self.dim_glob))
                Bz    = np.zeros((self.dim_glob))
            else:
                if field == 'rho'  : rho   = np.zeros((self.dim_glob))
                if field == 'E'    : E     = np.zeros((self.dim_glob))
                if field == 'rhoux': rhoux = np.zeros((self.dim_glob))
                if field == 'rhouy': rhouy = np.zeros((self.dim_glob))
                if field == 'rhouz': rhouz = np.zeros((self.dim_glob))
                if field == 'Bx'   : Bx    = np.zeros((self.dim_glob))
                if field == 'By'   : By    = np.zeros((self.dim_glob))
                if field == 'Bz'   : Bz    = np.zeros((self.dim_glob))

        if ghost:
            if (gather and inline):
                selx = [i+nbuf for i in range(self.nx)]
                sely = [i+nbuf for i in range(self.ny)]
                selz = [i+nbuf for i in range(self.nz)]
            else:
                selx = np.array([])
                sely = np.array([])
                selz = np.array([])
                for i in range(self.nxslice):
                    temp = np.linspace(i*self.nx+(2*i+1)*nbuf,(i+1)*self.nx+(2*i+1)*nbuf-1,self.nx)
                    selx = np.append(selx,temp)
                for i in range(self.nyslice):
                    temp = np.linspace(i*self.ny+(2*i+1)*nbuf,(i+1)*self.ny+(2*i+1)*nbuf-1,self.ny)
                    sely = np.append(sely,temp)
                for i in range(self.nzslice):
                    temp = np.linspace(i*self.nz+(2*i+1)*nbuf,(i+1)*self.nz+(2*i+1)*nbuf-1,self.nz)
                    selz = np.append(selz,temp)
                selx = selx.astype('int')
                sely = sely.astype('int')
                selz = selz.astype('int')
        elif inline:
            if (gather and not(ghost)):
                selx = [i for i in range(self.nx)]
                sely = [i for i in range(self.ny)]
                selz = [i for i in range(self.nz)]
            else:
                selx = [i+nbuf for i in range(self.nx)]
                sely = [i+nbuf for i in range(self.ny)]
                selz = [i+nbuf for i in range(self.nz)]
        elif gather:
            selx = [i for i in range(nxglob)]
            sely = [i for i in range(nyglob)]
            selz = [i for i in range(nzglob)]
        else:
            selx = [i+nbuf for i in range(nxglob)]
            sely = [i+nbuf for i in range(nyglob)]
            selz = [i+nbuf for i in range(nzglob)]

        #Get datacubes (rho,E,rhou,B)
        if inline:
            if (gather and not(ghost)):
                nex = nx; ney = ny; nez = nz
            else:
                nex = nx+2*nbuf; ney = ny+2*nbuf; nez = nz+2*nbuf
            if not field:
                rhor   = f.root.rho[:]
                rhor   = np.reshape(rhor, (npes, nez, ney, nex))
                rhor   = np.transpose(rhor, (0, 3, 2, 1))
                rhouxr = f.root.rho_vx[:]
                rhouxr = np.reshape(rhouxr, (npes, nez, ney, nex))
                rhouxr = np.transpose(rhouxr, (0, 3, 2, 1))
                rhouyr = f.root.rho_vy[:]
                rhouyr = np.reshape(rhouyr, (npes, nez, ney, nex))
                rhouyr = np.transpose(rhouyr, (0, 3, 2, 1))
                rhouzr = f.root.rho_vz[:]
                rhouzr = np.reshape(rhouzr, (npes, nez, ney, nex))
                rhouzr = np.transpose(rhouzr, (0, 3, 2, 1))
                Er     = f.root.E[:]
                Er     = np.reshape(Er, (npes, nez, ney, nex))
                Er     = np.transpose(Er, (0, 3, 2, 1))
                Bxr    = f.root.Bx[:]
                Bxr    = np.reshape(Bxr, (npes, nez, ney, nex))
                Bxr    = np.transpose(Bxr, (0, 3, 2, 1))
                Byr    = f.root.By[:]
                Byr    = np.reshape(Byr, (npes, nez, ney, nex))
                Byr    = np.transpose(Byr, (0, 3, 2, 1))
                Bzr    = f.root.Bz[:]
                Bzr    = np.reshape(Bzr, (npes, nez, ney, nex))
                Bzr    = np.transpose(Bzr, (0, 3, 2, 1))
            else:
                if field == 'rho':
                    rhor   = f.root.rho[:]
                    rhor   = np.reshape(rhor, (npes, nez, ney, nex))
                    rhor   = np.transpose(rhor, (0, 3, 2, 1))
                if field == 'rhoux':
                    rhouxr = f.root.rho_vx[:]
                    rhouxr = np.reshape(rhouxr, (npes, nez, ney, nex))
                    rhouxr = np.transpose(rhouxr, (0, 3, 2, 1))
                if field == 'rhouy':
                    rhouyr = f.root.rho_vy[:]
                    rhouyr = np.reshape(rhouyr, (npes, nez, ney, nex))
                    rhouyr = np.transpose(rhouyr, (0, 3, 2, 1))
                if field == 'rhouz':
                    rhouzr = f.root.rho_vz[:]
                    rhouzr = np.reshape(rhouzr, (npes, nez, ney, nex))
                    rhouzr = np.transpose(rhouzr, (0, 3, 2, 1))
                if field == 'E':
                    Er     = f.root.E[:]
                    Er     = np.reshape(Er, (npes, nez, ney, nex))
                    Er     = np.transpose(Er, (0, 3, 2, 1))
                if field == 'Bx':
                    Bxr    = f.root.Bx[:]
                    Bxr    = np.reshape(Bxr, (npes, nez, ney, nex))
                    Bxr    = np.transpose(Bxr, (0, 3, 2, 1))
                if field == 'By':
                    Byr    = f.root.By[:]
                    Byr    = np.reshape(Byr, (npes, nez, ney, nex))
                    Byr    = np.transpose(Byr, (0, 3, 2, 1))
                if field == 'Bz':
                    Bzr    = f.root.Bz[:]
                    Bzr    = np.reshape(Bzr, (npes, nez, ney, nex))
                    Bzr    = np.transpose(Bzr, (0, 3, 2, 1))

            for n in range(npes):
                xposition = n%nxslice
                yposition = n/(nxslice*nzslice)%nyslice
                zposition = n/(nxslice)%nzslice
                i0=xposition*nx ; i1=(xposition+1)*nx
                j0=yposition*ny ; j1=(yposition+1)*ny
                k0=zposition*nz ; k1=(zposition+1)*nz
                if not field:
                    rho[i0:i1,j0:j1,k0:k1]   = rhor[n][selx,:,:][:,sely,:][:,:,selz]
                    rhoux[i0:i1,j0:j1,k0:k1] = rhouxr[n][selx,:,:][:,sely,:][:,:,selz]
                    rhouy[i0:i1,j0:j1,k0:k1] = rhouyr[n][selx,:,:][:,sely,:][:,:,selz]
                    rhouz[i0:i1,j0:j1,k0:k1] = rhouzr[n][selx,:,:][:,sely,:][:,:,selz]
                    E[i0:i1,j0:j1,k0:k1]     = Er[n][selx,:,:][:,sely,:][:,:,selz]
                    Bx[i0:i1,j0:j1,k0:k1]    = Bxr[n][selx,:,:][:,sely,:][:,:,selz]
                    By[i0:i1,j0:j1,k0:k1]    = Byr[n][selx,:,:][:,sely,:][:,:,selz]
                    Bz[i0:i1,j0:j1,k0:k1]    = Bzr[n][selx,:,:][:,sely,:][:,:,selz]
                else:
                    if field == 'rho'  : 
                        rho[i0:i1,j0:j1,k0:k1]   = rhor[n][selx,:,:][:,sely,:][:,:,selz]
                    if field == 'rhoux': 
                        rhoux[i0:i1,j0:j1,k0:k1] = rhouxr[n][selx,:,:][:,sely,:][:,:,selz]
                    if field == 'rhouy': 
                        rhouy[i0:i1,j0:j1,k0:k1] = rhouyr[n][selx,:,:][:,sely,:][:,:,selz]
                    if field == 'rhouz': 
                        rhouz[i0:i1,j0:j1,k0:k1] = rhouzr[n][selx,:,:][:,sely,:][:,:,selz]
                    if field == 'E'    : 
                        E[i0:i1,j0:j1,k0:k1]     = Er[n][selx,:,:][:,sely,:][:,:,selz]
                    if field == 'Bx'   : 
                        Bx[i0:i1,j0:j1,k0:k1]    = Bxr[n][selx,:,:][:,sely,:][:,:,selz]
                    if field == 'By'   : 
                        By[i0:i1,j0:j1,k0:k1]    = Byr[n][selx,:,:][:,sely,:][:,:,selz]
                    if field == 'Bz'   : 
                        Bz[i0:i1,j0:j1,k0:k1]    = Bzr[n][selx,:,:][:,sely,:][:,:,selz]



            if not field:        
                self.rho  = rho
                self.rhou = np.array([rhoux, rhouy, rhouz])
                self.E    = E
                self.B    = np.array([Bx, By, Bz])
            else:
                if field == 'rho'  : self.rho   = rho
                if field == 'E'    : self.E     = E
                if field == 'rhoux': self.rhoux = rhoux
                if field == 'rhouy': self.rhouy = rhouy
                if field == 'rhouz': self.rhouz = rhouz
                if field == 'Bx'   : self.Bx    = Bx
                if field == 'By'   : self.By    = By
                if field == 'Bz'   : self.Bz    = Bz
        else:
            if not field:
                rho = np.transpose(f.root.rho[:,:,:],(2,1,0))
                self.rho = rho[selx,:,:][:,sely,:][:,:,selz]
                E = np.transpose(f.root.E[:,:,:],(2,1,0))
                self.E = E[selx,:,:][:,sely,:][:,:,selz]
                rhoux = np.transpose(f.root.rho_vx[:,:,:],(2,1,0))
                self.rhou[:,:,:,0] = rhoux[selx,:,:][:,sely,:][:,:,selz]
                rhouy=np.transpose(f.root.rho_vy[:,:,:],(2,1,0))
                self.rhou[:,:,:,1] = rhouy[selx,:,:][:,sely,:][:,:,selz]
                rhouz=np.transpose(f.root.rho_vz[:,:,:],(2,1,0))
                self.rhou[:,:,:,2] = rhouz[selx,:,:][:,sely,:][:,:,selz]
                Bx=np.transpose(f.root.Bx[:,:,:],(2,1,0))
                self.B[:,:,:,0] = Bx[selx,:,:][:,sely,:][:,:,selz]
                By=np.transpose(f.root.By[:,:,:],(2,1,0))
                self.B[:,:,:,1] = By[selx,:,:][:,sely,:][:,:,selz]
                Bz=np.transpose(f.root.Bz[:,:,:],(2,1,0))
                self.B[:,:,:,2] = Bz[selx,:,:][:,sely,:][:,:,selz]
            else:
                if field == 'rho'  : 
                    rho = np.transpose(f.root.rho[:,:,:],(2,1,0))
                    self.rho = rho[selx,:,:][:,sely,:][:,:,selz]
                if field == 'E'    : 
                    E = np.transpose(f.root.E[:,:,:],(2,1,0))
                    self.E = E[selx,:,:][:,sely,:][:,:,selz]
                if field == 'rhoux': 
                    rhoux = np.transpose(f.root.rho_vx[:,:,:],(2,1,0))
                    self.rhoux[:,:,:] = rhoux[selx,:,:][:,sely,:][:,:,selz]
                if field == 'rhouy': 
                    rhouy = np.transpose(f.root.rho_vy[:,:,:],(2,1,0))
                    self.rhouy[:,:,:] = rhouy[selx,:,:][:,sely,:][:,:,selz]
                if field == 'rhouz': 
                    rhouz = np.transpose(f.root.rho_vz[:,:,:],(2,1,0))
                    self.rhouz[:,:,:] = rhouz[selx,:,:][:,sely,:][:,:,selz]
                if field == 'Bx'   : 
                    Bx = np.transpose(f.root.Bx[:,:,:],(2,1,0))
                    self.Bx[:,:,:] = Bx[selx,:,:][:,sely,:][:,:,selz]
                if field == 'By'   : 
                    By = np.transpose(f.root.By[:,:,:],(2,1,0))
                    self.By[:,:,:] = By[selx,:,:][:,sely,:][:,:,selz]
                if field == 'Bz'   : 
                    Bz = np.transpose(f.root.Bz[:,:,:],(2,1,0))
                    self.Bz[:,:,:] = Bz[selx,:,:][:,sely,:][:,:,selz]
        #Close file
        f.close()

class DumsesDataCollective_nc:
    def __init__(self,idump=1,filedir='./',nvar=8,nbuf=3,ghost=True,inline=False,field=None):
        #Dump directory name
        dir=filedir+'output_'+str_suffix(idump)+'/'

        #Read first output to get basic info
        filename=dir+'slices.000000'

        #Open file
        f=ncdf.Dataset(filename)

        #Dataset "para_real"
        self.time=f.variables['para_real'][0]
        self.dt  =f.variables['para_real'][1]
        self.dx  =f.variables['para_real'][2]
        self.dy  =f.variables['para_real'][3]
        self.dz  =f.variables['para_real'][4]

        #Dataset "para_int"
        self.ndump  =f.variables['para_int'][0]
        self.nhist  =f.variables['para_int'][1]
        self.nspec  =f.variables['para_int'][2]
        self.nx     =f.variables['para_int'][3]
        self.ny     =f.variables['para_int'][4]
        self.nz     =f.variables['para_int'][5]
        self.nxslice=f.variables['para_int'][6]
        self.nyslice=f.variables['para_int'][7]
        self.nzslice=f.variables['para_int'][8]
        npes = self.nxslice*self.nyslice*self.nzslice
        nx = self.nx; nxslice = self.nxslice
        ny = self.ny; nyslice = self.nyslice
        nz = self.nz; nzslice = self.nzslice

        self.dim    =f.variables['para_int'][3:6]
        self.slice  =f.variables['para_int'][6:9]
        self.dim_glob=self.dim*self.slice

        #Dataset "boxSize"
        self.xmin=f.variables['boxSize'][0]
        self.xmax=f.variables['boxSize'][1]
        self.ymin=f.variables['boxSize'][2]
        self.ymax=f.variables['boxSize'][3]
        self.zmin=f.variables['boxSize'][4]
        self.zmax=f.variables['boxSize'][5]

        #Dataset "x", "y" and "z"
        self.x=np.linspace(self.xmin,self.xmax,self.dim_glob[0])
        self.y=np.linspace(self.xmin,self.xmax,self.dim_glob[1])
        self.z=np.linspace(self.xmin,self.xmax,self.dim_glob[2])

        #Get global sizes
        nxglob=self.dim_glob[0]
        nyglob=self.dim_glob[1]
        nzglob=self.dim_glob[2]
        ndim=1
        if (nyglob>1):
            ndim=ndim+1
            if (nzglob>1):
                ndim=ndim+1

        #Define arrays
        self.rhou=np.zeros((nxglob,nyglob,nzglob,ndim))
        self.B=np.zeros((nxglob,nyglob,nzglob,ndim))
        if inline:
            if not field:
                rho   = np.zeros((self.dim_glob))
                rhoux = np.zeros((self.dim_glob))
                rhouy = np.zeros((self.dim_glob))
                rhouz = np.zeros((self.dim_glob))
                E     = np.zeros((self.dim_glob))
                Bx    = np.zeros((self.dim_glob))
                By    = np.zeros((self.dim_glob))
                Bz    = np.zeros((self.dim_glob))
            else:
                if field == 'rho'  : rho   = np.zeros((self.dim_glob))
                if field == 'E'    : E     = np.zeros((self.dim_glob))
                if field == 'rhoux': rhoux = np.zeros((self.dim_glob))
                if field == 'rhouy': rhouy = np.zeros((self.dim_glob))
                if field == 'rhouz': rhouz = np.zeros((self.dim_glob))
                if field == 'Bx'   : Bx    = np.zeros((self.dim_glob))
                if field == 'By'   : By    = np.zeros((self.dim_glob))
                if field == 'Bz'   : Bz    = np.zeros((self.dim_glob))

        if ghost:
            selx = np.array([])
            sely = np.array([])
            selz = np.array([])
            for i in range(self.nxslice):
                temp = np.linspace(i*self.nx+(2*i+1)*nbuf,(i+1)*self.nx+(2*i+1)*nbuf-1,self.nx)
                selx = np.append(selx,temp)
            for i in range(self.nyslice):
                temp = np.linspace(i*self.ny+(2*i+1)*nbuf,(i+1)*self.ny+(2*i+1)*nbuf-1,self.ny)
                sely = np.append(sely,temp)
            for i in range(self.nzslice):
                temp = np.linspace(i*self.nz+(2*i+1)*nbuf,(i+1)*self.nz+(2*i+1)*nbuf-1,self.nz)
                selz = np.append(selz,temp)
            selx = selx.astype('int')
            sely = sely.astype('int')
            selz = selz.astype('int')
        elif inline:
            selx = [i+nbuf for i in range(self.nx)]
            sely = [i+nbuf for i in range(self.ny)]
            selz = [i+nbuf for i in range(self.nz)]
        else:
            selx = [i+nbuf for i in range(nxglob)]
            sely = [i+nbuf for i in range(nyglob)]
            selz = [i+nbuf for i in range(nzglob)]

        #Get datacubes (rho,E,rhou,B)
        if inline:
            if not field:
                rhor   = f.variables['rho'][:]
                rhor   = np.reshape(rhor, (npes, nz+2*nbuf, ny+2*nbuf, nx+2*nbuf))
                rhor   = np.transpose(rhor, (0, 3, 2, 1))
                rhouxr = f.variables['rho_vx'][:]
                rhouxr = np.reshape(rhouxr, (npes, nz+2*nbuf, ny+2*nbuf, nx+2*nbuf))
                rhouxr = np.transpose(rhouxr, (0, 3, 2, 1))
                rhouyr = f.variables['rho_vy'][:]
                rhouyr = np.reshape(rhouyr, (npes, nz+2*nbuf, ny+2*nbuf, nx+2*nbuf))
                rhouyr = np.transpose(rhouyr, (0, 3, 2, 1))
                rhouzr = f.variables['rho_vz'][:]
                rhouzr = np.reshape(rhouzr, (npes, nz+2*nbuf, ny+2*nbuf, nx+2*nbuf))
                rhouzr = np.transpose(rhouzr, (0, 3, 2, 1))
                Er     = f.variables['E'][:]
                Er     = np.reshape(Er, (npes, nz+2*nbuf, ny+2*nbuf, nx+2*nbuf))
                Er     = np.transpose(Er, (0, 3, 2, 1))
                Bxr    = f.variables['Bx'][:]
                Bxr    = np.reshape(Bxr, (npes, nz+2*nbuf, ny+2*nbuf, nx+2*nbuf))
                Bxr    = np.transpose(Bxr, (0, 3, 2, 1))
                Byr    = f.variables['By'][:]
                Byr    = np.reshape(Byr, (npes, nz+2*nbuf, ny+2*nbuf, nx+2*nbuf))
                Byr    = np.transpose(Byr, (0, 3, 2, 1))
                Bzr    = f.variables['Bz'][:]
                Bzr    = np.reshape(Bzr, (npes, nz+2*nbuf, ny+2*nbuf, nx+2*nbuf))
                Bzr    = np.transpose(Bzr, (0, 3, 2, 1))
            else:
                if field == 'rho':
                    rhor   = f.variables['rho'][:]
                    rhor   = np.reshape(rhor, (npes, nz+2*nbuf, ny+2*nbuf, nx+2*nbuf))
                    rhor   = np.transpose(rhor, (0, 3, 2, 1))
                if field == 'rhoux':
                    rhouxr = f.variables['rho_vx'][:]
                    rhouxr = np.reshape(rhouxr, (npes, nz+2*nbuf, ny+2*nbuf, nx+2*nbuf))
                    rhouxr = np.transpose(rhouxr, (0, 3, 2, 1))
                if field == 'rhouy':
                    rhouyr = f.variables['rho_vy'][:]
                    rhouyr = np.reshape(rhouyr, (npes, nz+2*nbuf, ny+2*nbuf, nx+2*nbuf))
                    rhouyr = np.transpose(rhouyr, (0, 3, 2, 1))
                if field == 'rhouz':
                    rhouzr = f.variables['rho_vz'][:]
                    rhouzr = np.reshape(rhouzr, (npes, nz+2*nbuf, ny+2*nbuf, nx+2*nbuf))
                    rhouzr = np.transpose(rhouzr, (0, 3, 2, 1))
                if field == 'E':
                    Er     = f.variables['E'][:]
                    Er     = np.reshape(Er, (npes, nz+2*nbuf, ny+2*nbuf, nx+2*nbuf))
                    Er     = np.transpose(Er, (0, 3, 2, 1))
                if field == 'Bx':
                    Bxr    = f.variables['Bx'][:]
                    Bxr    = np.reshape(Bxr, (npes, nz+2*nbuf, ny+2*nbuf, nx+2*nbuf))
                    Bxr    = np.transpose(Bxr, (0, 3, 2, 1))
                if field == 'By':
                    Byr    = f.variables['By'][:]
                    Byr    = np.reshape(Byr, (npes, nz+2*nbuf, ny+2*nbuf, nx+2*nbuf))
                    Byr    = np.transpose(Byr, (0, 3, 2, 1))
                if field == 'Bz':
                    Bzr    = f.variables['Bz'][:]
                    Bzr    = np.reshape(Bzr, (npes, nz+2*nbuf, ny+2*nbuf, nx+2*nbuf))
                    Bzr    = np.transpose(Bzr, (0, 3, 2, 1))

            for n in range(npes):
                xposition = n%nxslice
                yposition = n/(nxslice*nzslice)%nyslice
                zposition = n/(nxslice)%nzslice
                i0=xposition*nx ; i1=(xposition+1)*nx
                j0=yposition*ny ; j1=(yposition+1)*ny
                k0=zposition*nz ; k1=(zposition+1)*nz
                if not field:
                    rho[i0:i1,j0:j1,k0:k1]   = rhor[n][selx,:,:][:,sely,:][:,:,selz]
                    rhoux[i0:i1,j0:j1,k0:k1] = rhouxr[n][selx,:,:][:,sely,:][:,:,selz]
                    rhouy[i0:i1,j0:j1,k0:k1] = rhouyr[n][selx,:,:][:,sely,:][:,:,selz]
                    rhouz[i0:i1,j0:j1,k0:k1] = rhouzr[n][selx,:,:][:,sely,:][:,:,selz]
                    E[i0:i1,j0:j1,k0:k1]     = Er[n][selx,:,:][:,sely,:][:,:,selz]
                    Bx[i0:i1,j0:j1,k0:k1]    = Bxr[n][selx,:,:][:,sely,:][:,:,selz]
                    By[i0:i1,j0:j1,k0:k1]    = Byr[n][selx,:,:][:,sely,:][:,:,selz]
                    Bz[i0:i1,j0:j1,k0:k1]    = Bzr[n][selx,:,:][:,sely,:][:,:,selz]
                else:
                    if field == 'rho':
                        rho[i0:i1,j0:j1,k0:k1]   = rhor[n][selx,:,:][:,sely,:][:,:,selz]
                    if field == 'rhoux':
                        rhoux[i0:i1,j0:j1,k0:k1] = rhouxr[n][selx,:,:][:,sely,:][:,:,selz]
                    if field == 'rhouy':
                        rhouy[i0:i1,j0:j1,k0:k1] = rhouyr[n][selx,:,:][:,sely,:][:,:,selz]
                    if field == 'rhouz':
                        rhouz[i0:i1,j0:j1,k0:k1] = rhouzr[n][selx,:,:][:,sely,:][:,:,selz]
                    if field == 'E':
                        E[i0:i1,j0:j1,k0:k1]     = Er[n][selx,:,:][:,sely,:][:,:,selz]
                    if field == 'Bx':
                        Bx[i0:i1,j0:j1,k0:k1]    = Bxr[n][selx,:,:][:,sely,:][:,:,selz]
                    if field == 'By':
                        By[i0:i1,j0:j1,k0:k1]    = Byr[n][selx,:,:][:,sely,:][:,:,selz]
                    if field == 'Bz':
                        Bz[i0:i1,j0:j1,k0:k1]    = Bzr[n][selx,:,:][:,sely,:][:,:,selz]

            if not field:
                self.rho  = rho
                self.rhou = np.array([rhoux, rhouy, rhouz])
                self.E    = E
                self.B    = np.array([Bx, By, Bz])
            else:
                if field == 'rho'  : self.rho   = rho
                if field == 'E'    : self.E     = E
                if field == 'rhoux': self.rhoux = rhoux
                if field == 'rhouy': self.rhouy = rhouy
                if field == 'rhouz': self.rhouz = rhouz
                if field == 'Bx'   : self.Bx    = Bx
                if field == 'By'   : self.By    = By
                if field == 'Bz'   : self.Bz    = Bz

        else:
            if not field:
                rho = np.transpose(f.variables['rho'][:,:,:],(2,1,0))
                self.rho = rho[selx,:,:][:,sely,:][:,:,selz]
                E = np.transpose(f.variables['E'][:,:,:],(2,1,0))
                self.E = E[selx,:,:][:,sely,:][:,:,selz]
                rhoux = np.transpose(f.variables['rho_vx'][:,:,:],(2,1,0))
                self.rhou[:,:,:,0] = rhoux[selx,:,:][:,sely,:][:,:,selz]
                rhouy=np.transpose(f.variables['rho_vy'][:,:,:],(2,1,0))
                self.rhou[:,:,:,1] = rhouy[selx,:,:][:,sely,:][:,:,selz]
                rhouz=np.transpose(f.variables['rho_vz'][:,:,:],(2,1,0))
                self.rhou[:,:,:,2] = rhouz[selx,:,:][:,sely,:][:,:,selz]
                Bx=np.transpose(f.variables['Bx'][:,:,:],(2,1,0))
                self.B[:,:,:,0] = Bx[selx,:,:][:,sely,:][:,:,selz]
                By=np.transpose(f.variables['By'][:,:,:],(2,1,0))
                self.B[:,:,:,1] = By[selx,:,:][:,sely,:][:,:,selz]
                Bz=np.transpose(f.variables['Bz'][:,:,:],(2,1,0))
                self.B[:,:,:,2] = Bz[selx,:,:][:,sely,:][:,:,selz]
            else:
                if field == 'rho'  : 
                    rho = np.transpose(f.variables['rho'][:,:,:],(2,1,0))
                    self.rho = rho[selx,:,:][:,sely,:][:,:,selz]
                if field == 'E'    : 
                    E = np.transpose(f.variables['E'][:,:,:],(2,1,0))
                    self.E = E[selx,:,:][:,sely,:][:,:,selz]
                if field == 'rhoux': 
                    rhoux = np.transpose(f.variables['rho_vx'][:,:,:],(2,1,0))
                    self.rhoux[:,:,:] = rhoux[selx,:,:][:,sely,:][:,:,selz]
                if field == 'rhouy': 
                    rhouy=np.transpose(f.variables['rho_vy'][:,:,:],(2,1,0))
                    self.rhouy[:,:,:] = rhouy[selx,:,:][:,sely,:][:,:,selz]
                if field == 'rhouz': 
                    rhouz=np.transpose(f.variables['rho_vz'][:,:,:],(2,1,0))
                    self.rhouz[:,:,:] = rhouz[selx,:,:][:,sely,:][:,:,selz]
                if field == 'Bx'   : 
                    Bx=np.transpose(f.variables['Bx'][:,:,:],(2,1,0))
                    self.Bx[:,:,:] = Bx[selx,:,:][:,sely,:][:,:,selz]
                if field == 'By'   : 
                    By=np.transpose(f.variables['By'][:,:,:],(2,1,0))
                    self.By[:,:,:] = By[selx,:,:][:,sely,:][:,:,selz]
                if field == 'Bz'   : 
                    Bz=np.transpose(f.variables['Bz'][:,:,:],(2,1,0))
                    self.Bz[:,:,:] = Bz[selx,:,:][:,sely,:][:,:,selz]

        #Close file
        f.close()


class DumsesFile:

    def __init__(self,filename='slices.000000',nvar=8,nbuf=3,field=None):

        #Open file
        f=open(filename,'rb')

        #Get data dimensions from file
        a    =get_array(f,5,'f8') ; [self.time,dt,dx,dy,dz]=a
        a    =get_array(f,3,'i4')
        dim  =get_array(f,3,'i4')  ; [nx,ny,nz]=dim
        self.slice=get_array(f,3,'i4')  ; self.dim_glob=dim*self.slice
        
        #Get ndim and allocate arrays
        ndim=1
        self.x=np.zeros(nx+2*nbuf) ; self.y=np.zeros(ny) ; self.z=np.zeros(nz)
        shape=(nvar,1,1,nx+2*nbuf)
        if ny>1:
            ndim=2
            self.y=np.zeros(ny+2*nbuf)
            shape=(nvar,1,ny+2*nbuf,nx+2*nbuf)
            if nz>1:
                ndim=3
                self.z=np.zeros(nz+2*nbuf)
                shape=(nvar,nz+2*nbuf,ny+2*nbuf,nx+2*nbuf)
                
        #Get grid position and state variables
        if ndim==1:
            ntot=(nx+2*nbuf)*ny*nz*nvar
            positions=get_array(f,nx+ny+nz+2*nbuf,'f8')
            self.x  =positions[0           :nx      +2*nbuf]
            self.y  =positions[nx   +2*nbuf:nx+ny   +2*nbuf]
            self.z  =positions[nx+ny+2*nbuf:nx+ny+nz+2*nbuf]
            self.xyz=get_array(f,ntot              ,'f8').reshape(shape)
        if ndim==2:
            ntot=(nx+2*nbuf)*(ny+2*nbuf)*nz*nvar
            positions=get_array(f,nx+ny+nz+4*nbuf,'f8')
            self.x  =positions[0           :nx      +2*nbuf]
            self.y  =positions[nx   +2*nbuf:nx+ny   +4*nbuf]
            self.z  =positions[nx+ny+4*nbuf:nx+ny+nz+4*nbuf]
            self.xyz=get_array(f,ntot              ,'f8').reshape(shape)
        if ndim==3:
            ntot=(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)*nvar
            positions=get_array(f,nx+ny+nz+6*nbuf,'f8')
            self.x  =positions[0           :nx      +2*nbuf]
            self.y  =positions[nx   +2*nbuf:nx+ny   +4*nbuf]
            self.z  =positions[nx+ny+4*nbuf:nx+ny+nz+6*nbuf]
            if not field:
                self.xyz=get_array(f,ntot              ,'f8').reshape(shape)
            else:
                sizevar = ntot/nvar
                startpt = f.tell()
                shape   = shape[1:]
                sbytvar = sizevar*8
                if field == 'rho'  : 
                    self.xyz=seek_array(f, startpt + 0*sbytvar, sizevar, 'f8').reshape(shape)
                if field == 'E'    : 
                    self.xyz=seek_array(f, startpt + 4*sbytvar, sizevar, 'f8').reshape(shape)
                if field == 'rhoux': 
                    self.xyz=seek_array(f, startpt + 1*sbytvar, sizevar, 'f8').reshape(shape)
                if field == 'rhouy': 
                    self.xyz=seek_array(f, startpt + 2*sbytvar, sizevar, 'f8').reshape(shape)
                if field == 'rhouz': 
                    self.xyz=seek_array(f, startpt + 3*sbytvar, sizevar, 'f8').reshape(shape)
                if field == 'Bx'   : 
                    self.xyz=seek_array(f, startpt + 5*sbytvar, sizevar, 'f8').reshape(shape)
                if field == 'By'   : 
                    self.xyz=seek_array(f, startpt + 6*sbytvar, sizevar, 'f8').reshape(shape)
                if field == 'Bz'   : 
                    self.xyz=seek_array(f, startpt + 7*sbytvar, sizevar, 'f8').reshape(shape)

        #Close file
        f.close()

class DumsesFile_h5:

    def __init__(self,filename='slices.000000',nvar=8,nbuf=3,field=None):
        """ Read HDF5 data output from Heracles."""

        #Open file
        f=tables.openFile(filename)

        #Dataset "para_real"
        self.time=f.root.para_real[0]
        self.dt  =f.root.para_real[1]
        self.dx  =f.root.para_real[2]
        self.dy  =f.root.para_real[3]
        self.dz  =f.root.para_real[4]

        #Dataset "para_int"
        self.ndump  =f.root.para_int[0]
        self.nhist  =f.root.para_int[1]
        self.nspec  =f.root.para_int[2]
        self.nx     =f.root.para_int[3]
        self.ny     =f.root.para_int[4]
        self.nz     =f.root.para_int[5]
        self.nxslice=f.root.para_int[6]
        self.nyslice=f.root.para_int[7]
        self.nzslice=f.root.para_int[8]

        self.dim    =f.root.para_int[3:6]
        self.slice  =f.root.para_int[6:9]        
        self.dim_glob=self.dim*self.slice

        #Dataset "x", "y" and "z
        self.x=f.root.x[:]
        self.y=f.root.y[:]
        self.z=f.root.z[:]

        #Dataset "uin"
        if not field:
            self.xyz=f.root.uin[:,:,:,:]
        else:
            if field == 'rho'  : self.xyz = f.root.uin[0,:,:,:]
            if field == 'E'    : self.xyz = f.root.uin[4,:,:,:]
            if field == 'rhoux': self.xyz = f.root.uin[1,:,:,:]
            if field == 'rhouy': self.xyz = f.root.uin[2,:,:,:]
            if field == 'rhouz': self.xyz = f.root.uin[3,:,:,:]
            if field == 'Bx'   : self.xyz = f.root.uin[5,:,:,:]
            if field == 'By'   : self.xyz = f.root.uin[6,:,:,:]
            if field == 'Bz'   : self.xyz = f.root.uin[7,:,:,:]

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

def seek_array(fid, nb_seek, nb_count, type):
    fid.seek(nb_seek)
    pad   = np.fromfile(fid, count=1, dtype='i4')
    array = np.fromfile(fid, count=nb_count, dtype=type)
    return array

def GetSubArray3d(array,ivar,nx,ny,nz,nbuf,ndim):
    if ndim==1:
        SubArray3d=np.transpose(array[ivar,0:nz,0:ny,nbuf:nx+nbuf],(2,1,0))  
    if ndim==2:
        SubArray3d=np.transpose(array[ivar,0:nz,nbuf:ny+nbuf,nbuf:nx+nbuf],(2,1,0))  
    if ndim==3:
        SubArray3d=np.transpose(array[ivar,nbuf:nz+nbuf,nbuf:ny+nbuf,nbuf:nx+nbuf],(2,1,0))  
    
    return SubArray3d

def GetArray3d(array,nx,ny,nz,nbuf,ndim):
    if ndim==1:
        SubArray3d=np.transpose(array[0:nz,0:ny,nbuf:nx+nbuf],(2,1,0))  
    if ndim==2:
        SubArray3d=np.transpose(array[0:nz,nbuf:ny+nbuf,nbuf:nx+nbuf],(2,1,0))  
    if ndim==3:
        SubArray3d=np.transpose(array[nbuf:nz+nbuf,nbuf:ny+nbuf,nbuf:nx+nbuf],(2,1,0))  
    
    return SubArray3d

def GetRoot(D,M,E,gamma,guess):
    #solve Q-P-E=0 with NR to get Q E=E-D Q=Q-D
    x=guess ; epsilon=1.
    while (epsilon>1.e-6):
        root=Getf(x,D,M,E,gamma)
        droot=Getfprime(x,D,M,E,gamma)
        x=x-root/droot
        epsilon=abs(root/droot/x)
    return x

def Getf(x,D,M,E,gamma):
    u2=M**2/((x+D)**2-M**2) ; lor=(1+u2)**(0.5)
    Xsi=(x-u2/(lor+1)*D)/lor**2
    p=(gamma-1.0)/gamma*Xsi
    return x-p-E

def Getfprime(x,D,M,E,gamma):
    u2=M**2/((x+D)**2-M**2) ; lor=(1+u2)**(0.5)
    dpdx=(gamma-1.)/gamma
    dpdx=dpdx*(1.+M**2/(x+D)**2*(1.-D*lor/(x+D)))
    return 1.-dpdx


