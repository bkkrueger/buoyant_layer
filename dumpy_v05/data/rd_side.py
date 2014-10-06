import tables 
import numpy as np 
import pylab as pl
import glob

class DumsesSides:
    """ Extract data from the sides of DUMSES simulation 
        Usage: data = DumsesSides(ispec, type=type)
        with:
         - ispec: index of the 'special' output
         - type : 'sequential', 'collective_hdf5' or 'collective_mpi' as for the 'store_sides' variable in DUMSES

        data is an object with:
         - data.side[a][i], with a in (x,y,z) and i in (1,2)
        each data.side[a][i] contains the 4 variables: rho, E, rhou & B

        Optional arguments:
         - filedir: path of the directory containing the 'special' directories
         - nbuf   : number of ghost cells"""

    def __init__(self, ispec=1, filedir='./', nvar=8, nbuf=3, type='sequential'):
        self.ispec   = ispec
        self.filedir = filedir
        if(type=='sequential'):
            # Define files to read and retrieve basic informations
            fsidex = fsideSeq(ispec=ispec, filedir=filedir, direction='x')
            fsidey = fsideSeq(ispec=ispec, filedir=filedir, direction='y')
            fsidez = fsideSeq(ispec=ispec, filedir=filedir, direction='z')
            self.npes    = int(np.sqrt(fsidex.nsq*fsidey.nsq*fsidez.nsq))
            self.nxslice = self.npes/fsidex.nsq
            self.nyslice = self.npes/fsidey.nsq
            self.nzslice = self.npes/fsidez.nsq

            nxslice = self.nxslice 
            nyslice = self.nyslice
            nzslice = self.nzslice
            
            # Retrieve color and key values in each direction to organize thread files
            colorx = []; colory = []; colorz = []
            keyx   = []; keyy   = []; keyz   = []
            for mype in range(self.npes):
                colorx.append(mype%nxslice)
                colory.append((mype/(nxslice*nzslice))%nyslice)
                colorz.append((mype/nxslice)%nzslice)
                keyx.append(mype/nxslice)
                keyy.append(mype%(nxslice*nzslice))
                keyz.append(colory[mype]*nxslice \
                                 + (mype-colorz[mype]*nxslice)%nxslice)

            tfx = tables.openFile(fsidex.listtop[0])
            self.ny = tfx.root.rho.shape[1] - 2*nbuf
            self.nz = tfx.root.rho.shape[0] - 2*nbuf
            tfx.close()
            tfy = tables.openFile(fsidey.listtop[0])
            self.nx = tfy.root.rho.shape[1] - 2*nbuf
            tfy.close()
      
            # Create sides in each direction
            self.sidex1 = sideSeq(fsidex.listtop, self.ny, self.nz \
                                   , nyslice, nzslice, keyx, 'x', nbuf)
            self.sidex2 = sideSeq(fsidex.listbot, self.ny, self.nz \
                                   , nyslice, nzslice, keyx, 'x', nbuf)
            self.sidey1 = sideSeq(fsidey.listtop, self.nx, self.nz \
                                   , nxslice, nzslice, keyy, 'y', nbuf)
            self.sidey2 = sideSeq(fsidey.listbot, self.nx, self.nz \
                                   , nxslice, nzslice, keyy, 'y', nbuf)
            self.sidez1 = sideSeq(fsidez.listtop, self.nx, self.ny \
                                   , nxslice, nyslice, keyz, 'z', nbuf)
            self.sidez2 = sideSeq(fsidez.listbot, self.nx, self.ny \
                                   , nxslice, nyslice, keyz, 'z', nbuf)
            
        elif(type=='collective_hdf5'):
            filerep = filedir + 'special_%06d/' %ispec
            fside   = tables.openFile(filerep + 'sides.000000')

            # Retrieve basic informations
            self.nx      = fside.root.para_int[0]
            self.ny      = fside.root.para_int[1]
            self.nz      = fside.root.para_int[2]
            self.nxslice = fside.root.para_int[3]
            self.nyslice = fside.root.para_int[4]
            self.nzslice = fside.root.para_int[5]
            self.dim     = fside.root.para_int[0:3]
            self.slice   = fside.root.para_int[3:6]
            
            nglob = self.dim*self.slice

            # Create sides in each direction
            self.sidex1 = sideColHDF(fside, nglob, side='x', nside=1, nbuf=nbuf)
            self.sidex2 = sideColHDF(fside, nglob, side='x', nside=2, nbuf=nbuf)
            self.sidey1 = sideColHDF(fside, nglob, side='y', nside=1, nbuf=nbuf)
            self.sidey2 = sideColHDF(fside, nglob, side='y', nside=2, nbuf=nbuf)
            self.sidez1 = sideColHDF(fside, nglob, side='z', nside=1, nbuf=nbuf)
            self.sidez2 = sideColHDF(fside, nglob, side='z', nside=2, nbuf=nbuf)

            fside.close()

        elif(type=='collective_mpi'):
            filerep = filedir + 'special_%06d/' %ispec
            fsidex1 = tables.openFile(filerep + 'sides.000001')
            fsidex2 = tables.openFile(filerep + 'sides.000002')
            fsidey1 = tables.openFile(filerep + 'sides.000003')
            fsidey2 = tables.openFile(filerep + 'sides.000004')
            fsidez1 = tables.openFile(filerep + 'sides.000005')
            fsidez2 = tables.openFile(filerep + 'sides.000006')

            # Create sides in each direction
            self.sidex1 = sideColMPI(fsidex1)
            self.sidex2 = sideColMPI(fsidex2)
            self.sidey1 = sideColMPI(fsidey1)
            self.sidey2 = sideColMPI(fsidey2)
            self.sidez1 = sideColMPI(fsidez1)
            self.sidez2 = sideColMPI(fsidez2)

            fsidex1.close(); fsidex2.close()
            fsidey1.close(); fsidey2.close()
            fsidez1.close(); fsidez2.close()

        else:
            print("Wrong type!\n'type' should be equal to 'sequential', 'collective_hdf5' or 'collective_mpi'.")
            return

# Classes for sequential dump of the sides
class fsideSeq:
    def __init__(self, ispec=1, filedir='./', direction='x'):
        # Define directory for the given direction
        filerep = filedir + 'special_%06d/' %ispec
        if(direction=='x'): filerep = filerep + 'sidex/'
        if(direction=='y'): filerep = filerep + 'sidey/'
        if(direction=='z'): filerep = filerep + 'sidez/'

        # Create lists of the files of the sides in the given direction
        self.listtop = np.sort(glob.glob(filerep + 'side*1.*'))
        self.listbot = np.sort(glob.glob(filerep + 'side*2.*'))
        self.nsq     = len(self.listtop)

class sideSeq:
    def __init__(self, listf, nx, ny, nxslice, nyslice, key, direction, nbuf):
        # Initialize arrays
        rho   = np.zeros((nx*nxslice, ny*nyslice))
        rhoux = np.zeros((nx*nxslice, ny*nyslice))
        rhouy = np.zeros((nx*nxslice, ny*nyslice))
        rhouz = np.zeros((nx*nxslice, ny*nyslice))
        E     = np.zeros((nx*nxslice, ny*nyslice))
        Bx    = np.zeros((nx*nxslice, ny*nyslice))
        By    = np.zeros((nx*nxslice, ny*nyslice))
        Bz    = np.zeros((nx*nxslice, ny*nyslice))

        for i in range(nxslice*nyslice):
            mype = int(listf[i].split('.')[-1])
            pos  = key[mype]
            # Index in the global array for each thread file
            if(direction=='x'):
                i0 = (pos/nyslice)*nx; i1 = (pos/nyslice + 1)*nx
                j0 = (pos%nyslice)*ny; j1 = (pos%nyslice + 1)*ny
            else:
                i0 = (pos%nxslice)*nx; i1 = (pos%nxslice + 1)*nx
                j0 = (pos/nxslice)*ny; j1 = (pos/nxslice + 1)*ny

            # Open the file and retrieve data
            f = tables.openFile(listf[i])
            rho[i0:i1,j0:j1] = np.transpose(f.root.rho[nbuf:ny+nbuf,nbuf:nx+nbuf])
            rhoux[i0:i1,j0:j1] = np.transpose(f.root.rho_vx[nbuf:ny+nbuf,nbuf:nx+nbuf])
            rhouy[i0:i1,j0:j1] = np.transpose(f.root.rho_vy[nbuf:ny+nbuf,nbuf:nx+nbuf])
            rhouz[i0:i1,j0:j1] = np.transpose(f.root.rho_vz[nbuf:ny+nbuf,nbuf:nx+nbuf])
            Bx[i0:i1,j0:j1] = np.transpose(f.root.Bx[nbuf:ny+nbuf,nbuf:nx+nbuf])
            By[i0:i1,j0:j1] = np.transpose(f.root.By[nbuf:ny+nbuf,nbuf:nx+nbuf])
            Bz[i0:i1,j0:j1] = np.transpose(f.root.Bz[nbuf:ny+nbuf,nbuf:nx+nbuf])
            E[i0:i1,j0:j1] = np.transpose(f.root.E[nbuf:ny+nbuf,nbuf:nx+nbuf])
            f.close()

        self.rho  = rho
        self.E    = E
        self.rhou = np.transpose(np.array([rhoux,rhouy,rhouz]), (1,2,0))
        self.B    = np.transpose(np.array([Bx,By,Bz]), (1,2,0))

# Classes for collective dump of the sides using PHDF5
class sideColHDF:
    def __init__(self, fside, dim, side='x', nside=1, nbuf=3):
        # Define the selection depending on the selected side
        if(side == 'x'):
            rside = fside.root.side_x
            nxglob = dim[1]
            nyglob = dim[2]
        elif(side == 'y'):
            rside = fside.root.side_y
            nxglob = dim[0]
            nyglob = dim[2]
        elif(side == 'z'):
            rside = fside.root.side_z
            nxglob = dim[0]
            nyglob = dim[1]
        selx = [i+nbuf for i in range(nxglob)]
        sely = [i+nbuf for i in range(nyglob)]

        ndim = 3
        self.rhou = np.zeros((nxglob,nyglob,ndim))
        self.B    = np.zeros((nxglob,nyglob,ndim))

        # Retrieve data
        rho   = np.transpose(rside.rho[nside-1])
        self.rho         = rho[selx,:][:,sely]
        E     = np.transpose(rside.E[nside-1])
        self.E           = E[selx,:][:,sely]
        rhoux = np.transpose(rside.rho_vx[nside-1])
        self.rhou[:,:,0] = rhoux[selx,:][:,sely]
        rhouy = np.transpose(rside.rho_vy[nside-1])
        self.rhou[:,:,1] = rhouy[selx,:][:,sely]
        rhouz = np.transpose(rside.rho_vz[nside-1])
        self.rhou[:,:,2] = rhouz[selx,:][:,sely]
        Bx    = np.transpose(rside.Bx[nside-1])
        self.B[:,:,0]    = Bx[selx,:][:,sely]
        By    = np.transpose(rside.By[nside-1])
        self.B[:,:,1]    = By[selx,:][:,sely]
        Bz    = np.transpose(rside.Bz[nside-1])
        self.B[:,:,2]    = Bz[selx,:][:,sely]

# Classes for collective dump of the sides using MPI
class sideColMPI:
    def __init__(self, fside):
        # Retrieve data
        self.rho  = np.transpose(fside.root.rho[:])
        self.E    = np.transpose(fside.root.E[:])
        rhoux     = np.transpose(fside.root.rho_vx[:])
        rhouy     = np.transpose(fside.root.rho_vy[:])
        rhouz     = np.transpose(fside.root.rho_vz[:])
        self.rhou = np.transpose([rhoux,rhouy,rhouz], (1,2,0))
        Bx        = np.transpose(fside.root.Bx[:])
        By        = np.transpose(fside.root.By[:])
        Bz        = np.transpose(fside.root.Bz[:])
        self.B    = np.transpose([Bx,By,Bz], (1,2,0))
