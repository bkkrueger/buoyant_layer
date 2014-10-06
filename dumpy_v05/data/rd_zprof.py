import tables

class zData:

    def __init__(self,filename='zprof_000001.h5'):
        """ Read HDF5 data output from Heracles."""

        #Open file
        f=tables.openFile(filename)

        #Dataset "para_real"
        self.time=f.root.para_real[0]

        #Dataset "zglob"
        self.z=f.root.zglob[:]

        #Dataset "zuin"
        self.uin=f.root.zuin[:,:]

        #Dataset "rhovxvy"
        self.rhovxvy=f.root.rhovxvy[:]

        #Dataset "rhovx"
        self.rhovx=f.root.rhovx[:]

        #Dataset "rhovy"
        self.rhovy=f.root.rhovy[:]

        #Dataset "rhovz"
        self.rhovz=f.root.rhovz[:]

        #Dataset "maxwell"
        self.maxwell=f.root.maxwell[:]

        #Close file
        f.close()

    def get_array(self,type='rho'):
        if type=='rho':
            slice=self.uin[0,:]
        elif type=='vx':
            slice=self.uin[1,:]/self.uin[0,:]
        elif type=='vy':
            slice=self.uin[2,:]/self.uin[0,:]
        elif type=='vz':
            slice=self.uin[3,:]/self.uin[0,:]
        elif type=='E':
            slice=self.uin[4,:]
        elif type=='bx':
            slice=self.uin[5,:]
        elif type=='by':
            slice=self.uin[6,:]
        elif type=='bz':
            slice=self.uin[7,:]
        elif type=='alfven':
            slice=self.uin[7,:]/(self.uin[0,:])**(0.5)
        elif type=='maxwell':
            slice=self.maxwell[:]
        elif type=='reynolds':
            slice=self.rhovxvy[:]
        elif type=='rhovx':
            slice=self.rhovx[:]
        elif type=='rhovy':
            slice=self.rhovy[:]
        elif type=='rhovz':
            slice=self.rhovz[:]
        return slice
