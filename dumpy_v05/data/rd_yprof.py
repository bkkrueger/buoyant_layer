import tables

class zData:

    def __init__(self,filename='zprof_000001.h5'):
        """ Read HDF5 data output from Heracles."""

        #Open file
        f=tables.openFile(filename)

        #Dataset "para_real"
        self.time=f.root.para_real[0]

        #Dataset "yglob"
        self.y=f.root.yglob[:]

        #Dataset "yuin"
        self.uin=f.root.yuin[:,:]

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
        return slice
