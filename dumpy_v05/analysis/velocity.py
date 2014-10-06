from numpy import *
from dumpy_v05.data import rd_dumses

def vel_disp(istart,iend):
    c0=1.e-3
    dvx=0. ; dvy=0. ; dvz=0.
    for iloop in range(iend-istart+1):
	print "dump file number:",istart+iloop
	idump=istart+iloop
	data=rd_dumses.read_dumses_mpi(idump)
	if iloop==0:
	    (nvar,nz,ny,nx)=shape(data)
	    ncell=nx*ny*nz
	dvx=(dvx*iloop+sum(data[1,:,:,:]**2)/ncell)/(iloop+1)
	dvz=(dvz*iloop+sum(data[2,:,:,:]**2)/ncell)/(iloop+1)
	dvy=(dvy*iloop+sum(data[3,:,:,:]**2)/ncell)/(iloop+1)
	print "     ",sqrt(dvx)/c0,sqrt(dvy)/c0,sqrt(dvz)/c0
    print "Fluctuations:"
    print "    vx=",sqrt(dvx)/c0
    print "    vy=",sqrt(dvy)/c0
    print "    vz=",sqrt(dvz)/c0
