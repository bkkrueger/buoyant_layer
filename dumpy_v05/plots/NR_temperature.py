import numpy as np
from dumpy_v05.plots import cut
from dumpy_v05.plots import thermo
from dumpy_v05.data import rd_dumses
from pylab import *


def NR(x0,alpha,rho):

    x=x0
    alpha=alpha*10**2

    while (abs(x**(5./2.)-x0**4*x**(-3./2.)-alpha*rho) > 10**(-5)):

        df=5./2.*x**(3./2.)+3./2.*x**(-5./2.)*x0**4
        x=-(x**(5./2.)-x0**4*x**(-3./2.)-alpha*rho)/df+x
    
    return x


def fit_temperature(T0,alpha,rho):

    i=0
    a,=rho.shape
    f=np.arange(a)
    while i<a :
#        f[i]=NR(T0,alpha[i],rho[i])
        f[i]=NR(T0,alpha,rho[i])
        i += 1
    return f    
    
def profile(data,save=0,T0=1,gamma=1.00001):
#    x,alpha=cut.GetAlpha(data)
#    alpha=0.01
#    rho=cut.yz_mean(data.rho)
#    f=fit_temperature(T0,alpha,rho)
    t=thermo.pressure(data,gamma=gamma)#/data.rho[:,:,:]*100*gamma
    t=cut.yz_mean(t)-0.01*data.x**(-1.5-1.)/gamma#-data.x**(-1.5)*0.01/1.4
    
#    plot(data.x,f,label='fit')
    plot(data.x,t,label='exp')
    xlabel('radius')
    ylabel('T')
    legend(loc=1)
    #show()

    if save==1:
        savefig('temperature.ps')

    return
