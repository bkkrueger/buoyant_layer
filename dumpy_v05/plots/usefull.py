from pylab import *
import numpy as np
from dumpy_v05.plots import history

def alpha_compare(**repertoires):
    for cle,valeur in repertoires.items():
        data=history.Hist(filename='./'+valeur+'/history.txt')
        plot(data.time/(2.*pi),data.alpha, label=str(cle))
        xlabel('Time (in orbits)')
        ylabel('Alpha')
    legend(loc=2)
