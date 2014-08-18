# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 17:36:14 2014

@author: pedro
"""

from libphys import *

def realtimeplot():
    figure()
#    pause(1e-2)
#    x, y = rand(), rand()
    ax1, = plot([0,1],[0,1])
    x = array([])
    y = array([])
    for i in range(0,int(1e4)):
        x = append(x,rand())
        y = append(y,rand())
        ax1.set_data(x,y)
#        draw()
        pause(1e-6)
#        waitforbuttonpress(1e-2)
#        draw()
#    return 0

realtimeplot()
