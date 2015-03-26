# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 18:38:53 2014

@author: lab
"""

from libDAQmx import *
import time

parameter_space=np.ogrid[0.:0.:10j,0.695:0.695:10j,0.159:0.159:10j]

for x in np.nditer(parameter_space[0]):
    for y in np.nditer(parameter_space[1]):
        for z in np.nditer(parameter_space[2]):
            # Sequence
#            analog_out("Dev1/ao0", x)
#            analog_out("Dev1/ao4", x) # X
            analog_out("Dev2/ao0", y) # Y
            analog_out("Dev2/ao1", z) # Z
#            analog_out("Dev1/ao0", 7.2)
#            analog_out("Dev1/ao3", 1.7)
            ramp_analog_out("Dev1/ao5", 2.4, 2.4, 1+1)
            time.sleep(3)
#            ramp_analog_out("Dev1/ao0", 7.2, 8., 1+0)
#            ramp_analog_out("Dev1/ao3", 1.7, 1.5, 1+0)
            ramp_analog_out("Dev1/ao5", 2.4, 0., 1+1)
            time.sleep(.5)
#            digital_output("Dev1/port0/line0")
            print x,y,z
#            time.sleep(1.)

