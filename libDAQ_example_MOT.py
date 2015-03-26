# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 18:16:19 2014

@author: lab
"""

from libDAQmx import *
import time
# Sequence

# Set the frequency for trap beams (7.2 V for 2.7 gamma detuning)
analog_out("Dev1/ao1", 1.61)

# Set the intensity for trap beams (5.5 for maximum intensity)
analog_out("Dev1/ao0", 5.5)

# Set the intensity for rep beams (saturation value at 3.8 V)
analog_out("Dev1/ao2", 3.8)

# Set the frequency for trap beams
#analog_out("Dev1/ao3", 1.7)

# Set the frequency for repump beams
analog_out("Dev1/ao7", 1.475)

#probe int
analog_out("Dev1/ao3", 0)

#probe det
analog_out("Dev1/ao4", 2.65)

# Set the compensation uniform B field
#turn_on_compensation_coils()
#turn_off_compensation_coils()
#analog_out("Dev1/ao7", 0.) # X
analog_out("Dev2/ao0", 0.695) # Y
analog_out("Dev2/ao1", 0.159) # Z

# Set the gradient B field
#ramp_analog_out("Dev2/ao0", 2., 2., 1+1)
#time.sleep(3.)
ramp_analog_out("Dev1/ao5", 2.4, 2.4, 1+1)
time.sleep(.4)
#digital_output("Dev1/port0/line0")

# Read from photodiode
pd = analog_in("Dev2/ai0")
print pd.mean()