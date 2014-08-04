# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 17:41:02 2014

@author: Pedro
"""

from PyDAQmx import *
import numpy as np

def analog_in(port):
    analog_input = Task()
    read = int32()
    data = np.zeros((1000,), dtype=np.float64)
    
    # DAQmx Configure Code
    analog_input.CreateAIVoltageChan(port,"",DAQmx_Val_Cfg_Default,-10.0,10.0,DAQmx_Val_Volts,None)
    analog_input.CfgSampClkTiming('',1000.,DAQmx_Val_Rising,DAQmx_Val_FiniteSamps,1000)
    # DAQmx Start Code
    analog_input.StartTask()
    # DAQmx Read Code
    analog_input.ReadAnalogF64(1000,10.0,DAQmx_Val_GroupByChannel,data,1000,byref(read),None)
    # DAQmx stop code    
    analog_input.StopTask()
    # Returns    
    print "Acquired %d points"%read.value
    return data

def analog_out(port,value):
    analog_output = Task()
    write = int32()
    data = np.array([value], dtype=np.float64)
    
    # DAQmx Configure Code
    analog_output.CreateAOVoltageChan(port,"",-10.0,10.0,DAQmx_Val_Volts,None)
#    analog_output.CfgSampClkTiming("",10000.0,DAQmx_Val_Rising,DAQmx_Val_FiniteSamps,1000)
    # DAQmx Start Code
    analog_output.StartTask()
    # DAQmx Read Code
    analog_output.WriteAnalogF64(1,1,10.0,DAQmx_Val_GroupByChannel,data,None,None)
    # DAQmx Stop Code
    analog_output.StopTask()
    # Returns
    print "Set voltage %f"%data
