# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 14:31:23 2014

@author: pedro
"""

from scipy import ndimage

#Generating Picture for covers
#pump1 = normal
#probe = normal
mayavi.mlab.clf()
#probe3 = probe *.50
pump3 = pump1 *(0.5)# + 20
#probe5=ndimage.gaussian_filter(probe3, sigma=1)
pump5=ndimage.gaussian_filter(pump3, sigma=2.)
#probep = mayavi.mlab.surf(probe5,colormap='OrRd',opacity=1,vmax=np.amax(probe5)*.5)
pumpp = mayavi.mlab.surf(pump5,colormap='jet',representation='surface',line_width=2.,opacity=1,vmax=np.amax(pump3)*0.4)
#mayavi.mlab.savefig('titlepage.tiff', size=(1500,2500), magnification='auto')
#lutprobe = probep.module_manager.scalar_lut_manager.lut.table.to_array()
#lutpump = pumpp.module_manager.scalar_lut_manager.lut.table.to_array()
#lutprobe[:,-1]=np.linspace(0,255,256)
#lutpump[:,-1]=np.linspace(0,255,256)
#probep.module_manager.scalar_lut_manager.lut.table = lutprobe
#pumpp.module_manager.scalar_lut_manager.lut.table = lutpump
#pumpp.module_manager.scalar_lut_manager.reverse_lut = True
mayavi.mlab.draw()
#mayavi.mlab.savefig('researchimages-1.tiff', size=(2500,int(sqrt(2)*2500)), magnification='auto')
