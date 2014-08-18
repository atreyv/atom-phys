# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 01:41:16 2014

@author: pedro
"""
from sympy import *

x=symbols('x')
mod_max = 4*pi
pmax = 4*pi
pmin = -pi

#input field
#p = line([(0,-1),(mod_max,-1)],color='blue', linestyle='-')
#p+= line([(0,-1.3),(mod_max,-1.3)],color='blue', linestyle='--')
p= arrow((mod_max/2,-pi),(mod_max/2,-pi+1.5), color = 'blue', arrowsize=6, width=4)

#length scale dimension
p+= arrow((2*pi-1,-0.3),(pi,-0.3),color='black',width=0.5, arrowsize=3)
p+= arrow((2*pi+1,-0.3),(3*pi,-0.3),color='black',width=0.5, arrowsize=3)
p+= text('$\Lambda$',(mod_max/2,-0.4), color='black', rotation='vertical', fontsize=20)

#medium
#p+= line([(0,0),(mod_max,0)],color='black')
#p+= line([(0,2),(mod_max,2)],color='black')
#p+= line([(0,2),(0,0)],color='black')
#p+= line([(mod_max,2),(mod_max,0)],color='black')

#axis x
p+= arrow((0,0),(mod_max*1.2,0),color='black',width=0.5)
p+= text('$x$',(mod_max*1.2,-.5), color='black', rotation='vertical', fontsize=20)

#medium density and refractive index
p+= plot(0.1*cos(x)+1,(x,0,mod_max), axes=False, ymax=pmax,ymin=pmin, linestyle='-', color='green', legend_label='n(x)')
#p+= plot(cos(x+pi)+1,(x,0,mod_max), linestyle='-', color='red', legend_label='density grating')

#PM field
p+= plot(cos(x)+3,(x,0,mod_max), linestyle='--', color='blue', legend_label='Phase mod')
p+= arrow((pi,2),(pi,2.5),color='blue',arrowsize=5,width=1)
p+= arrow((3*pi,2),(3*pi,2.5),color='blue',arrowsize=5,width=1)

#feedback AM field 
p+= plot(cos(x+pi)+3,(x,0,mod_max), linestyle='-', color='blue', legend_label='Amp mod')
p+= arrow((pi,4),(pi,3.5),color='blue',arrowsize=5,width=1)
p+= arrow((3*pi,4),(3*pi,3.5),color='blue',arrowsize=5,width=1)

#Main beam and sidebands
p+= arrow((mod_max/2,2),(mod_max/2+.5,2+1.5),color='blue',width=2,arrowsize=5)
p+= arrow((mod_max/2,2),(mod_max/2,2+sqrt(2.5)),color='blue',arrowsize=5, width = 3)
p+= arrow((mod_max/2,2),(mod_max/2-.5,2+1.5),color='blue',width=2,arrowsize=5)

#Mirror
#p+= line([(0,pmax/2+1),(mod_max,pmax/2+1)],color='black')
#p+= sum([line([(1.2*y,pmax/2+1),(1.2*y+.5,pmax/2+1.5)], color='black') for y in [0..mod_max/1.2]])

#AM field
p+= plot(cos(x)+pmax-1,(x,0,mod_max), linestyle='-', color='blue')
p+= arrow((pi,pmax-2),(pi,pmax-1.5),color='blue',arrowsize=5,width=1)
p+= arrow((3*pi,pmax-2),(3*pi,pmax-1.5),color='blue',arrowsize=5,width=1)

#d dimensions
p+= arrow((mod_max+1, pmax/2+1),(mod_max+1,2), color='black', width=0.5, arrowsize=3)
p+= arrow((mod_max+1, 2),(mod_max+1,pmax/2+1), color='black', width=0.5, arrowsize=3)
p+= text('$d$',(mod_max+1.5,pmax-pmax/4+.5), color='black', rotation='vertical', fontsize=20)

p+= arrow((mod_max+1, pmax/2+1),(mod_max+1,pmax), color='black', width=0.5, arrowsize=3)
p+= arrow((mod_max+1, pmax),(mod_max+1,pmax/2+1), color='black', width=0.5, arrowsize=3)
p+= text('$d$',(mod_max+1.5,2+pmax/4-.5), color='black', rotation='vertical', fontsize=20)

#L dimension
#p+= arrow((mod_max+1, 0),(mod_max+1,2), color='black', width=0.5, arrowsize=3)
#p+= arrow((mod_max+1, 2),(mod_max+1,0), color='black', width=0.5, arrowsize=3)
#p+= text('$L$',(mod_max+1.5,1), color='black', rotation='horizontal', fontsize=20)

p.save('/home/pedro/Dropbox/PhD at Strathclyde/PF drafts/figures/talbot.png',aspect_ratio=1,dpi=150,show_legend=False)
p.show(aspect_ratio=1,dpi=150,show_legend=False)