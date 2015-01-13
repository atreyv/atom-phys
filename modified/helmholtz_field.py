# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 10:08:58 2014

@author: lab
"""


from scipy.constants import *
from libphys import *
#import mayavi.mlab
#from sympy import *


#define params
x = 0
N_turns = 8
I_cur=  np.array([0.91])   # A
a = 16.5   # c m ,   f u l l   s i d e   l e n g t h
b = 16.5   # c m ,   f u l l   s i d e   l e n g t h
d_H = 9.75   # c m, distance to center
d_aH = 10
mu_0/(4*pi)*10**4*10**2
d_sep = d_H
B_1 = 1/2.*8/5.*N_turns*I_cur*a*b*(a**2+b**2+2*(d_sep-x)**2)/((a**2+b**2)*(b**2+(d_sep-x)**2)*sqrt(a**2 + b**2 + (d_sep-x)**2))
B_2 = 1/2.*8/5.*N_turns*I_cur*a*b*(a**2+b**2+2*(d_sep+x)**2)/((a**2+b**2)*(b**2+ (d_sep + x)**2)*sqrt(a**2 + b**2 + (d_sep + x)**2)) 
B_H=B_1+B_2
B_H

print 'B_H = %f G'%B_H
#xmin, xmax=-1, 1
#plot(B_H.subs({d_sep:d_H}), (x, xmin , xmax))#, linestyle = '-', legend_label = 'H' , axes_labels = [
#'x(cm)',' B(G)'])
# p   + =   p l o t ( B _ a H ( d _ s e p = d _ H ) , ( x , x m i n , x m a x ) ,   l i n e s t y l e   = ' - - ' ,   l e g e n d _ l a b e l = ' a n t i H ' , c o l o r = ' r e d ' )
# s l o p e   =   t a y l o r ( B _ a H ( d _ s e p = d _ H ) , x , 0 , 1 )
# p   + =   p l o t ( s l o p e , ( x , x m i n , x m a x ) ,   l i n e s t y l e = ' - . ' , c o l o r = ' g r e e n ' )
