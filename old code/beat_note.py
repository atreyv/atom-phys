from libphys import *
from scipy.optimize import leastsq as spleastsq

def lorentz(x,p):
    """Calculates the Lorentz distribution for a given set"""
    return p[0]/(1 + pow(2*(x-p[1])/p[2],2)) + p[3]

def residuals_lorentz(p,y,x):
    """Calculates the array of residuals from Lorentz distribution"""
    lmax, lx0, lfw, loffset = p
    err = y - lmax/(1 + pow(2*(x-lx0)/lfw,2)) - loffset
    return err

dname = '/home/pedro/Downloads/LAB/DATA/2012/Beatnote/'

ion()


fwhm, pos = [], []

i=-1
j=-1
while (i!=0):
    if(i==-1):
        i=input('positive to move forward, negative to backwards, anything else to quit\n')
    j=j+i
    if(j<10):
        f=sp.loadtxt(dname+'ms140212000'+str(j)+'.dat',skiprows=6)
    elif(j>=10 and j<100):
        f=sp.loadtxt(dname+'ms14021200'+str(j)+'.dat',skiprows=6)
    else:
        f=sp.loadtxt(dname+'ms1402120'+str(j)+'.dat',skiprows=6)
    f=sp.transpose(f)
    f[0]/=1000000.
    #uncomment for log scale
#    f[1]=pow(10,f[1]/10.)
    plt.clf()
    plt.xlabel('Beat note (MHz)')
    plt.ylabel('')
    lmax=sp.amax(f[1])
    largmax=f[0][sp.argmax(f[1])]
    loffset=sp.amin(f[1])
    a=spleastsq(residuals_lorentz,[lmax,largmax,1,loffset],args=(f[1],f[0]))
    l=[]
    for k in range(0,len(f[0])):
        l+=[lorentz(f[0][k],a[0])]
    plt.scatter(f[0],f[1])
    plt.plot(f[0],l)
    fwhm += [a[0][2]]
    pos += [a[0][1]]
    print 'FWHM = '+str(a[0][2])
    print 'Pos = '+str(a[0][1])
    print 'Current file is '+str(j)
    waitforbuttonpress(timeout=.5) #shitty hack!....
    i=input('positive to move forward, negative to backwards, anything else to quit\n')

print '<FWHM> = ' + str(sp.average(fwhm))
print 'stdev_FWHM = ' + str(sp.std(fwhm))
print 'POSmax - POSmin = ' + str(sp.amax(pos)-sp.amin(pos))
print 'stdev_POS = ' + str(sp.std(pos))
