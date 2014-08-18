from libphys import *

lamb = 780.
dv_fsr = 1.e9
c = 3.e17
dl_fsr = -1 / c * lamb**2 * dv_fsr

filename = 'TDS-FP-20111214-000'
points_start = 13
points_end = 15

middle = 4000
offset = 50
fsr_start = 0
fsr_end = 10000

i1 = [96.78, 96.57, 96.23, 96.04, 95.80, 95.57]
i2 = [93.00, 92.92, 92.71, 92.59, 92.36, 92.27, 92.04, 91.93]
i3 = [89.19, 89.04, 88.94, 88.77, 88.67, 88.55, 88.41, 88.31, 88.21, 87.99, 87.84]
i31= [87.84, 87.68, 87.35, 87.15]

plt.ion()
plt.clf()

fpfile,fp=[],[]
res=[]

#DATA pre-processing
for i in range(points_start, points_end+1):
    if (i>9): filename = 'TDS-FP-20111214-00'
    fpfile += [scp.loadtxt(filename + str(i) + '.dat')]
    fp += [scp.transpose(fpfile[i-points_start][fsr_start:fsr_end])]
    plt.plot(fp[i-points_start][0],fp[i-points_start][1],label=str(i))
    plt.legend()
if (input('1 to continue: ')== 1):
    for i in range(0, points_end + 1 - points_start):
        zero = scp.average(fp[i][1][middle-offset:middle+offset])
        first_peak = scp.amax(fp[i][1][:middle])
        first_peak_pos = scp.argmax(fp[i][1][:middle])
        sec_peak = scp.amax(fp[i][1][middle:])
        sec_peak_pos = scp.argmax(fp[i][1][middle:]) + middle
        first_width_point = scp.argmin(scp.absolute( 0.5*(first_peak-zero) - (fp[i][1][:scp.argmax(fp[i][1][:middle])] - zero)))
        sec_width_point = scp.argmin(scp.absolute(0.5*(first_peak-zero) - (fp[i][1][scp.argmax(fp[i][1][:middle]):middle]-zero))) + scp.argmax(fp[i][1][:middle])
        third_width_point = scp.argmin(scp.absolute(0.5*(sec_peak-zero) - (fp[i][1][middle:middle+scp.argmax(fp[i][1][middle:])] - zero))) + middle
        forth_width_point = scp.argmin(scp.absolute(0.5*(sec_peak-zero) - (fp[i][1][middle+scp.argmax(fp[i][1][middle:]):] - zero))) + middle + scp.argmax(fp[i][1][middle:])
        first_fwhm = scp.absolute( fp[i][0][first_width_point] - fp[i][0][sec_width_point] )
        sec_fwhm = scp.absolute( fp[i][0][third_width_point] - fp[i][0][forth_width_point] )
        fsr = scp.absolute( fp[i][0][first_peak_pos] - fp[i][0][sec_peak_pos] )
        res += [[fsr, fp[i][0][first_peak_pos], first_peak, first_fwhm, fp[i][0][sec_peak_pos], sec_peak, sec_fwhm]]
res = scp.transpose(res)


#DATA Analysis

#First peak shift with descendant current
res[1] = res[1] - res[1][0]
first_shift = []
for i in range(0, points_end - points_start + 1):
    first_shift += [res[1][i] / res[0][i]]
dl, dv = [], []
for i in range(0, points_end - points_start + 1):
    dl += [-1 * lamb**2 / c * dv_fsr * first_shift[i]]
    dv += [dv_fsr * first_shift[i]]
plt.clf()
plt.xlabel('I (mA)')
plt.ylabel('Freq Detuning (Hz)')
plt.scatter(i1,dv)
f = open('TDS-FP-20111214-i1', 'w')
f.write('I(mA)\tDetune(Hz)\n')
for i in range(0, points_end - points_start + 1):
    f.write(str(i1[i])+'\t'+'%e'%dv[i]+'\n')
f.write('# FSR conditions middle, offset from middle, fsr_start and fsr_end\n# ')
f.write(str(middle)+' '+str(offset)+' '+str(fsr_start)+' '+str(fsr_end)+'\n')
f.close()
