#(1010,621)
import urllib.request as url
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import os
from pathlib import Path
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from skimage import exposure
import skimage
import skimage.morphology as morph
from skimage.feature import peak_local_max
from skimage import data, img_as_float
from scipy import ndimage as ndi
from matplotlib.colors import LogNorm
#import the files in the directory as a list
d1 = Path('26Proserpina/26Proserpina/20120125/Day 4')
day1 = os.listdir(d1)
headerd1 = fits.open(d1/ day1[0])

#getting data and median for the Proserpina data
for i in range(0,int(len(day1))):
    zd1 = fits.getdata(d1/ day1[i])
    if i == 0:
        yd1 = zd1
    else:
        yd1 = np.dstack([yd1,zd1])
xd1 = yd1.transpose(2,0,1)
xd1 = np.median(xd1,axis=0)

img = fits.getdata(d1/ day1[0])

img = fits.getdata(d1/ day1[0])
img = img[::-1]
plt.imshow(img,cmap='gray_r',origin = "lower", norm=LogNorm(), vmin=3)
plt.savefig("d4_1.pdf",bbox_inches = 'tight')
plt.show()
plt.imshow(img,cmap='gray',origin = "lower", norm=LogNorm(), vmin=50)


#import the files in the directory as a list for darks
d1_dark = Path('26Proserpina/26Proserpina/20120125/Dark day 4')
day1_dark = os.listdir(d1_dark)

#getting data and median for the dark data
for i in range(0,int(len(day1_dark))):
    zd1_dark = fits.getdata(d1_dark/ day1_dark[i])
    if len(zd1_dark) == 2048:
        
        if i == 0:
            yd1_dark = zd1_dark
        else:
            yd1_dark = np.dstack([yd1_dark,zd1_dark])
xd1_dark = yd1_dark.transpose(2,0,1)
xd1_dark = np.median(xd1_dark,axis=0)



#import the files in the directory as a list for bias
d1_bias = Path('26Proserpina/26Proserpina/20120125/Bias day 4')
day1_bias = os.listdir(d1_bias)

#getting data and median for the bias data
for i in range(0,int(len(day1_bias))):
    zd1_bias = fits.getdata(d1_bias/ day1_bias[i])
    if i == 0:
        yd1_bias = zd1_bias
    else:
        yd1_bias = np.dstack([yd1_bias,zd1_bias])
xd1_bias = yd1_bias.transpose(2,0,1)
xd1_bias = np.median(xd1_bias,axis=0)

 
 
#import the files in the directory as a list for flat
d1_flat = Path('26Proserpina/26Proserpina/20120125/Flats day 4')
day1_flat = os.listdir(d1_flat)

#getting data and median for the flat data
for i in range(0,int(len(day1_flat))):
    zd1_flat = fits.getdata(d1_flat/ day1_flat[i])
    if i == 0:
        yd1_flat = zd1_flat
    else:
        yd1_flat = np.dstack([yd1_flat,zd1_flat])
xd1_flat = yd1_flat.transpose(2,0,1)
xd1_flat = np.median(xd1_flat,axis=0) 



#subtracting bias from dark and then the value from raw data
dd1 = np.abs(np.subtract(xd1,xd1_dark))
dd2 = np.abs(np.subtract(xd1_flat,xd1_bias))
#dividing the flat from the numerator
day1data = np.divide(dd1,dd2) 
day1data = day1data[::-1]
#non stacked data
raw = fits.getdata(d1/ day1[0])
rawdark = fits.getdata(d1_dark/ day1_dark[10])
rawbias = fits.getdata(d1_bias/ day1_bias[0])
rawflat = fits.getdata(d1_flat/ day1_flat[0])



listpeak =[]
listpeaks =[]
import skimage
from skimage.feature import peak_local_max
peaks = peak_local_max(day1data,indices=True,num_peaks = 100) #locate the peaks of intensity
a = np.transpose(peaks)
#plot the rows and columns of pixels
row = a[1]
col = a[0]
plt.axis((0,2045,0,2045))
plt.title("26/01/2012")
plt.xlabel("Pixel")
plt.ylabel("Pixel")
plt.axis((0,2045,0,2045))
#find the centroids
centroids_xpos = []
centroids_ypos = []
for i in range(100):
    (x, y) = (row[i], col[i])
    numeratorx = 0
    numeratory = 0
    total_intensity = 0
#create a 6 pixel long and high square around 
    for j in range(-3, 4): 
        for k in range(-3, 4):
            (xi, yi) = (x+k, y+j)
#do not forget to exclude the last square as it will be out of range
            if (xi > 2045) or (yi > 2045):
                pass
            if (xi == 0) or (yi == 0):
                pass
            else:
#use the equations of the centroids
                curr_intensity = day1data[xi][yi]
                total_intensity += curr_intensity
                numeratorx += xi*curr_intensity
                numeratory += yi*curr_intensity
    xpos = numeratorx / total_intensity
    ypos = numeratory / total_intensity
    centroids_xpos.append(xpos)
    centroids_ypos.append(ypos)
#plot the centroids
plt.scatter(centroids_xpos, centroids_ypos, s=50,facecolor='none',edgecolor='orange')
plt.savefig("d4_2.pdf",bbox_inches = 'tight')
plt.show()


#find the USNO stars
def usno(radeg,decdeg,fovam): # RA/Dec in decimal degrees/J2000.0 FOV in arc min. 
    
    #url format for USNO
    str1 = 'http://webviz.u-strasbg.fr/viz-bin/asu-tsv/?-source=USNO-B1' 
    str2 = '&-c.ra={:4.6f}&-c.dec={:4.6f}&-c.bm={:4.7f}/{:4.7f}&-out.max=unlimited'.format(radeg,decdeg,fovam,fovam)
    #final URL: make sure it does not have any spaces or carriage returns/line feeds when copy-pasting
    sr = str1+str2 

    # Read from the webpage, storing the page's contents in 's'. 
    f = url.urlopen(sr)
    s = f.read()
    f.close()

    #column interpretation compass
    namecol, RAcol, DECcol, rband = 0, 1, 2, 12
    null1, null2 = '     ', ''
    
    #split webpage string into lines
    sl = s.splitlines() 
    sl = sl[45:-1] # get rid of header 
    name = np.array([]) 
    rad = np.array([]) # RA in degrees 
    ded = np.array([]) # DEC in degrees 
    rmag = np.array([]) # rmage 
    #get data from each line
    for k in sl:
        kw = k.decode().split('\t')
        if kw[0] != '':  
            name = np.append(name,kw[namecol])
            rad = np.append(rad,float(kw[RAcol])) 
            ded = np.append(ded,float(kw[DECcol]))
            # deal with case where no mag is reported
            if (kw[rband] != null1) and (kw[rband] != null2):
                rmag = np.append(rmag,float(kw[rband]))    
            else:
                rmag = np.append(rmag,np.nan)
    #return data
           
    return name,rad,ded,rmag

#main function
if __name__ == '__main__':

    #get stars in a 5'x5' square around (RA, DEC = 30, 30)
    ras = headerd1[0].header['RA']
    des = headerd1[0].header['DEC']
    radeg = 15*(float(ras[0:2])+float(ras[3:5])/60.+float(ras[6:])/3600)
    dsgn = np.sign(float(des[0:3]))
    dedeg = float(des[0:3])+dsgn*float(des[4:6])/60.+dsgn*float(des[7:])/3600.
    fovam = 37
    name,rad,ded,rmag = usno(radeg,dedeg,fovam)
    print(rmag) #there's missing mags, but you don't need mags anyways.
    #mask determines the rmag maximum so the image won't be too crowded.
    mask = np.where(rmag<15)[0]
    rad = rad[mask]
    ded = ded[mask]
    
    #plot positions
    plt.title("USNO-B1 Stars")
    plt.scatter(rad*100, ded*100, marker='x',c='red')
    plt.xlabel("RA [deg]")
    plt.ylabel("DEC [deg]")
    plt.tight_layout()
    plt.savefig("d4_3.pdf",bbox_inches = 'tight')
    plt.show()



alist = []
#convert the units to radians
rad_radians = np.radians(rad)
ded_radians = np.radians(ded)
alphanot = np.radians(radeg)
deltanot = np.radians(dedeg)
p = 0.018
X_values = []
#first equation apply on the RA and DEC to pixels of x and y
for i,j in zip(rad_radians,ded_radians):
        numX_dividedby_denX=-(np.divide((np.cos(j)*np.sin(i-alphanot)),(np.cos(deltanot)*np.cos(j)*np.cos(i-alphanot) + np.sin(j)*np.sin(deltanot))))
        X_values.append(numX_dividedby_denX)
       
X_values_array = np.asarray(X_values)
x_pixels = 3454*(X_values_array/p)+1024

Y_values = []
for i,j in zip(rad_radians,ded_radians):
        numY_dividedby_denY=-(np.divide((np.sin(deltanot)*np.cos(j)*np.cos(i-alphanot)-np.cos(deltanot)*np.sin(j)),(np.cos(deltanot)*np.cos(j)*np.cos(i-alphanot) + np.sin(j)*np.sin(deltanot))))
        Y_values.append(numY_dividedby_denY)
       
Y_values_array = np.asarray(Y_values)
y_pixels = 3454*(Y_values_array/p)+1024
#plot the USNO pixel positions
plt.title("Pixel Positions of USNO-B1 Stars")
plt.scatter(x_pixels, y_pixels, marker='x',c='black')
plt.savefig("d4_4.pdf",bbox_inches = 'tight')
plt.show()


#overlay
plt.scatter(centroids_xpos, centroids_ypos, marker='+',c='blue',label='CCD')
plt.scatter(x_pixels, y_pixels, marker='x',c='red',label='USNO')
plt.legend()
plt.savefig("d4_5.pdf",bbox_inches = 'tight')
plt.show()


x_star=[]
y_star=[]
xusno = []
yusno= []
xoffset =[]
yoffset =[]
#finding the separation between the points selected with distance lower than 20 pixels
length_rad = len(rad)
#range is 100 based on number of points of CCD image
for i in range(100):
    raccd = centroids_xpos[i]
    deccd = centroids_ypos[i]
#range is determined by the number of USNO points
    for j in range(len(rad)):
        rausno = x_pixels[j]
        decusno = y_pixels[j]
#separation is the distance between the points
        separation = np.sqrt(((raccd-rausno)**2)+((deccd-decusno)**2)) 
        if separation < 20:
            x_star.append(raccd)
            y_star.append(deccd)
            xusno.append(rausno)
            yusno.append(decusno)
            xoffset.append(raccd-rausno)
            yoffset.append(deccd-decusno)

xusno = np.array(xusno)
yusno=np.array(yusno)
X_values = (((xusno-1024))*0.018)/3454
Y_values = (((yusno-1024))*0.018)/3454     
X_values = np.array(X_values)
Y_values = np.array(Y_values)

#plotting the callibrated points that match with a maximum offset of 20
plt.plot(xusno,yusno,'o',c='red', marker='x',label='USNO')
plt.plot(x_star,y_star,'o',c='blue',marker='+',label='CCD')
plt.legend(loc='best')
plt.legend()
plt.savefig("d4_6.pdf",bbox_inches = 'tight')
plt.show()
plt.figure(figsize=(9,4))
#plotting the residuals or the distance between the callibrated points
plt.title("Residuals")
plt.xlabel("x or y [pixel]")
plt.ylabel("Pixel offset")
plt.scatter(xusno, xoffset,c='red', marker='v',label='   x')
plt.scatter(yusno, yoffset,c='blue', marker='^',label='   y')
plt.legend(loc='best')
plt.savefig("d4_7.pdf",bbox_inches = 'tight')
plt.show()


#T1
fp = 3454.00/0.018 
ax = np.zeros(len(X_values))
#since the shapes don't match, a loop has to be done
for i in range(len(x_star)):
    ax[i] = x_star[i]
#finding the B matrix
#B is already transposed according to its shape (matrix_B.shape)
matrix_Bt = np.array([fp*X_values,fp*Y_values,np.ones(len(X_values))])
matrix_Bt_ax = np.dot(matrix_Bt,ax)
matrix_B = np.transpose(matrix_Bt)
Bt_times_B = np.dot(matrix_Bt,matrix_B)
Bt_B_inverse = np.linalg.inv(Bt_times_B)
cx = np.dot(Bt_B_inverse,matrix_Bt_ax)
#for y
ay = np.zeros(len(Y_values))
for i in range(len(y_star)):
    ay[i] = y_star[i]
#finding the B matrix
matrix_Bty= np.array([fp*X_values,fp*Y_values,np.ones(len(X_values))])
matrix_Bt_ay = np.dot(matrix_Bty,ay)
matrix_By= np.transpose(matrix_Bty)
Bt_times_By= np.dot(matrix_Bty,matrix_By)
Bt_B_inversey= np.linalg.inv(Bt_times_By)
cy = np.dot(Bt_B_inversey,matrix_Bt_ay)
#plate constants
print('Plate Constants')
print('a_11, a_12, x_0 = ',cx)
print('a_21, a_22, y_0 = ',cy)
#T matrix
T = [[fp*(cx[0]), fp*(cx[1]) , cx[2]], [fp*(cy[0]) , fp*(cy[1]) , cy[2]], [0, 0, 1]]
print('T Matrix = ',T)

#find f from T1 using:
fp2=(np.sqrt(np.linalg.det(T)))
fp3=fp2
#3217.724312443244
#big 2
#checking the values through
import math
"""
bigx =(( 1037-1024)/3454)*p
bigy =(( 991-1024)/3454)*p
tann = -(bigx/((np.cos(deltanot))-(bigy*np.sin(deltanot))))
tann = math.atan(tann)
alpha=tann+ alphanot
"""

"""
#T2

matrix_Bt2= np.array([fp2*X_values_array,fp2*Y_values_array,ones])
matrix_Bt_ax2= np.dot(matrix_Bt2,ax)
matrix_B2= np.transpose(matrix_Bt2)
Bt_times_B2= np.dot(matrix_Bt2,matrix_B2)
Bt_B_inverse2= np.linalg.inv(Bt_times_B2)
cx2= np.dot(Bt_B_inverse2,matrix_Bt_ax2)


matrix_Bty2= np.array([fp2*X_values_array,fp2*Y_values_array,ones])
matrix_Bt_ay2= np.dot(matrix_Bty2,ay)
matrix_By2=np.transpose(matrix_Bty2)
Bt_times_By2=np.dot(matrix_Bty2,matrix_By2)
Bt_B_inversey2=np.linalg.inv(Bt_times_By2)
cy2= np.dot(Bt_B_inversey2,matrix_Bt_ay2)

print('Plate Constants')
print('a_11, a_12, x_0 = ',cx2)
print('a_21, a_22, y_0 = ',cy2)

#T matrix
T2= [[fp*cx2[0] , fp*cx2[1] , cx2[2]], [fp*cy2[0] , fp*cy2[1] , cy2[2]], [0, 0, 1]]
print('T Matrix = ',T2)
"""
#T2
matrix_Bt3= np.array([fp3*X_values,fp3*Y_values,np.ones(len(Y_values))])
matrix_Bt_ax3= np.dot(matrix_Bt3,ax)
matrix_B3= np.transpose(matrix_Bt3)
Bt_times_B3= np.dot(matrix_Bt3,matrix_B3)
Bt_B_inverse3= np.linalg.inv(Bt_times_B3)
cx3= np.dot(Bt_B_inverse3,matrix_Bt_ax3)

matrix_Bty3= np.array([fp3*X_values,fp3*Y_values,np.ones(len(Y_values))])
matrix_Bt_ay3= np.dot(matrix_Bty3,ay)
matrix_By3=np.transpose(matrix_Bty3)
Bt_times_By3=np.dot(matrix_Bty3,matrix_By3)
Bt_B_inversey3=np.linalg.inv(Bt_times_By3)
cy3= np.dot(Bt_B_inversey3,matrix_Bt_ay3)

print('Plate Constants')
print('a_11, a_12, x_0 = ',cx3)
print('a_21, a_22, y_0 = ',cy3)

#T matrix
T3= [[fp*(cx3[0]) , fp*(cx3[1]) , cx3[2]], [fp*(cy3[0]) , fp*(cy3[1]) , cy3[2]], [0, 0, 1]]
print('T Matrix = ',T3)
T_inv = np.linalg.inv(T)
#finding the big X and big Y using T2
T3_inv=np.linalg.inv(T3)
#finding residuals after using T matrix to multiply and find a new set of standard coordinates
"""
position=[]
arrayy = np.array([centroids_xpos,centroids_ypos,ones])
arrayy = np.transpose(arrayy)
for star in arrayy:
    okk = np.dot(T_inv,star)
    position.append(okk)
lit = np.array(position)
"""
arrayy=[]
okkk=[]
okkk3=[]
x_star=np.array(x_star)
y_star=np.array(y_star)
arrayy = np.array([x_star,y_star,np.ones(len(x_star))])
arrayy = np.transpose(arrayy)
for i in arrayy:
    okk = np.dot(T_inv,i)
    okkk.append(okk)
okkk = np.transpose(okkk)
xCCD_residual = okkk[0]
yCCD_residual = okkk[1]

for i in arrayy:
    okk3 = np.dot(T3_inv,i)
    okkk3.append(okk3)
okkk3 = np.transpose(okkk3)
xCCD_residual3 = okkk3[0]
yCCD_residual3 = okkk3[1]

xCCD=[]
yCCD=[]
xUSNO = []
yUSNO= []
xOFFSET =[]
yOFFSET =[]
raccdx=[]
deccdy=[]
rausnox=[]
decusnoy=[]

xCCD3=[]
yCCD3=[]
xUSNO3 = []
yUSNO3= []
xOFFSET3 =[]
yOFFSET3 =[]
raccdx3=[]
deccdy3=[]
rausnox3=[]
decusnoy3=[]

xCCDpix = 3454*(xCCD_residual/p)+1024
yCCDpix = 3454*(yCCD_residual/p)+1024
xUSNOpix = 3454*(X_values/p)+1024
yUSNOpix = 3454*(Y_values/p)+1024
xCCDpix = np.array(xCCDpix)
yCCDpix = np.array(yCCDpix)
xUSNOpix = np.array(xUSNOpix)
yUSNOpix = np.array(yUSNOpix)

xCCDpix3 = 3454*(xCCD_residual3/p)+1024
yCCDpix3 = 3454*(yCCD_residual3/p)+1024
xUSNOpix3 = 3454*(X_values/p)+1024
yUSNOpix3 = 3454*(Y_values/p)+1024
xCCDpix3 = np.array(xCCDpix3)
yCCDpix3 = np.array(yCCDpix3)
xUSNOpix3 = np.array(xUSNOpix3)
yUSNOpix3  = np.array(yUSNOpix3)

#range is 100 based on number of points of CCD image
for i in range(len(xCCD_residual)):
    raccdx = xCCDpix[i]
    deccdy = yCCDpix[i]
#range is determined by the number of USNO points
    for j in range(len(X_values)):
        rausnox = xUSNOpix[j]
        decusnoy = yUSNOpix[j]
#separation is the distance between the points
        sep = np.sqrt(((raccdx-rausnox)**2)+((deccdy-decusnoy)**2)) 
        if sep < 40:
            xCCD.append(raccdx)
            yCCD.append(deccdy)
            xUSNO.append(rausnox)
            yUSNO.append(decusnoy)
            xOFFSET.append(raccdx-rausnox)
            yOFFSET.append(deccdy-decusnoy)

xCCD = np.array(xCCD)
yCCD = np.array(yCCD)
xUSNO = np.array(xUSNO)
yUSNO = np.array(yUSNO)
xOFFSET = np.array(xOFFSET)
yOFFSET = np.array(yOFFSET)


#range is 100 based on number of points of CCD image
for i in range(len(xCCD_residual3)):
    raccdx3 = xCCDpix3[i]
    deccdy3 = yCCDpix3[i]
#range is determined by the number of USNO points
    for j in range(len(X_values)):
        rausnox3 = xUSNOpix3[j]
        decusnoy3 = yUSNOpix3[j]
#separation is the distance between the points
        sep3 = np.sqrt(((raccdx3-rausnox3)**2)+((deccdy3-decusnoy3)**2)) 
        if sep3 < 10:
            xCCD3.append(raccdx3)
            yCCD3.append(deccdy3)
            xUSNO3.append(rausnox3)
            yUSNO3.append(decusnoy3)
            xOFFSET3.append(raccdx3-rausnox3)
            yOFFSET3.append(deccdy3-decusnoy3)

xCCD3 = np.array(xCCD3)
yCCD3 = np.array(yCCD3)
xUSNO3 = np.array(xUSNO3)
yUSNO3 = np.array(yUSNO3)
xOFFSET3 = np.array(xOFFSET3)
yOFFSET3 = np.array(yOFFSET3)


#plotting the residuals or the distance between the callibrated points
plt.figure(figsize=(9,4))
#plt.scatter(xUSNO, xOFFSET,c='red', marker='v',label='   x')
#plt.scatter(yUSNO, yOFFSET,c='blue', marker='^',label='   y')
plt.scatter(xCCD, xOFFSET,c='red', marker='v')
plt.scatter(yCCD, yOFFSET,c='blue', marker='^')
plt.title("Residuals after T1")
plt.xlabel("x or y [pixel]")
plt.ylabel("Pixel offset")
plt.legend()
#plt.axis((0,2045,-5,5))
plt.savefig("d4_8",bbox_inches = 'tight')
plt.show()

plt.figure(figsize=(9,4))
plt.title("Residuals after T2")
plt.xlabel("x or y [pixel]")
plt.ylabel("Pixel offset")
#plt.scatter(xUSNO3, xOFFSET3,c='red', marker='v',label='   x')
#plt.scatter(yUSNO3, yOFFSET3,c='blue', marker='^',label='   y')
plt.scatter(xCCD3, xOFFSET3,c='red', marker='v')
plt.scatter(yCCD3, yOFFSET3,c='blue', marker='^')
#getting rid of the offsets that are large
#plt.axis((0,2045,-100,100))
plt.legend()
plt.savefig("d4_9.pdf",bbox_inches = 'tight')
plt.show()
#finding the big X and big Y using T1

#location from pictures
smallxy =  [1634.51,1103.37,1]
#finding big X and big Y standard coordinates
bigxy1 =  np.dot(T_inv,smallxy)
#transforming units from standard to radians
tan1 = -(bigxy1[0]/((np.cos(deltanot))-(bigxy1[1]*np.sin(deltanot))))
tan1 = np.arctan(tan1)
alpha1=tan1+ alphanot
delta1 = np.sqrt(((bigxy1[0])**2)+((bigxy1[1])**2)+1)
delta1 =( np.arcsin(((np.sin(deltanot))+(((bigxy1[1])*np.cos(deltanot))))/delta1))
print('alpha and delta are',alpha1,delta1)

"""
#finding the big X and big Y using T2
T2_inv=np.linalg.inv(T2)
bigxy2 =  np.dot(T2_inv,smallxy)

tan2=-(bigxy2[0]/((np.cos(deltanot))-(bigxy2[1]*np.sin(deltanot))))
tan2= np.arctan(tan2)
alpha2=tan2+ alphanot
delta2 =(np.arcsin((np.sin(deltanot)+((bigxy2[1])*np.cos(deltanot)))/(np.sqrt((1)+((bigxy2[1])**2)+((bigxy2[0])**2)))))
print('alpha and delta are',alpha2,delta2)
"""

bigxy3 =  np.dot(T3_inv,smallxy)
tan3=-(bigxy3[0]/((np.cos(deltanot))-(bigxy3[1]*np.sin(deltanot))))
tan3= np.arctan(tan3)
alpha3=tan3+ alphanot
delta3 =(np.arcsin((np.sin(deltanot)+((bigxy3[1])*np.cos(deltanot)))/(np.sqrt((1)+((bigxy3[1])**2)+((bigxy3[0])**2)))))
print('alpha and delta are',alpha3,delta3)

#finding the error
#since B and c are different shapes, a loop has to be formed to find the dot product
#finding Chi squared
#DOF
#DOF = (len(xOFFSET)) - 3
#T1 x
bdotc = []
matrixbarray = np.array(matrix_B)
for i in matrix_B:
    dotx1 = np.dot(i,cx)
    bdotc.append(dotx1)
bdotc = np.array(bdotc)
aminusbc = np.subtract(ax,bdotc)
aminusbctrans = np.transpose(aminusbc)
chi_T1_x = (np.dot(aminusbc,aminusbctrans))
#T1 y
bdotcy = []
matrixbarrayy = np.array(matrix_B)
for i in matrix_B:
    doty1 = np.dot(i,cy)
    bdotcy.append(doty1)
bdotcy = np.array(bdotcy)
aminusbcy = np.subtract(ay,bdotcy)
aminusbctransy = np.transpose(aminusbcy)
chi_T1_y = (np.dot(aminusbcy,aminusbctransy))
DOF2 = (len(xOFFSET3)) - 3
#T2 x
bdotct2 = []
matrixbarrayt2 = np.array(matrix_B3)
for i in matrix_B3:
    dotx1t2 = np.dot(i,cx3)
    bdotct2.append(dotx1t2)
bdotct2 = np.array(bdotct2)
aminusbct2 = np.subtract(ax,bdotct2)
aminusbctranst2 = np.transpose(aminusbct2)
chi_T2_x = (np.dot(aminusbct2,aminusbctranst2))
#T2 y
bdotcyt2 = []
matrixbarrayyt2 = np.array(matrix_B3)
for i in matrix_B3:
    doty1t2 = np.dot(i,cy3)
    bdotcyt2.append(doty1t2)
bdotcyt2 = np.array(bdotcyt2)
aminusbcyt2 = np.subtract(ay,bdotcyt2)
aminusbctransyt2 = np.transpose(aminusbcyt2)
chi_T2_y = (np.dot(aminusbcyt2,aminusbctransyt2))
#plotting the overlay to check between new CCD after T and USNO data
plt.figure(figsize=(9,4))
plt.plot(xUSNO,yUSNO,'o',label='USNO')
plt.plot(xCCD,yCCD,'x',c='red',label='CCD after T')
plt.legend()
plt.savefig("d4_10.pdf",bbox_inches = 'tight')
plt.show()

#RMS
meanxT1 = (np.sum(xOFFSET))/((len(xOFFSET)))
for i in xOFFSET:
    RMSxT1 = np.sqrt(((np.sum(((i-(meanxT1)))**2)))/len(xOFFSET))
meanyT1 = (np.sum(yOFFSET))/(len(yOFFSET))
for i in yOFFSET:
    RMSyT1 = np.sqrt(((np.sum(((i-(meanyT1)))**2)))/len(xOFFSET))
    
meanxT2 = (np.sum(xOFFSET3))/((len(xOFFSET3)))
for i in xOFFSET3:
    RMSxT2 = np.sqrt(((np.sum(((i-(meanxT2)))**2)))/len(xOFFSET3))
meanyT2 = (np.sum(yOFFSET3))/((len(yOFFSET3)))
for i in yOFFSET3:
    RMSyT2 = np.sqrt(((np.sum(((i-(meanyT2)))**2)))/len(xOFFSET3))
#error in arcsec
effectivescale = fp3
RMS_alpha = (RMSxT2/effectivescale)*206265
RMS_delta = (RMSyT2/effectivescale)*206265
error_arcsec = (np.sqrt((RMS_alpha**2)+(RMS_delta**2)))
#1.043
#error in degrees
RMS_xDeg = 0.000277778*RMS_alpha
RMS_yDeg = 0.000277778*RMS_delta