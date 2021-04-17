import urllib.request as url
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

fo = fits.open('NGC7331-S001-R001-C001-r.fts')

fo[0].header

f_d=fo[0].data

hdu = fo[0]
print(type(f_d))
print(f_d.shape)

#plt.figure(11,figsize=(15,15)) #shape of CCD
plt.imshow(f_d,cmap='gray',origin = "lower",vmin=8000,vmax=12000) #plotting with high brightness to viw image
plt.title('CCD image')
plt.colorbar()
plt.show()
plt.imshow(f_d,cmap='gray_r',origin = "lower",vmin=8000,vmax=12000) #plotting with high brightness to viw image
plt.title('CCD image')
plt.colorbar()
plt.show()
dark = fits.open('Dark-S001-R003-C003-B2.fts')
dark[0].header
dark_d = dark[0].data
flat = fits.open('combflatr.fits')
flat[0].header
flat_d = flat[0].data
numert = np.abs(np.subtract(f_d,dark_d))
newdata = np.divide(numert,flat_d)

listpeak =[]
listpeaks =[]
import skimage
from skimage.feature import peak_local_max
f_d = f_d[::-1]
peaks = peak_local_max(f_d,indices=True,num_peaks = 500,min_distance=100)
a = np.transpose(peaks)
plt.scatter(a[1],a[0],marker='+')
row = a[1]
col = a[0]
plt.title("NGC7331")
plt.xlabel("Pixel")
plt.ylabel("Pixel")
plt.tight_layout()




centroids_xpos = []
centroids_ypos = []
for i in range(500):
    (x, y) = (row[i], col[i])
    numeratorx = 0
    numeratory = 0
    total_intensity = 0
    for j in range(-3, 4): 
        for k in range(-3, 4):
            (xi, yi) = (x+k, y+j)
            if (xi > 2045) or (yi > 2045):
                pass
            else:
                curr_intensity = f_d[xi][yi]
                total_intensity += curr_intensity
                numeratorx += xi*curr_intensity
                numeratory += yi*curr_intensity
    xpos = numeratorx / total_intensity
    ypos = numeratory / total_intensity
    centroids_xpos.append(xpos)
    centroids_ypos.append(ypos)
    
plt.scatter(centroids_xpos, centroids_ypos, alpha=0.1,c='orange')
plt.show()


        
# =============================================================================
#         average = np.sum(np.array([i+col[j]]))*(np.array([row[j]]))*np.sum(np.array([col[j]+i]))
#         sumofintense = np.array([i+col[j]])
#         x = np.divide(average,sumofintense)
# =============================================================================

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
    ras = fo[0].header['RA']
    des = fo[0].header['DEC']
    radeg = 15*(float(ras[0:2])+float(ras[3:5])/60.+float(ras[6:])/3600)
    dsgn = np.sign(float(des[0:3]))
    dedeg = float(des[0:3])+dsgn*float(des[4:6])/60.+dsgn*float(des[7:])/3600.
    fovam = 37
    name,rad,ded,rmag = usno(radeg,dedeg,fovam)
    print(rmag) #there's missing mags, but you don't need mags anyways.
    mask = np.where(rmag<13)[0]
    rad = rad[mask]
    ded = ded[mask]
    
    #plot positions
    plt.title("USNO-B1 Stars")
    plt.scatter(rad, ded, marker='x',c='red')
    plt.xlabel("RA [deg]")
    plt.ylabel("DEC [deg]")
    plt.tight_layout()

    plt.show()

alist = []
rad_radians = np.radians(rad)
ded_radians = np.radians(ded)
alphanot = np.radians(radeg)
deltanot = np.radians(dedeg)
p = 0.018
X_values = []
for i,j in zip(rad_radians,ded_radians):
        numX_dividedby_denX = (np.divide((np.cos(j)*np.sin(i-alphanot)),(np.cos(deltanot)*np.cos(j)*np.cos(i-alphanot) + np.sin(j)*np.sin(deltanot))))
        X_values.append(numX_dividedby_denX)
       
X_values_array = np.asarray(X_values)
x_pixels = 3454*(X_values_array/p)+1024

Y_values = []
for i,j in zip(rad_radians,ded_radians):
        numY_dividedby_denY = -(np.divide((np.sin(deltanot)*np.cos(j)*np.cos(i-alphanot)-np.cos(deltanot)*np.sin(j)),(np.cos(deltanot)*np.cos(j)*np.cos(i-alphanot) + np.sin(j)*np.sin(deltanot))))
        Y_values.append(numY_dividedby_denY)
       
Y_values_array = np.asarray(Y_values)
y_pixels = 3454*(Y_values_array/p)+1024

plt.title("Pixel Positions of USNO-B1 Stars")
plt.scatter(x_pixels, y_pixels, marker='x',c='red')
plt.show()

#overlay
plt.scatter(centroids_xpos, centroids_ypos, marker='+')
plt.scatter(x_pixels, y_pixels, marker='x',c='red')
plt.show()


x_star=[]
y_star=[]
xusno = []
yusno= []
xoffset =[]
yoffset =[]
length_rad = len(rad)
for i in range(500):
    raccd = centroids_xpos[i]
    deccd = centroids_ypos[i]
    for j in range(len(rad)):
        rausno = x_pixels[j]
        decusno = y_pixels[j]
        separation = np.sqrt(((raccd-rausno)**2)+((deccd-decusno)**2)) 
        if separation < 10:
            x_star.append(raccd)
            y_star.append(deccd)
            xusno.append(rausno)
            yusno.append(decusno)
            xoffset.append(raccd-rausno)
            yoffset.append(deccd-decusno)



plt.plot(xusno,yusno,'o',c='red', marker='x',label='USNO')
plt.plot(x_star,y_star,'o',c='blue',marker='+',label='CCD')
plt.legend(loc='best')
plt.show()

plt.title("Residuals")
plt.xlabel("x or y [pixel]")
plt.ylabel("Pixel offset")
plt.scatter(xusno, xoffset,c='red', marker='v',label='   x')
plt.scatter(yusno, yoffset,c='blue', marker='^',label='   y')

plt.legend(loc='best')
plt.show()


fp = 3454.00/0.018      # = 190080
ax = np.zeros(len(X_values_array))
for i in range(len(x_star)):
    ax[i] = x_star[i]
ones = np.ones([len(X_values_array)])
matrix_Bt = np.array([190080*X_values_array,190080*Y_values_array,ones])

matrix_Bt_ax = np.dot(matrix_Bt,ax)
matrix_B = np.transpose(matrix_Bt)
Bt_times_B = np.dot(matrix_Bt,matrix_B)
Bt_B_inverse = np.linalg.inv(Bt_times_B)
cx = np.dot(Bt_B_inverse,matrix_Bt_ax)


ay = np.zeros(len(Y_values_array))
for i in range(len(y_star)):
    ay[i] = y_star[i]
ones = np.ones([len(Y_values_array)])
matrix_Bt = np.array([190080*X_values_array,190080*Y_values_array,ones])
#matrix_Bt = np.transpose(matrix_B)
matrix_Bt_ay = np.dot(matrix_Bt,ay)
matrix_B = np.transpose(matrix_Bt)
Bt_times_B = np.dot(matrix_Bt,matrix_B)
Bt_B_inverse = np.linalg.inv(Bt_times_B)
cy = np.dot(Bt_B_inverse,matrix_Bt_ay)
print('Plate Constants')
print('a_11, a_12, x_0 = ',cx)
print('a_21, a_22, y_0 = ',cy)

#T matrix
T = [[fp*cx[0] , fp*cx[1] , cx[2]], [fp*cy[0] , fp*cy[1] , cy[2]], [0, 0, 1]]
print('T Matrix = ',T)
