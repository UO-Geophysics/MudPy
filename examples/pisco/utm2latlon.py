'''
Convert insar data from UTM to lat/lon
'''

from pyproj import Proj
from numpy import genfromtxt,c_,savetxt
        
#Make projection object
p = Proj(proj='utm',zone='18K',ellps='WGS84')
#Read data
los=genfromtxt('/Users/dmelgar/Pisco2007/InSAR/DataFromSAR.csv',delimiter=',')
#Convert (inverse=True goes from UTM to lat/lon)
LL=p(los[:,0],los[:,1],inverse=True)
#Write to file
out=c_[LL[0],LL[1],los[:,2]]
savetxt('/Users/dmelgar/Pisco2007/InSAR/DataFromSAR.csv',out,fmt='%12.6f\t%12.6f\t%8.4f