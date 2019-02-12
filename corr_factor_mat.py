#Vmax > 10 km/s center big galaxy
import numpy as np
import sys, os

if len(sys.argv) > 5 and len(sys.argv) < 3:
    print("Usage: python {0} <lstn> <Vpeak or Vmax> [M_t = 10] [host = 1]".format(sys.argv[0]))
    sys.exit()

#Info: [type](SDSS,DES,CLA<classic>), [Name], [absolute magnitude],[actual distance from the sun](kpc), 
#      [detect range of SDSS](kpc), [detect range of DES](kpc)
Satellite = [["SDSS","Bootes_I",-6.3,60,352,769],["SDSS","Bootes_II",-2.7,43,67,116],["SDSS","Canes_Venatici_I",-8.6,224,1017,2575],["SDSS",'Canes_Venatici_II',-4.9,151,185,369],['SDSS','Coma',-4.1,44,128,243],['SDSS','Hercules',-6.6,138,405,901],['SDSS,','Leo_IV',-5.0,158,194,389],['SDSS','Leo_T',-8.0,417,771,1879],['SDSS','Segue_1',-1.5,23,39,62],['SDSS','Ursa_Major_I',-5.5,106,244,506],['SDSS','Ursa_Major_II',-4.2,32,134,256],['SDSS','Willman_1',-2.7,38,67,116],['DES','Columba_1',-4.5,182,154,299],['DES','Grus_1',-3.4,120,93,168],['DES','Grus_2',-3.9,53,117,218],['DES','Horologium_1',-3.4,79,93,168],['DES','Horologium_2',-2.6,78,64,110],['DES','Indus_2',-4.3,214,140,269],['DES','Phoenix_2',-2.8,83,70,123],['DES','Pictoris_1',-3.1,115,81,143],['DES','Reticulum_2',-2.7,30,67,116],['DES','Reticulum_3',-3.3,92,89,159],['DES','Tucana_2',-4.3,58,140,269],['DES','Tucana_3',-2.4,25,58,99],['DES','Tucana_4',-3.5,48,97,177],['DES','Tucana_5',-1.6,55,40,65],['CLA','Carina',-9.4,94,1469,3919],['CLA','Draco',-9.4,79,1469,3919],['CLA','Fornax',-13.1,138,8074,27340],['CLA','LMC',-18.5,49,97076,465586],['CLA','Leo_U',-11.9,270,4646,14561],['CLA','Leo_II', -10.1,205,2028,5660],['CLA','Ursa_Minor',-8.9,69,1167,3014],['CLA','SMC',-17.1,63,50946,223254],['CLA','Sculptor',-9.8,88,1767,4835],['CLA','Sextans',-9.5,86,1539,4130],['CLA','Sagittarius',-15,28,19369,74131],['DES','Eridanus_II',-7.4,330,585,1371],['DES','Eridanus_III',-2,87,49,81],['SDSS','Leo_V',-5.2,178,212,432],['SDSS','Segue_II',-2.5,35,61,105],['MagLiteS','Pictor_II',-3.2,45,84,151],['SDSS','Pegasus_III',-4.1,205,127,242],['SDSS','Bootes_III',-5.8,46,280,592],['Gaia_ESO','Triagullum_II',-1.8,30,44,73],["Pan_Starrs",'Draco_II',-2.9,67,74,129],['CFHT','Munoz_I',-0.4,45,23,35],['SDSS','Pisces_II',-5,182,194,389],['DES','DES_J0225',-1.2,24,34,53],['DES','Hydra_II',-5.1,151,203,410],['HSSCSSP','Virgo_I',-0.8,87,30,43]]
#51
data_files = [i for i in os.listdir('../ELVIS/PairedTrees/') if "README" not in i]

if len(sys.argv) > 3:
    M_t = int(sys.argv[3])
else:
    M_t = 10

if len(sys.argv) > 4:
    MW = int(sys.argv[4])
    Andr = 1
    fname = "_mask_small"
else:
    MW = 1
    Andr = 0
    fname = "_mask_big"

pnts = 10 

lstn = int(sys.argv[1])
use_vmax = "max" in sys.argv[2].lower()

DES_cone_len = Satellite[lstn][5]/1000.
SDSS_cone_len = Satellite[lstn][4]/1000.  #t2 3

DES_cone_size = 0.7259
SDSS_coneS_size = 0.7259  #t2
SDSS_coneB_size = 0.9929   #t3    

corr_factors = np.zeros((len(data_files),pnts))
count_mas = np.zeros((len(data_files),pnts))

if use_vmax:
    fname = Satellite[lstn][1]+fname+"_CF_Vmax.txt"
else:
    fname = Satellite[lstn][1]+fname+"_CF_Vp{}.txt".format(str(int(M_t)))

def cart2sph(x,y,z):
    XsqPlusYsq = np.square(x) + np.square(y)
    r = np.sqrt(XsqPlusYsq + np.square(z))         # r
    elev = np.arctan2(z,np.sqrt(XsqPlusYsq))       # phi
    az = np.arctan2(y,x)                           # theta
    return r, elev, az   

def calc_open_angle(phi_x,the_x,phi_t,the_t):
    return np.arccos(np.cos(phi_t) * np.cos(phi_x) + np.sin(phi_t) * np.sin(phi_x) * np.cos(the_x - the_t)) 

def calc_corr_factors(Data_X, Data_Y, Data_Z, Data_M, Data_Rvir, pnts, MW, Andr):
    # Rvir converted in Mpc (was kpc) 
    RO2 = Data_Rvir[Andr]/1000.

    # in Mpc
    z0_cens = np.vstack([Data_X,Data_Y,Data_Z]).T
    pair_dist = z0_cens[Andr] - z0_cens[MW]   # First two are the hosts (M31 & MW)
    coord_dist = z0_cens[2:] - z0_cens[MW]

    dis_g = np.linalg.norm(pair_dist)

    angle_BG = np.arctan(RO2/dis_g)
    
    _, phi_x, the_x = cart2sph(coord_dist[:,0],coord_dist[:,1],coord_dist[:,2])
    _, phi_t, the_t = cart2sph(-pair_dist[0],-pair_dist[1],-pair_dist[2])

    ang_dis = calc_open_angle(phi_x,the_x,phi_t,the_t)

    Andr_gal = ang_dis <= angle_BG # Mask things that are in the same area as Andr.
    Andr_gal = Andr_gal[np.newaxis].T & np.ones(pnts,dtype = 'bool') # Make into N x pnts matrix

    #random  
    #np.random.seed(777)
    u = np.random.uniform(0,1,pnts)
    v = np.random.uniform(0,1,pnts)

    the_t = 2. * np.pi * u
    phi_t = np.arccos(2. * v - 1)

    the_t2 = the_t - 0.09
    phi_t2 = phi_t + 1.18
    the_t3 = the_t + 3.1
    phi_t3 = phi_t + 1.72    

    the_t2[ the_t2 < 0 ] += 2. * np.pi
    phi_t2[ phi_t2 > np.pi ] -= np.pi
    the_t3[ the_t3 > 2*np.pi ] -= 2*np.pi
    phi_t3[ phi_t3 > np.pi ] -= np.pi

    s = cart2sph(coord_dist[:,0],coord_dist[:,1],coord_dist[:,2])

    phi_x = (s[1])[np.newaxis].T
    the_x = (s[2])[np.newaxis].T
    r = s[0]

    dis = calc_open_angle(phi_x,the_x,phi_t,the_t)
    dis2 = calc_open_angle(phi_x,the_x,phi_t2,the_t2)
    dis3 = calc_open_angle(phi_x,the_x,phi_t3,the_t3)

    r_msk = (r[np.newaxis].T < np.ones(pnts)) & (~Andr_gal)

    r_SDSS = (r[np.newaxis].T < np.ones(pnts) * SDSS_cone_len) & r_msk
    s_SDSS = (dis2 <= SDSS_coneS_size) & r_SDSS
    b_SDSS = (dis3 <= SDSS_coneB_size) & r_SDSS

    in_SDSS = s_SDSS | b_SDSS
    in_DES = (r[np.newaxis].T < np.ones(pnts) * DES_cone_len) & (dis <= DES_cone_size) & r_msk

    in_cones = in_SDSS | in_DES

    count_mas = np.sum(in_cones[Data_M[2:] >= M_t], axis = 0 ) * 1.
    count_tot = np.sum(r_msk[Data_M[2:] >= M_t], axis = 0 ) * 1.

    return count_tot,count_mas


for file_name in xrange(len(data_files)):
    # file.write(" NEW " + file_name+' ')
    if use_vmax:
        Data_M = np.loadtxt('../ELVIS/PairedTrees/'+data_files[file_name]+'/Vmax.txt')[:,0]
    else:
        Data_M = np.amax(np.loadtxt('../ELVIS/PairedTrees/'+data_files[file_name]+'/Vmax.txt'),axis=1)
    Data_X = np.loadtxt('../ELVIS/PairedTrees/'+data_files[file_name]+'/X.txt')[:,0]
    Data_Y = np.loadtxt('../ELVIS/PairedTrees/'+data_files[file_name]+'/Y.txt')[:,0]
    Data_Z = np.loadtxt('../ELVIS/PairedTrees/'+data_files[file_name]+'/Z.txt')[:,0]
    Data_Rvir = np.loadtxt('../ELVIS/PairedTrees/'+data_files[file_name]+'/Rvir.txt')[:,0]

    corr_factors[file_name],count_mas[file_name] = calc_corr_factors(Data_X, Data_Y, Data_Z, Data_M, Data_Rvir, pnts, MW, Andr)

np.savetxt(fname,corr_factors.T)
np.savetxt("in_cone.txt",count_mas.T)  






