#Vmax > 10 km/s center big galaxy
lstn = 0

import matplotlib
matplotlib.use('Agg')

import numpy as np
import math as m
import matplotlib.pyplot as plt
import scipy.interpolate
import os
from mpl_toolkits.mplot3d import Axes3D

#Info: [type](SDSS,DES,CLA<classic>), [Name], [absolute magnitude],[actual distance from the sun](kpc), 
#      [detect range of SDSS](kpc), [detect range of DES](kpc)
Satellite = [["SDSS","Bootes_I",-6.3,60,352,769],["SDSS","Bootes_II",-2.7,43,67,116],["SDSS","Canes_Venatici_I",-8.6,224,1017,2575],["SDSS",'Canes_Venatici_II',-4.9,151,185,369],['SDSS','Coma',-4.1,44,128,243],['SDSS','Hercules',-6.6,138,405,901],['SDSS,','Leo_IV',-5.0,158,194,389],['SDSS','Leo_T',-8.0,417,771,1879],['SDSS','Segue_1',-1.5,23,39,62],['SDSS','Ursa_Major_I',-5.5,106,244,506],['SDSS','Ursa_Major_II',-4.2,32,134,256],['SDSS','Willman_1',-2.7,38,67,116],['DES','Columba_1',-4.5,182,154,299],['DES','Grus_1',-3.4,120,93,168],['DES','Grus_2',-3.9,53,117,218],['DES','Horologium_1',-3.4,79,93,168],['DES','Horologium_2',-2.6,78,64,110],['DES','Indus_2',-4.3,214,140,269],['DES','Phoenix_2',-2.8,83,70,123],['DES','Pictoris_1',-3.1,115,81,143],['DES','Reticulum_2',-2.7,30,67,116],['DES','Reticulum_3',-3.3,92,89,159],['DES','Tucana_2',-4.3,58,140,269],['DES','Tucana_3',-2.4,25,58,99],['DES','Tucana_4',-3.5,48,97,177],['DES','Tucana_5',-1.6,55,40,65],['CLA','Carina',-9.4,94,1469,3919],['CLA','Draco',-9.4,79,1469,3919],['CLA','Fornax',-13.1,138,8074,27340],['CLA','LMC',-18.5,49,97076,465586],['CLA','Leo_U',-11.9,270,4646,14561],['CLA','Leo_II', -10.1,205,2028,5660],['CLA','Ursa_Minor',-8.9,69,1167,3014],['CLA','SMC',-17.1,63,50946,223254],['CLA','Sculptor',-9.8,88,1767,4835],['CLA','Sextans',-9.5,86,1539,4130],['CLA','Sagittarius',-15,28,19369,74131],['DES','Eridanus_II',-7.4,330,585,1371],['DES','Eridanus_III',-2,87,49,81],['SDSS','Leo_V',-5.2,178,212,432],['SDSS','Pisces_II',-5,182,194,389],['SDSS','Segue_II',-2.5,35,61,105],['MagLiteS','Pictor_II',-3.2,45,84,151],['SDSS','Pegasus_III',-4.1,205,127,242],['SDSS','Bootes_III',-5.8,46,280,592],['Gaia_ESO','Triagullum_II',-1.8,30,44,73],["Pan_Starrs",'Draco_II',-2.9,67,74,129],['CFHT','Munoz_I',-0.4,45,23,35],['SDSS','Pisces_II',-5,182,194,389],['DES','DES_J0225',-1.2,24,34,53],['DES','Hydra_II',-5.1,151,203,410],['HSSCSSP','Virgo_I',-0.8,87,30,43]]
#51

M_t = 10
DES_cone_len = (Satellite[lstn][5])/1000.
SDSS_cone_len = (Satellite[lstn][4])/1000.  #t2 3

DES_cone_size = 0.7259
SDSS_coneS_size = 0.7259  #t2
SDSS_coneB_size = 0.9929   #t3

def cart2sph(x,y,z):
    XsqPlusYsq = x**2 + y**2
    r = m.sqrt(XsqPlusYsq + z**2)               # r
    elev = m.atan2(z,m.sqrt(XsqPlusYsq))     # phi
    az = m.atan2(y,x)                           # theta
    return r, elev, az

def spherical_to_cartesian(theta,phi):
    d = np.cos(theta) * np.sin(phi)
    g = np.sin(phi) * np.sin(theta)
    h =np.cos(phi)
    return d,g,h       

file = open(Satellite[lstn][1]+"_CF_Vmax.txt", "w")
data_files = [i for i in os.listdir('../ELVIS/PairedTrees/') if "README" not in i]

fig = plt.figure()
for file_name in data_files:
    file.write(" NEW " + file_name+' ')
    Data_M = np.loadtxt('../ELVIS/PairedTrees/'+file_name+'/Vmax.txt')
    Data_X = np.loadtxt('../ELVIS/PairedTrees/'+file_name+'/X.txt')
    Data_Y = np.loadtxt('../ELVIS/PairedTrees/'+file_name+'/Y.txt')
    Data_Z = np.loadtxt('../ELVIS/PairedTrees/'+file_name+'/Z.txt')
    Data_Rvir = np.loadtxt('../ELVIS/PairedTrees/'+file_name+'/Rvir.txt')
    
    Plot_X = []
    Plot_X_1 = []

    #Big_gal: list needed to be blocked for big galaxy, aka the small glaxay
    Big_gal = []
    Big_gal_the = []
    Big_gal_phi = []
    Small_gal = []
    Small_gal_the = []
    Small_gal_phi = []

    # Rvir converted in Mpc (was kpc) 
    RO1 = (Data_Rvir[0][0])/1000.
    RO2 = (Data_Rvir[1][0])/1000.

    # in Mpc
    XO1 = (Data_X[0][0])
    YO1 = (Data_Y[0][0])
    ZO1 = (Data_Z[0][0])
    XO2 = (Data_X[1][0])
    YO2 = (Data_Y[1][0])
    ZO2 = (Data_Z[1][0])

    dis_g1 = cart2sph(XO1-XO2,YO1-YO2,ZO1-ZO2)
    dis_g = dis_g1[0]
    #angle of block for big galaxy
    angle_BG = np.arctan(RO2/dis_g)
    angle_SG = np.arctan(RO1/dis_g)
    B_gal = cart2sph(XO2-XO1,YO2-YO1,ZO2-ZO1)
    S_gal = cart2sph(XO1-XO2,YO1-YO2,ZO1-ZO2)
    
    for i in range(2,len(Data_M)):
        x1 = (Data_X[i][0]) - XO1
        y1 = (Data_Y[i][0]) - YO1
        z1 = (Data_Z[i][0]) - ZO1
        x2 = (Data_X[i][0]) - XO2
        y2 = (Data_Y[i][0]) - YO2
        z2 = (Data_Z[i][0]) - ZO2
    
        #consider for big gal
        s = cart2sph(x1,y1,z1)
        phi_x = s[1] #+ 0.5 * np.pi
        the_x = s[2] #+ np.pi
        a = cart2sph(XO2-XO1,YO2-YO1,ZO2-ZO1)
        phi_t = a[1]
        the_t = a[2]
        dis = np.arccos(np.cos(phi_t)*np.cos(phi_x)+np.sin(phi_t)*np.sin(phi_x)*np.cos(the_x-the_t))
        r = a[0]
        if dis <= angle_BG:
            Big_gal.append(i)
            Big_gal_the.append(the_x)
            Big_gal_phi.append(phi_x)


   
    #random
    pj_the = []
    pj_phi = []  
    np.random.seed(777) 
    u = np.random.uniform(0,1,10)
    v = np.random.uniform(0,1,10)
    for i in xrange(len(u)):
        u1 = u[i]
        v1 = v[i]
        the_t = np.pi*2*u1
        phi_t = np.arccos(2*v1-1)
    
        two_the = the_t - 0.09
        two_phi = phi_t + 1.18
    
        three_the = the_t + 3.1
        three_phi = phi_t + 1.72
    
        if two_the < 0:
            the_t2 = 6.2832 + two_the
        if two_the >= 0:
            the_t2 = two_the
        if two_phi <= np.pi:
            phi_t2 = two_phi
        if two_phi >np.pi:
            phi_t2 = two_phi - np.pi
    
        if three_the <= 2*np.pi:
            the_t3 = three_the
        if three_the > 2*np.pi:
            the_t3 = three_the - 2*np.pi
        if three_phi <= np.pi:
            phi_t3 = three_phi
        if three_phi > np.pi:
            phi_t3 = three_phi - np.pi
    
        #already got t 123
        pj_the.append(the_t)
        pj_phi.append(phi_t)
        count_tot = []
        count_mas = []
        #determin in cone
        f_theta = []
        f_phi = []
        to_phi = []
        t_theta = []
        for n in range(2,len(Data_X)):
            x1 = Data_X[n][0] - XO1
            y1 = Data_Y[n][0] - YO1
            z1 = Data_Z[n][0] - ZO1
            x2 = Data_X[n][0] - XO2
            y2 = Data_Y[n][0] - YO2
            z2 = Data_Z[n][0] - ZO2
            s = cart2sph(x1,y1,z1)
            phi_x = s[1]+0.5*np.pi
            the_x = s[2]+np.pi
            r = s[0]
            dis = np.arccos(np.cos(phi_t)*np.cos(phi_x)+np.sin(phi_t)*np.sin(phi_x)*np.cos(the_x-the_t))
            dis2 = np.arccos(np.cos(phi_t2)*np.cos(phi_x)+np.sin(phi_t2)*np.sin(phi_x)*np.cos(the_x-the_t2))
            dis3 = np.arccos(np.cos(phi_t3)*np.cos(phi_x)+np.sin(phi_t3)*np.sin(phi_x)*np.cos(the_x-the_t3))
            ma = Data_M[n][0]
            if n not in Big_gal and r <= 1:
            #if r <= 0.05:
                if r < SDSS_cone_len:
                    
                    if dis2 <= SDSS_coneS_size:
                        f_theta.append(the_x)
                        f_phi.append(phi_x)
                        if ma >= M_t and n not in count_mas:
                            count_mas.append(n)
                            if n not in count_tot:
                                count_tot.append(n)
                    elif dis3 <= SDSS_coneB_size:
                        f_theta.append(the_x)
                        f_phi.append(phi_x)
                        if ma >= M_t and n not in count_mas:
                            count_mas.append(n)
                            if n not in count_tot:
                                count_tot.append(n)
                    elif ma >= M_t and n not in count_tot:
                        if dis2 > SDSS_coneS_size and dis3 > SDSS_coneB_size:
                            t_theta.append(the_x)
                            to_phi.append(phi_x)
                            count_tot.append(n)
                            
                                


                
                if r < DES_cone_len:
                    if dis <= DES_cone_size:
                        f_theta.append(the_x)
                        f_phi.append(phi_x)
                        if ma >= M_t and n not in count_mas:
                            count_mas.append(n)
                            if n not in count_tot:
                                count_tot.append(n)
                    elif ma >= M_t and n not in count_tot:
                        t_theta.append(the_x)
                        to_phi.append(phi_x)
                        count_tot.append(n)
                        
                            

                if r >= SDSS_cone_len and r >= DES_cone_len:
                    if ma >= M_t and n not in count_tot:
                        count_tot.append(n)
        #res = (len(count_mas))/(len(count_tot))
        #Plot_X.append(res)
    
        if len(count_mas) != 0:
            # print "Pointing {0}:".format(i)
            # print "\t",len(count_tot),"  ",len(count_mas)
            res1 =  (1.*len(count_tot))/(len(count_mas))   
            Plot_X_1.append(res1)
        else:
            res1 = 0.0
        file.write(str(res1)+' ')
    
 
    #bar = np.linspace(0.0, 10, num=51)
    #plt.hist(Plot_X_1, bins=bar,histtype='step') 
           
file.close()    






