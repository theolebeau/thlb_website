import numpy as np
#import matplotlib
#matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
from scipy.io import FortranFile
from astropy import constants as const
import pickle
import functions as f
# from sklearn.utils import resample
import time
import sys
from scipy.optimize import curve_fit
import f_to_py as ftp

kb = const.k_B.to_value('J/K')
mp = const.m_p.to_value('kg')
mu = 0.60
g = const.G.to_value('m3/(kg s2)')
H = 67.7699966  # hubble constant in km/s/mpc
omega_m = 0.307114988565445E+00
omega_l = 0.692884981632233E+00
omega_k = 0.298023223876953E-07
omega_b = 0.450000017881393E-01
omega_c = omega_m - omega_b
c = const.c.to_value('m/s')
pc = const.pc.to_value('m')
rho_c = ((3 * (H * (1E3 / (pc * 1E6))) ** 2) / (8 * np.pi * g))
# *(omega_c/omega_m)
m_sun = const.M_sun.to_value('kg')
kev = 1.602176634e-16

unit_l = 0.227550518265197E+28
unit_d = 0.264088847033921E-29
unit_t = 0.455724529603677E+18

start=time.time()

#rbin = np.array([i * 10 + 5 for i in range(400)])
rlogbin = np.array([10 ** ((35 + i) * 0.05) for i in range(40)])
rlogbin_cen = np.array([10 ** ((35.5 + i) * 0.05) for i in range(40)])

#print("rlogbin",rlogbin)
#print("rlogbin_cen",rlogbin_cen)

#print("rlogbin",rlogbin)

#f.velocity_structure_function('./maps/high_res/map_high_19_MW_v_los.bin')
#f.ratios_dep()

#f.ellipticity()

#f.relaxedness()
#sys.exit()
#ratio=(proj_P_21-proj_P_19)/proj_P_19
#mean_ratio=np.mean(np.absolute(ratio))
#print("ratio",ratio)
#print("mean ratio",mean_ratio)

#sys.exit()


#f.gianfagna_overplot()
#sys.exit()
#f.toy_model()
#f.intmass_bias_relation()
#f.los_distrib("virgo_xyz_hydro_l19_MW_los.dat")
#f.los_distrib("virgo_xyz_hydro_l19_fil_los.dat")
#f.ratios_dep()
#f.gian_bias_overplot()
#f.offset_virgo()


####splashback radius project

#f.radial_gal_vel()
#f.estimate_r200()
#f.r_splb_DM()
#f.r_splb_all()
#f.r_sp_proj("y",1)
#f.r_sp_proj_stack("sd",1)
#f.r_splb_gas()

#f.vrad_r_phase_space()

#### turbulence project

#f.test_PK()
#f.vel_power_spectrum()
#f.read_fit_file()
#f.velocity_operators_2D()
#f.rsp_summary_column()
#f.decompose_vel_field_2D(4,"-")
#f.vel_pk_nD()
#f.decompose_vel_field()
#f.show_decomposed_vel_fields()
#f.test_interpolation()
#f.energy_spectrum_2D()
#f.energy_spectrum_2D_two_col("tot")
#f.energy_spectrum_2D_one_row()
#f.energy_spectrum_2D_vs_3D_core_test()
#f.convolution_test()
#f.energy_spectrum_2D_vs_3D_test_all_slices()
#f.check_gal_pos_in_slice()

##project 4

#f.VSF(1)
#f.projection('./maps/high_res/velocity/15f15_analysis/map_high_15f15_cen_map_P_mw_Tsup_0.2kev_15Mpc2.bin',"P","P_rad_log_l15_20b.npy","std_P_rad_log_l15_20b.npy",0,0,1)
f.projection('./maps/high_res/ne_maps/map_high_19_z_ne_los.bin',"ne","./vel_Virgo_core_project/deproj_prof/ne_deproj_prof_z.npy","./vel_Virgo_core_project/deproj_prof/ne_deproj_prof_z_std.npy",0,1,0)
#f.projection('./maps/high_res/velocity/15f15_analysis/map_high_15f15_cen_map_P_mw_Tsup_8e-3kev_15Mpc2_gc_1e8.5.bin',"P","P_rad_log_l15_20b.npy","std_P_rad_log_l15_20b.npy",0,0,1)

sys.exit()

#f.PDF_3D("/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l15_high.dat")
#f.PDF_2D()
#print('test')
#f.masses_estimation(3)
#f.vel_rad_prof()
#f.estimate_overdensity()

##big sim

#f.create_namelist()

#f.scaling_relation_2D_vs_3D_stats()

#f.Mach_and_Pnt(2)

#f.scaling_relation_2D_vs_3D_stats_random_projs()
#f.thermo_profs_univ_prof_overplot()

#sys.exit()

def load_lograd_dat(file):
    rad_data = np.load(file)
    Plogmed = rad_data[0, :]
    Pstdlogmed = rad_data[1, :]
    Plogmean = rad_data[2, :]
    Pstdlogmean = rad_data[3, :]
    Tlogmed = rad_data[4, :]
    Tstdlogmed = rad_data[5, :]
    Tlogmean = rad_data[6, :]
    Tstdlogmean = rad_data[7, :]
    nelogmed = rad_data[8, :]
    nestdlogmed = rad_data[9, :]
    nelogmean = rad_data[10, :]
    nestdlogmean = rad_data[11, :]
    nlog = rad_data[12, :]
    mlog = rad_data[13, :]
    ndm = rad_data[14, :]
    mdm = rad_data[15, :]
    inf = rad_data[16, :]
    middle = rad_data[17, :]
    sup = rad_data[18, :]
    print("Loading ",file,"complete")
    return Plogmed,Pstdlogmed,Plogmean,Pstdlogmean,Tlogmed,Tstdlogmed,Tlogmean,Tstdlogmean,nelogmed,nestdlogmed,nelogmean,nestdlogmean,nlog,mlog,ndm,mdm,inf,middle,sup
load_data=0
if load_data==1:
    Plogmed_c8,Pstdlogmed_c8,Plogmean_c8,Pstdlogmean_c8,Tlogmed_c8,Tstdlogmed_c8,Tlogmean_c8,Tstdlogmean_c8,nelogmed_c8,nestdlogmed_c8,nelogmean_c8,nestdlogmean_c8,nlog_c8,mlog_c8,ndm_c8,mdm_c8,inf_c8,middle_c8,sup_c8=load_lograd_dat("lograd_data_lvl21_cleanm1e8.5.npy")
    Plogmed_17,Pstdlogmed_17,Plogmean_17,Pstdlogmean_17,Tlogmed_17,Tstdlogmed_17,Tlogmean_17,Tstdlogmean_17,nelogmed_17,nestdlogmed_17,nelogmean_17,nestdlogmean_17,nlog_17,mlog_17,ndm_17,mdm_17,inf_17,middle_17,sup_17=load_lograd_dat("lograd_data_lvl17_cleanm1e9.npy")
    #Plogmed_17o,Pstdlogmed_17o,Plogmean_17o,Pstdlogmean_17o,Tlogmed_17o,Tstdlogmed_17o,Tlogmean_17o,Tstdlogmean_17o,nelogmed_17o,nestdlogmed_17o,nelogmean_17o,nestdlogmean_17o,nlog_17o,mlog_17o,ndm_17o,mdm_17o,inf_17o,middle_17o,sup_17o=load_lograd_dat("lograd_data_lvl17_cleanm1e9_KPam-3.npy")

    #Plogmed_c8l,Pstdlogmed_c8l,Plogmean_c8l,Pstdlogmean_c8l,Tlogmed_c8l,Tstdlogmed_c8l,Tlogmean_c8l,Tstdlogmean_c8l,nelogmed_c8l,nestdlogmed_c8l,nelogmean_c8l,nestdlogmean_c8l,nlog_c8l,mlog_c8l,ndm_c8l,mdm_c8l,inf_c8l,middle_c8l,sup_c8l=load_lograd_dat("lograd_data_lvl21_cleanm1e8.5_left.npy")
    #Plogmed_c8r,Pstdlogmed_c8r,Plogmean_c8r,Pstdlogmean_c8r,Tlogmed_c8r,Tstdlogmed_c8r,Tlogmean_c8r,Tstdlogmean_c8r,nelogmed_c8r,nestdlogmed_c8r,nelogmean_c8r,nestdlogmean_c8r,nlog_c8r,mlog_c8r,ndm_c8r,mdm_c8r,inf_c8r,middle_c8r,sup_c8r=load_lograd_dat("lograd_data_lvl21_cleanm1e8.5_right.npy")

    Plogmed_c8, Pstdlogmed_c8, Plogmean_c8, Pstdlogmean_c8, Tlogmed_c8, Tstdlogmed_c8, Tlogmean_c8, Tstdlogmean_c8, nelogmed_c8, nestdlogmed_c8, nelogmean_c8, nestdlogmean_c8, nlog_c8, mlog_c8, ndm_c8, mdm_c8, inf_c8, middle_c8, sup_c8 = load_lograd_dat("lograd_data_lvl21_cleanm1e8.5.npy")

    #print("plogmedc8",Plogmed_c8)
    #print("plogmedc8 test",Plogmed_c8_t)
    #sys.exit()
    #print(nlog_c8)
    #sys.exit()

    ne_urban=np.loadtxt('./obs_data/ne_prof+errb_urban.txt')
    P_urban=np.loadtxt('./obs_data/P_prof+errb_urban.txt')
    ne_planck=np.loadtxt('./obs_data/ne_prof+errb_planck.txt')
    P_planck=np.loadtxt('./obs_data/P_prof+errb_planck.txt')
    ne_ghiz=np.loadtxt('./obs_data/ne_prof_ghizzardi.txt')
    P_planck_deg=np.loadtxt('./obs_data/P_prof_deg_planck.txt')
    #print("deg",P_planck_deg[:,0])
    r_deg21=np.tan(P_planck_deg[:,0]*(np.pi/180))*21e3
    r_deg16=np.tan(P_planck_deg[:,0]*(np.pi/180))*16e3
    #print("rdeg",r_deg)

    m_sum = np.array([np.sum(mdm_c8[0:i + 1]) + np.sum(mlog_c8[0:i + 1]) for i in range(40)])

    #Plogmed_c8-=1.136842

#plt.scatter(rlogbin_cen, Plogmean_17, s=6, c='blue')
#plt.errorbar(rlogbin_cen, Plogmean_17, yerr=Pstdlogmean_17, ls='dotted',label='3D, 2048^3, lvl 17 (max cell length=5,625kpc)', alpha=0.7, c='blue')
def threeD_reso_test():
    p14=np.loadtxt("3D_test_P_lvl14.txt")
    p14_errb=np.loadtxt("3D_test_P_lvl14_errb.txt")
    plt.scatter(rlogbin_cen, p14, s=6, c='blue')
    plt.errorbar(rlogbin_cen, p14, yerr=p14_errb, ls='dotted',label='3D, 2048^3, lvl 14 (max cell length=45kpc)', alpha=0.7, c='orange')
    p15=np.loadtxt("3D_test_P_lvl15.txt")
    p15_errb=np.loadtxt("3D_test_P_lvl15_errb.txt")
    plt.scatter(rlogbin_cen, p15, s=6, c='blue')
    plt.errorbar(rlogbin_cen, p15, yerr=p15_errb, ls='dotted',label='3D, 2048^3, lvl 15 (max cell length=22,5kpc)', alpha=0.7, c='green')
    p16=np.loadtxt("3D_test_P_lvl16.txt")
    p16_errb=np.loadtxt("3D_test_P_lvl16_errb.txt")
    plt.scatter(rlogbin_cen, p16, s=6, c='blue')
    plt.errorbar(rlogbin_cen, p16, yerr=p16_errb, ls='dotted',label='3D, 2048^3, lvl 16 (max cell length=11,25kpc)', alpha=0.7, c='green')

#plt.scatter(rlogbin_cen, Plogmean_17, s=6, c='blue')
#plt.errorbar(rlogbin_cen, Plogmean_17, yerr=Pstdlogmean_17, ls='dotted',label='3D, 2048^3', alpha=0.7, c='blue') #'Smallest grid size: 5.62kpc'
#Plogmean, Tlogmean, nelogmean, nlog, mlog, nestdlogmean, Pstdlogmean, Tstdlogmean = f.hydro("/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l15_rescale_low_nline.dat", 1, 1)
#sys.exit()

#f.dm("virgo_xyz_dm_high_res.dat")

sort_data = 0
if sort_data == 1:
    #print("test")
    #Plogmean, Tlogmean, nelogmean, nlog, mlog, nestdlogmean, Pstdlogmean, Tstdlogmean = f.hydro("/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l15_low.dat", 1, 1)
    #Plogmean, Tlogmean, nelogmean, nlog, mlog, nestdlogmean, Pstdlogmean, Tstdlogmean = f.hydro("/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21_gal_clean_m1e8.5.dat", 1, 2, 0, 77)
    Plogmean, Tlogmean, nelogmean, nlog, mlog, nestdlogmean, Pstdlogmean, Tstdlogmean = f.hydro("/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l21.dat", 1, 1, 0, 40)
    Plogmean2, Tlogmean2, nelogmean2, nlog2, mlog2, nestdlogmean2, Pstdlogmean2, Tstdlogmean2= f.hydro("/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l21.dat", 1, 1, 0, 40)
    Plogmean3, Tlogmean3, nelogmean3, nlog3, mlog3, nestdlogmean3, Pstdlogmean3, Tstdlogmean3 = f.hydro("/data/cluster/tlebeau/virgo/virgo_xyz_files/virgo_xyz_hydro_l21.dat", 1, 1, 0, 40)
    #mdm, ndm = f.dm("/data/cluster/tlebeau/virgo/virgo_xyz_dm_high_res.dat")
    #d = FortranFile("rad_log_dm.dat", 'r')
    #mdm=d.read_reals()
    #mdm_sum=np.array([np.sum(mdm[0:i]) for i in range(60)])
    #dm=d.read_reals()

    #print('test')
    

    print("Plogmean",Plogmean)

    sys.exit()
    #plt.plot(rlogbin,mdm_sum)
    #plt.show()
    #mdm, ndm = f.dm("/data/cluster/tlebeau/virgo/virgo_xyz_dm_low_res.dat")
    #sys.exit()
    #T, n_e, P, mh, n, Plog, Tlog, nelog, nlog, mlog, nestd, nestdlog, Pstdlog, Tstdlog = f.hydro("/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21_gal_clean.dat",0, 0)
    # plt.plot(rlogbin,Plog,label='P')
    # print('p',Plog)
    # ptest=nelog+Tlog+np.log10(kb)
    # print('ptest',ptest)
    # plt.plot(rlogbin,nelog+Tlog+np.log10(kb),label='$n_ek_BT')
    # plt.show()
    #T, n_e, P, mh, n, Plog, Tlog, nelog, nlog, mlog, nestd, nestdlog, Pstdlog, Tstdlog = f.hydro("/data/cluster/tlebeau/virgo_gal_clean.dat", 1, 1)
    #Plogmean, Tlogmean, nelogmean, nlog, mlog, nestdlogmean, Pstdlogmean, Tstdlogmean = f.hydro("/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l19_low_gal_clean_m1e9.dat", 1, 1)
    #Plogmean, Tlogmean, nelogmean, nlog, mlog, nestdlogmean, Pstdlogmean, Tstdlogmean = f.hydro("/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21_gal_clean_m1e8.5.dat", 1, 1)
    #Plogmed, Tlogmed, nelogmed, nlog, mlog, nestdlogmed, Pstdlogmed, Tstdlogmed = f.hydro("/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l19_low_gal_clean_m1e9.dat", 1, 2)
    # den,den_norm,mtot,mdmtot=f.density(mlog,mdm,1,1024)
    # bxc,bszc=f.ratios(P,Plog,T,Tlog,n_e,nelog,mtot)
    #Plog2, Tlog2, nelog2, nlog2, mlog2, nestdlog2, Pstdlog2,Tstdlog2 = f.hydro("/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l19_low.dat", 0, 1)
    #den, den_norm, mtot, mdmtot = f.density(mlog, mdm, 1, 1024)
    # T3,n_e3,P3,mh3,n3,Plog3,Tlog3,nelog3,nlog3,mlog3,nestd3,nestdlog3,Pstdlog3,Tstdlog3=f.hydro("virgo_rmv_rvir.dat",0,0)
    # print("Tlog3",Tlog3)
    #T4,n_e4,P4,mh4,n4,Plog4,Tlog4,nelog4,nlog4,mlog4,nestd4,nestdlog4,Pstdlog4,Tstdlog4=f.hydro("virgo_hydro_l17.dat",2)
    # return mtot2
    #print("tlogmean",Tlogmean)
    #T_kev=np.array([10**(Tlogmean[i]) for i in range(40)])
    #print("t_kev",T_kev)

    #Plogmean, Tlogmean, nelogmean, nlog, mlog, nestdlogmean, Pstdlogmean, Tstdlogmean = f.hydro("/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l19_low.dat", 1, 1)

    #Plogmean, Tlogmean, nelogmean, nlog, mlog, nestdlogmean, Pstdlogmean, Tstdlogmean = f.hydro("/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l19.dat", 1, 1)
    #Plogmean, Tlogmean, nelogmean, nlog, mlog, nestdlogmean, Pstdlogmean, Tstdlogmean = f.hydro("/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l19_low_gal_clean_m1e9.dat", 1, 1)

    #Plogmean, Tlogmean, nelogmean, nlog, mlog, nestdlogmean, Pstdlogmean, Tstdlogmean = f.hydro_ellipsoid("/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l19_low_gal_clean_m1e9.dat", 1, 1)

    #Plogmean, Tlogmean, nelogmean, nlog, mlog, nestdlogmean, Pstdlogmean, Tstdlogmean = f.hydro("/data/cluster/tlebeau/virgo/virgo_xyz_hydro_l21.dat", 1, 1)


def rad_data_storage(file):
    rad_data = np.zeros((19,40))
    rad_data[0, :] = Plogmed
    rad_data[1, :] = Pstdlogmed
    rad_data[2, :] = Plogmean
    rad_data[3, :] = Pstdlogmean
    rad_data[4, :] = Tlogmed
    rad_data[5, :] = Tstdlogmed
    rad_data[6, :] = Tlogmean
    rad_data[7, :] = Tstdlogmean
    rad_data[8, :] = nelogmed
    rad_data[9, :] = nestdlogmed
    rad_data[10, :] = nelogmean
    rad_data[11, :] = nestdlogmean
    rad_data[12, :] = nlog2
    rad_data[13, :] = mlog2
    rad_data[14, :] = ndm
    rad_data[15, :] = mdm
    rad_data[16, :] = np.array([10**(35+i)*0.05 for i in range(0,40)])
    rad_data[18, :] = np.array([10**(36+i)*0.05 for i in range(0,40)])
    rad_data[17, :] = 0.5*(rad_data[18,:]-rad_data[16,:])+rad_data[16,:]

    np.save(file,rad_data)

#rad_data_storage("lograd_data_lvl19_cleanm1e9.npy")

#print("data storage complete")

#sys.exit()

def load_lograd_dat(file):
    rad_data = np.load(file)
    Plogmed = rad_data[0, :]
    Pstdlogmed = rad_data[1, :]
    Plogmean = rad_data[2, :]
    Pstdlogmean = rad_data[3, :]
    Tlogmed = rad_data[4, :]
    Tstdlogmed = rad_data[5, :]
    Tlogmean = rad_data[6, :]
    Tstdlogmean = rad_data[7, :]
    nelogmed = rad_data[8, :]
    nestdlogmed = rad_data[9, :]
    nelogmean = rad_data[10, :]
    nestdlogmean = rad_data[11, :]
    nlog = rad_data[12, :]
    mlog = rad_data[13, :]
    ndm = rad_data[14, :]
    mdm = rad_data[15, :]
    inf = rad_data[16, :]
    middle = rad_data[17, :]
    sup = rad_data[18, :]
    print("Loading ",file,"complete")
    return Plogmed,Pstdlogmed,Plogmean,Pstdlogmean,Tlogmed,Tstdlogmed,Tlogmean,Tstdlogmean,nelogmed,nestdlogmed,nelogmean,nestdlogmean,nlog,mlog,ndm,mdm,inf,middle,sup


#Plogmed,Pstdlogmed,Plogmean,Pstdlogmean,Tlogmed,Tstdlogmed,Tlogmean,Tstdlogmean,nelogmed,nestdlogmed,nelogmean,nestdlogmean,nlog,mlog,ndm,mdm,inf,middle,sup=load_lograd_dat("lograd_data_lvl21_full.npy")
#Plogmed_ca,Pstdlogmed_ca,Plogmean_ca,Pstdlogmean_ca,Tlogmed_ca,Tstdlogmed_ca,Tlogmean_ca,Tstdlogmean_ca,nelogmed_ca,nestdlogmed_ca,nelogmean_ca,nestdlogmean_ca,nlog_ca,mlog_ca,ndm_ca,mdm_ca,inf_ca,middle_ca,sup_ca=load_lograd_dat("lograd_data_lvl21_cleanall.npy")
#Plogmed_c9,Pstdlogmed_c9,Plogmean_c9,Pstdlogmean_c9,Tlogmed_c9,Tstdlogmed_c9,Tlogmean_c9,Tstdlogmean_c9,nelogmed_c9,nestdlogmed_c9,nelogmean_c9,nestdlogmean_c9,nlog_c9,mlog_c9,ndm_c9,mdm_c9,inf_c9,middle_c9,sup_c9=load_lograd_dat("lograd_data_lvl21_cleanm1e9.npy")
load_data=1
if load_data==1:
    #Plogmed_c8,Pstdlogmed_c8,Plogmean_c8,Pstdlogmean_c8,Tlogmed_c8,Tstdlogmed_c8,Tlogmean_c8,Tstdlogmean_c8,nelogmed_c8,nestdlogmed_c8,nelogmean_c8,nestdlogmean_c8,nlog_c8,mlog_c8,ndm_c8,mdm_c8,inf_c8,middle_c8,sup_c8=load_lograd_dat("lograd_data_lvl21_cleanm1e8.5.npy")
    Plogmed_17,Pstdlogmed_17,Plogmean_17,Pstdlogmean_17,Tlogmed_17,Tstdlogmed_17,Tlogmean_17,Tstdlogmean_17,nelogmed_17,nestdlogmed_17,nelogmean_17,nestdlogmean_17,nlog_17,mlog_17,ndm_17,mdm_17,inf_17,middle_17,sup_17=load_lograd_dat("lograd_data_lvl17_cleanm1e9.npy")
    #Plogmed_17o,Pstdlogmed_17o,Plogmean_17o,Pstdlogmean_17o,Tlogmed_17o,Tstdlogmed_17o,Tlogmean_17o,Tstdlogmean_17o,nelogmed_17o,nestdlogmed_17o,nelogmean_17o,nestdlogmean_17o,nlog_17o,mlog_17o,ndm_17o,mdm_17o,inf_17o,middle_17o,sup_17o=load_lograd_dat("lograd_data_lvl17_cleanm1e9_KPam-3.npy")

    #Plogmed_c8l,Pstdlogmed_c8l,Plogmean_c8l,Pstdlogmean_c8l,Tlogmed_c8l,Tstdlogmed_c8l,Tlogmean_c8l,Tstdlogmean_c8l,nelogmed_c8l,nestdlogmed_c8l,nelogmean_c8l,nestdlogmean_c8l,nlog_c8l,mlog_c8l,ndm_c8l,mdm_c8l,inf_c8l,middle_c8l,sup_c8l=load_lograd_dat("lograd_data_lvl21_cleanm1e8.5_left.npy")
    #Plogmed_c8r,Pstdlogmed_c8r,Plogmean_c8r,Pstdlogmean_c8r,Tlogmed_c8r,Tstdlogmed_c8r,Tlogmean_c8r,Tstdlogmean_c8r,nelogmed_c8r,nestdlogmed_c8r,nelogmean_c8r,nestdlogmean_c8r,nlog_c8r,mlog_c8r,ndm_c8r,mdm_c8r,inf_c8r,middle_c8r,sup_c8r=load_lograd_dat("lograd_data_lvl21_cleanm1e8.5_right.npy")

    Plogmed_c8, Pstdlogmed_c8, Plogmean_c8, Pstdlogmean_c8, Tlogmed_c8, Tstdlogmed_c8, Tlogmean_c8, Tstdlogmean_c8, nelogmed_c8, nestdlogmed_c8, nelogmean_c8, nestdlogmean_c8, nlog_c8, mlog_c8, ndm_c8, mdm_c8, inf_c8, middle_c8, sup_c8 = load_lograd_dat("lograd_data_lvl21_cleanm1e8.5.npy")

    Plogmed_19, Pstdlogmed_19, Plogmean_19, Pstdlogmean_19, Tlogmed_19, Tstdlogmed_19, Tlogmean_19, Tstdlogmean_19, nelogmed_19, nestdlogmed_19, nelogmean_19, nestdlogmean_19, nlog_19, mlog_19, ndm_19, mdm_19, inf_19, middle_19, sup_19 = load_lograd_dat("lograd_data_lvl19_cleanm1e9.npy")

    #print("plogmedc8",Plogmed_c8)
    #print("plogmedc8 test",Plogmed_c8_t)
    #sys.exit()
    #print(nlog_c8)
    #sys.exit()

    ne_urban=np.loadtxt('./obs_data/ne_prof+errb_urban.txt')
    P_urban=np.loadtxt('./obs_data/P_prof+errb_urban.txt')
    ne_planck=np.loadtxt('./obs_data/ne_prof+errb_planck.txt')
    P_planck=np.loadtxt('./obs_data/P_prof+errb_planck.txt')
    ne_ghiz=np.loadtxt('./obs_data/ne_prof_ghizzardi.txt')
    P_planck_deg=np.loadtxt('./obs_data/P_prof_deg_planck.txt')
    #print("deg",P_planck_deg[:,0])
    r_deg21=np.tan(P_planck_deg[:,0]*(np.pi/180))*21e3
    r_deg16=np.tan(P_planck_deg[:,0]*(np.pi/180))*16e3
    #print("rdeg",r_deg)

    m_sum = np.array([np.sum(mdm_c8[0:i + 1]) + np.sum(mlog_c8[0:i + 1]) for i in range(40)])
    #print("mdmc8",mdm_c8)
    #sys.exit()

    k_data=np.load("K_prof_high.npy")
    k_h=k_data[0,:]
    k_h_std=k_data[1,:]
    k_data = np.load("K_prof_low.npy")
    k_l = k_data[0, :]
    k_l_std = k_data[1, :]

    #Plogmed_c8-=1.136842

#f.bias_w_mc(Plogmean_c8,Pstdlogmean_c8,nelogmean_c8,nestdlogmean_c8,Tlogmean_c8,Tstdlogmean_c8,m_sum,10000)
#sys.exit()
#end=time.time()
#print("data read time=", end-start)

# Tlogmed-=np.log10(kb/1.60218e-16)
# Tstdlogmed-=np.log10(kb/1.60218e-16)
# Tlogmean-=np.log10(kb/1.60218e-16)
# Tstdlogmean-=np.log10(kb/1.60218e-16)

# T=(T*kb)/(1.6e-16)
# T2=(T2*kb)/(1.6e-16)
# T3=(T3*kb)/(1.6e-16)
# plt.scatter(rlogbin,Plog,s=6)
# plt.errorbar(rlogbin,Plog,yerr=Pstdlog,ls='solid',label='3D profile')
# plt.scatter(rlogbin,Plog2,s=6)
# plt.errorbar(rlogbin,Plog2,yerr=Pstdlog2,ls='dashed',label='All cells',alpha=0.8)
# plt.xscale('log')
# plt.xlabel('R (kpc)')
# plt.ylabel('$log_{10} (P[Pa])$')
# plt.legend()
# plt.show()

# end=time.time()

# print("temps:", end-start)

def compar_elena(Plog):
    print(Plogmed_17)
    P = np.array([10 ** Plogmed_17[i] for i in range(40)])
    #P /= (1.60218e-10)
    print(P)

    Pelena = np.loadtxt('./elena/virgo_press_prof_elena.txt')
    Pobs = np.loadtxt('./elena/virgo_obs_elena.txt')
    #print(Pelena.T)

    #plt.plot(rlogbin, Plogmed_17, '.', label='RAMSES (Theo)')
    plt.plot(Pelena.T[1], np.log10(Pelena.T[0]), '.', label='GADGET (Elena)')
    plt.plot(Pobs.T[0], np.log10(Pobs.T[1]), '.', label='Observation (Elena)')
    #plt.legend()
    #plt.show()


#compar_elena(Plogmed_c8)

#sys.exit()
def bias():
    # den,den_norm,mtot2,mdmtot2=f.density(mlog2,mdm,1,1024)
    # den,den_norm,mtot3,mdmtot3=f.density(mlog3,mdm,1,1024)
    # bx7,bsz7=f.ratios(P3,Plog3,T3,Tlog3,n_e3,nelog3,mtot3)
    #Plogp = Plog + Pstdlog
    #Plogm = Plog - Pstdlog
    #Tlogp = Tlog + Tstdlog
    #Tlogm = Tlog - Tstdlog
    #nelogp = nelog + nestdlog
    #nelogm = nelog - nestdlog

    #bxc, bszc = f.ratios(Plog, Tlog, nelog, mtot2)
    #bxcp, bszcp = f.ratios(Plog, Tlogp, nelogp, mtot2)
    #bxcm, bszcm = f.ratios(Plogm, Tlogm, nelogm, mtot2)

    #Plogpm = Plogmean + Pstdlogmean
    #Plogmm = Plogmean - Pstdlogmean
    #Tlogpm = Tlogmean + Tstdlogmean
    #Tlogmm = Tlogmean - Tstdlogmean
    #nelogpm = nelogmean + nestdlogmean
    #nelogmm = nelogmean - nestdlogmean

    r500=1087
    #print("rlogbin_cen",rlogbin_cen[25])
    Plogmean_c8_unlog=10**Plogmean_c8
    P500=Plogmean_c8_unlog[25]


    print("P_500",P500)

    P500=3.65e-3

    print("P_500",P500)

    pratio_log=np.log10(Plogmean_c8_unlog/P500)

    pratio=Plogmean_c8_unlog/P500



    #sys.exit()
    #popt,pcov=curve_fit(f.upplog,rlogbin_cen/r500,pratio_log,p0=[6.41,1.81,1.33,4.13,0.31],maxfev=10000)#bounds=([0,0,0,0,-10],[1000,20,10,100,10]),method='trf',)
    popt, pcov = curve_fit(f.upp, rlogbin_cen / r500, pratio,p0=[6.41,1.81,1.33,4.13,0.31],maxfev=100000)#,bounds=([0,0,0,0,-10],[600,20,10,50,10]),method='trf',maxfev=10000)
    #popt, pcov = curve_fit(f.upplog, rlogbin_cen, Plogmean_c8, p0=[6.41, 1.81, 1.33, 4.13, 0.31],bounds=([0, 0, 0, 0, -10], [600, 20, 10, 15, 10]), method='trf')
    print("popt p",popt)
    print("sigma p",np.sqrt(np.diagonal(pcov)))

    plt.plot(rlogbin_cen/r500, pratio, label="3D, $8192^3$", ls='dashed', marker='.', color='orange')
    plt.plot(rlogbin_cen/r500, f.upp(rlogbin_cen/r500,*popt), label='gNFW fit', color='black')

    Plogmean = np.load("P_prof_l21.npy")
    r500 = 1087
    rbin = rlogbin_cen / r500
    n=40
    p500_1 = 1.66e-3
    p500_2 = 3.652e-3
    Plogmean_1 = 10 ** Plogmean / p500_1
    Plogmean_2 = 10 ** Plogmean / p500_2
    p13 = [f.upp(rbin[i], 6.41, 1.81, 1.33, 4.13, 0.31) for i in range(0, n)]
    a10 = [f.upp(rlogbin_cen[i] / r500, 8.40, 1.18, 1.05, 5.49, 0.31) for i in range(0, n)]
    p7 = [f.upp(rlogbin_cen[i] / r500, 3.36, 1.18, 1.08, 4.30, 0.31) for i in range(0, n)]
    a85p = [f.upp(rlogbin_cen[i] / r500, 5.99, 0.02, 0.48, 14.97, 0.31) for i in range(0, n)]
    ghir19 = [f.upp(rlogbin_cen[i] / r500, 5.68, 1.49, 1.33, 4.05, 0.29) for i in range(0, n)]

    #plt.plot(rlogbin_cen / r500, Plogmean_1, label=r'Virgo, $P_{500}=3.65\times 10^{-3}~keV~cm^{-3}$' + '\n (mass-w pressure within $R_{500}$)', marker='.', color='black')
    #plt.plot(rlogbin_cen / r500, Plogmean_2, label=r'Virgo, $P_{500}=1.66e\times 10^{-3}~keV~cm^{-3}$' + '\n (Ghirardini+19 scaling relation)', marker='.', color='grey')
    #plt.plot(rlogbin_cen / r500, ghir19, label='Ghirardini+ 19', color='red')
    #plt.plot(rlogbin_cen / r500, a10, label='Arnaud+ 10', color='blue')
    # plt.plot(rlogbin_cen / r500, p13,label='Planck+ 13')
    plt.plot(rlogbin_cen / r500, p7, label='PACT 21', color='forestgreen')

    plt.legend()
    plt.yscale('log')
    plt.xscale('log')
    plt.show()

    sys.exit()

    popt_ne,pcov_ne=curve_fit(f.log_beta_model,rlogbin_cen,nelogmean_c8)
    print("popt ne",popt_ne)
    print("sigma ne",np.sqrt(np.diagonal(pcov_ne)))

    alpha = np.load('alpha.npy')

    popt_alpha,pcov_alpha=curve_fit(f.alpha_model,rlogbin_cen,alpha)

    #plt.plot(rlogbin_cen,alpha)
    #plt.plot(rlogbin_cen,f.alpha_model(rlogbin_cen,*popt_alpha))
    #plt.xscale("log")
    #plt.legend()
    #plt.show()
    #sys.exit()




    def plot_pressure_fit():
        a10=[8.40, 1.18, 1.05, 5.49, 0.31]
        plt.plot(rlogbin_cen/r500, pratio, label="3D, $8192^3$", ls='dashed', marker='.',color='orange')
        plt.plot(rlogbin_cen/r500,f.upplog(rlogbin_cen/r500,*popt),label='gNFW fit',color='black')
        plt.plot(rlogbin_cen / r500, f.upplog(rlogbin_cen / r500, *a10), label='a10', color='blue')
        plt.xscale('log')
        plt.xlabel('$R$/$R_{500}$', size=16)
        plt.ylabel("lo$\mathrm{g_{10}}$($P$[keV.c$\mathrm{m^{-3}}$])", size=16)
        plt.ylabel('$P$/$P_{500}$')
        plt.axvline(x=1, color='grey')  # label='$R_{500}$')
        plt.axvline(x=2147/1087, color='grey', ls='dashed')  # label='$R_{Vir}$')
        plt.text(1, 1, "$R_{500}$", rotation=90, size=16)
        plt.text(2147/1087, 1, "$R_{vir}$", rotation=90, size=16)
        plt.legend()
        plt.show()
        sys.exit()

    def plot_ne_fit():
        plt.plot(rlogbin_cen, nelogmean_c8, label="3D, $8192^3$", ls='dashed', marker='.',color='orange')
        plt.plot(rlogbin_cen,f.log_beta_model(rlogbin_cen,*popt_ne),label='$beta$-model fit',color='black')
        plt.xscale('log')
        plt.xlabel('R (kpc)', size=16)
        #plt.ylabel("lo$\mathrm{g_{10}}$($P$[keV.c$\mathrm{m^{-3}}$])", size=16)
        #plt.ylabel('$P$/$P_{500}$')
        plt.ylabel("lo$\mathrm{g_{10}}$($n_{\mathrm{e}}$[c$\mathrm{m^{-3}}$])", size=16)
        plt.axvline(x=1087, color='grey')  # label='$R_{500}$')
        plt.axvline(x=2147, color='grey', ls='dashed')  # label='$R_{Vir}$')
        plt.text(1087, -3.5, "$R_{500}$", rotation=90, size=16)
        plt.text(2147, -3.5, "$R_{vir}$", rotation=90, size=16)
        plt.legend()
        plt.show()
        sys.exit()

    plot_pressure_fit()
    #plot_ne_fit()

    Plogmean_c8_fit=np.log10(10**(f.upplog(rlogbin_cen/r500,*popt))*P500)
    nelogmean_c8_fit=f.log_beta_model(rlogbin_cen,*popt_ne)
    alpha_fit=f.alpha_model(rlogbin_cen,*popt_alpha)
    #print(Plogmean_c8_fit)
    #sys.exit()
    #plt.plot(rlogbin_cen/r500,Plogmean_c8,label="3D")
    #plt.plot(rlogbin_cen/r500,Plogmean_c8_fit,label='fit')
    #plt.plot(rlogbin_cen,nelogmean_c8,label="3D")
    #plt.plot(rlogbin_cen,nelogmean_c8_fit,label='fit')
    #plt.legend()
    #plt.xscale('log')
    #plt.show()
    #sys.exit()

    m_sum = np.array([np.sum(mdm_c8[0:i + 1]) + np.sum(mlog_c8[0:i + 1]) for i in range(40)])
    #print("m500 dm",np.sum(mdm_c8[0:32]))
    print("mtot",m_sum)

    m_cumul_dm = np.load("m_cumul_dm.npy")
    m_cumul_ba = np.load("m_cumul_ba.npy")
    m_cumul_sum = m_cumul_dm + m_cumul_ba

    m_sum=m_cumul_sum


    print("rlogbincen",rlogbin_cen[31])
    print("m_cumul_sum",m_cumul_sum[31])
    print("m_cumul_dm",m_cumul_dm[31])
    print("m_cumul_ba",m_cumul_ba[31])

    #sys.exit()



    #plt.plot(rlogbin_cen,m_cumul_sum,'.', label="m_cumul_sum")
    #plt.plot(rlogbin_cen,m_sum,'.', label="m sum old")
    #plt.legend()
    #plt.show()
    #sys.exit()

    #np.savetxt("radial_mtot_virgo.txt",m_sum)
    #np.savetxt("radial_mhydro_virgo.txt",mlog_c8)
    #np.savetxt("radial_mdm_virgo.txt", mdm_c8)

    print("rlogbincen",rlogbin_cen)
    print("msum",m_sum)
    #sys.exit()
    #print("rlogbincen",rlogbin_cen[25])
    #sys.exit()
    #plt.plot(rlogbin_cen,m_sum)
    #plt.show()
    #sys.exit()
    #m_suml = np.array([np.sum(mdm_c8l[0:i + 1]) + np.sum(mlog_c8l[0:i + 1]) for i in range(40)])
    #m_sumr = np.array([np.sum(mdm_c8r[0:i + 1]) + np.sum(mlog_c8r[0:i + 1]) for i in range(40)])
    #m_sum17 = np.array([np.sum(mdm_17[0:i + 1]) + np.sum(mlog_17[0:i + 1]) for i in range(40)])
    #m_sum17o = np.array([np.sum(mdm_17o[0:i + 1]) + np.sum(mlog_17o[0:i + 1]) for i in range(40)])
    print("no fit")
    bx, bsz, smx, smsz, bcorr = f.ratios(Plogmean_c8, Tlogmean_c8, nelogmean_c8, m_sum,1,Pstdlogmean_c8,nestdlogmean_c8,Tstdlogmean_c8,alpha)
    print("fit p")
    bxfit, bszfit, smxfit, smszfit, bcorrfit = f.ratios(Plogmean_c8_fit, Tlogmean_c8, nelogmean_c8, m_sum, 1, Pstdlogmean_c8, nestdlogmean_c8,Tstdlogmean_c8,alpha)
    print("fit p + ne")
    bxfit2, bszfit2, smxfit2, smszfit2, bcorrfit2 = f.ratios(Plogmean_c8_fit, Tlogmean_c8, nelogmean_c8_fit, m_sum, 1, Pstdlogmean_c8, nestdlogmean_c8, Tstdlogmean_c8,alpha_fit)
    #bxl, bszl, smxl, smszl = f.ratios(Plogmean_c8l, Tlogmean_c8l, nelogmean_c8l, m_sum, 1,Pstdlogmean_c8l,nestdlogmean_c8l,Tstdlogmean_c8l)
    #bxr, bszr, smxr, smszr = f.ratios(Plogmean_c8r, Tlogmean_c8r, nelogmean_c8r, m_sum, 1,Pstdlogmean_c8r,nestdlogmean_c8r,Tstdlogmean_c8r)
    #plt.show()
    #bxt, bszt = f.ratios(Plogmean_c8_t, Tlogmean_c8, nelogmean_c8, m_sum, 1)
    #bx17, bsz17 = f.ratios(Plogmean_17, Tlogmean_17, nelogmean_17, m_sum17,1)
    #bx17med, bsz17med = f.ratios(Plogmed_17, Tlogmed_17, nelogmed_17, m_sum17, 1)
    #bx17o, bsz17o = f.ratios(Plogmean_17o, Tlogmean_17o, nelogmean_17o, m_sum17o, 0)
    #bxp, bszp = f.ratios(Plogp, Tlogp, nelogp, mtot2)
    #bxm, bszm = f.ratios(Plogm, Tlogm, nelogm, mtot2)
    #print("plogmean",Plogmean_c8)
    #print("nelogmean*tlogmean",Tlogmean_c8*nelogmean_c8)
    #print("ratio",Plogmean_c8/(Tlogmean_c8*nelogmean_c8))

    #sys.exit()

    #bxmed, bszmed = f.ratios(Plogmed, Tlogmed, nelogmed, mtot2)

    #np.savetxt("fitted_3D_bias",bszfit2)
    #sys.exit()

    #np.save("bsz_virgo_3D_h.npy",bsz)
    #print("3D bias saved")

    p_rad_test = np.load('p_rad_test.npy')
    ne_rad_test = np.load('ne_rad_test.npy')

    print("plogmeanc8",Plogmean_c8)
    print("p_rad_test",p_rad_test)

    print("nelogmeanc8", nelogmean_c8)
    print("ne_rad_test", ne_rad_test)

    bx, bsz_test, smx, smsz, bcorr = f.ratios(p_rad_test, Tlogmean_c8, ne_rad_test, m_sum, 1, Pstdlogmean_c8, nestdlogmean_c8, Tstdlogmean_c8, alpha)

    print("bsz",bsz)
    print("bsz_test",bsz_test)


    plt.plot(rlogbin_cen,bsz,label='original')
    plt.plot(rlogbin_cen,bsz_test,label='test')
    plt.legend()
    plt.xscale('log')
    plt.xlim(200, 4500)
    plt.ylim(0.75, 2.45)
    plt.axvline(x=1087, color='grey')  # label='$R_{500}$')
    plt.axvline(x=2147, color='grey', ls='dashed')  # label='$R_{Vir}$')
    plt.text(1087, 0.83, "$R_{500}$", rotation=90, size=16)
    plt.text(2147, 0.83, "$R_{vir}$", rotation=90, size=16)
    plt.axhline(y=1)
    plt.show()


    sys.exit()

    #np.save("fitted_3D_bias.npy",bszfit2)

    end = time.time()

    #print("time=", end - start)

    # r=[100*(i+1) for i in range(39)]
    # r=np.array(r)
    # plt.plot(rlogbin,bx,label='from T, all cells, lvl 17')
    # plt.plot(rlogbin,bsz,label='from P, all cells,lvl 13 ')
    # plt.plot(rlogbin,bx5,label='from T, $T>10^5$, lvl 17')
    # plt.plot(rlogbin,bsz5,label='from P,$T>10^5$, lvl 13')
    #print("plogmed17",Plogmed_17)
    #print("plogmed17o",Plogmed_17o)
    #print(bx)
    #print(bsz)
    #print(bx17)
    #print(bx17o)
    #print("ratio",bszt/bxt)
    #plt.plot(rlogbin, bxt, label='From T and $n_e$ profile with weighted mean')
    #plt.plot(rlogbin, bszt, label='From P profile with weighted mean')

    #plt.plot(rlogbin_cen, bsz,label="From P profile with weighted mean",color='blue')
    def fit_plot():
        r500=1087
        un = [1 for i in range(40)]
        zerohuit = [0.8 for i in range(40)]
        un = np.array(un)
        zerohuit = np.array(zerohuit)

        f, (a, b) = plt.subplots(2, 1)

        a.plot(rlogbin_cen, bsz, label="3D, no fit", color='orange',marker='.')
        a.plot(rlogbin_cen, bszfit, label="3D, fitted pressure profile (gNFW)", color='darkgreen',marker='.',ls='dashed')
        a.plot(rlogbin_cen, bszfit2, label="3D, fitted pressure (gNFW) + electron density(beta-model)", color='darkgreen',marker='.',ls='dotted')
        a.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        a.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)
        #a.plot(rlogbin, zerohuit, color='grey')
        #a.plot(rlogbin, un, color='black')
        a.set_xscale('log')
        a.legend(prop={'size': 16})
        a.set_ylabel("$(1-\mathrm{b})=rac{M_{\mathrm{HE}}}{M_{\mathrm{tot}}}$")

        a.axvline(x=1087, color='grey')  # label='$R_{500}$')
        a.axvline(x=2147, color='grey', ls='dashed')  # label='$R_{Vir}$')
        a.text(1087, 0.83, "$R_{500}$", rotation=90, size=16)
        a.text(2147, 0.83, "$R_{vir}$", rotation=90, size=16)
        a.plot(rlogbin, zerohuit, color='grey')
        a.plot(rlogbin, un, color='black')
        a.set_xlim(200, 4500)
        a.set_ylim(0.75, 2.45)



        b.plot(rlogbin_cen, bsz, label="3D, no fit", color='orange', marker='.')
        b.plot(rlogbin_cen, bcorr, label="Non-thermal pressure correction, no fit",color='firebrick',marker='.',ls='dashed')
        b.plot(rlogbin_cen, bcorrfit2, label="Non-thermal pressure correction, fit ", color='firebrick',ls='dotted',marker='.')
        b.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        b.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)
        #b.plot(rlogbin, zerohuit, color='grey')
        #b.plot(rlogbin, un, color='black')
        b.set_xscale('log')
        b.legend(prop={'size': 16})
        b.set_ylabel("$(1-\mathrm{b})=rac{M_{\mathrm{HE}}}{M_{\mathrm{tot}}}$")
        b.set_xlabel("$R [\mathrm{kpc}]$")

        b.axvline(x=1087, color='grey')  # label='$R_{500}$')
        b.axvline(x=2147, color='grey', ls='dashed')  # label='$R_{Vir}$')
        b.text(1087, 0.83, "$R_{500}$", rotation=90, size=16)
        b.text(2147, 0.83, "$R_{vir}$", rotation=90, size=16)
        b.plot(rlogbin, zerohuit, color='grey')
        b.plot(rlogbin, un, color='black')
        b.set_xlim(200, 4500)
        b.set_ylim(0.75, 2.45)

        # plt.axvline(x=950, color='grey', ls='dotted')# label='$R_{bump}$')
        # plt.text(950, 1.1, "$R_{bump}$", rotation=90, size=16)


        plt.subplots_adjust(wspace=0, hspace=0)
        plt.show()
        print("end")
        sys.exit()
    #plt.plot(rlogbin_cen, bszl, label="Only cells with x<x_Virgo", color='orange')
    #plt.plot(rlogbin_cen, bszr, label="Only cells with x>x_Virgo", color='green')
    #plt.errorbar(rlogbin_cen, bsz, yerr=smsz, label="From P profile with weighted mean")
    #plt.fill_between(rlogbin_cen,bsz-smsz,bsz+smsz,color="blue",alpha=0.1)
    #plt.plot(rlogbin_cen, bx,label="From T and ne profiles with weighted mean",alpha=0.5,marker='+')
    #plt.fill_between(rlogbin_cen,bx-smx,bx+smx,color='orange',alpha=0.2)
    #plt.errorbar(rlogbin_cen, bx, yerr=smx, label="From T and ne profiles with weighted mean")

    #plt.plot(rlogbin, bxg, label='From T and $n_e$ profile with weighted mean and median respectively')
    #lt.plot(rlogbin, bszg / 1.136818, label='From P profile with median')
    #plt.plot(rlogbin, bx17, label='From T and $n_e$ profiles with weighted mean, low res',c="blue",alpha=0.8)

    #plt.plot(rlogbin, bxl, label='From T and $n_e$ profile with weighted mean, only cells with x<x_virgo')
    #plt.plot(rlogbin, bxr, label='From T and $n_e$ profile with weighted mean, only cells with x>x_virgo')
    #plt.plot(rlogbin, bx17o, label='From T and $n_e$ profiles with weighted mean, low res K Pa m-3 units',ls="dashed",alpha=0.8,c="orange")
    #plt.plot(rlogbin, bsz17/1.136818, label='From P profile with weighted mean, low res', c="purple", alpha=0.8)
    #plt.plot(rlogbin, bsz17o/1.136818, label='From P profile with weighted mean, low res K Pa m-3 units', ls="dashed",alpha=0.8, c="red")
    #plt.plot(rlogbin, bsz17, label='From P profile with weighted mean, low res')


    #m_dmsum = np.array([np.sum(mdm_c8[0:i]) for i in range(40)])
    #plt.plot(rlogbin,m_sum,label="total mass: DM+baryons")
    #plt.plot(rlogbin, m_sum17, label="total mass: DM+baryons low res")
    #plt.plot(rlogbin, m_dmsum, label="total mass: DM only")
    # plt.plot(rlogbin,bxc,label='From T and $n_e$ profiles, median for ne, mean for T')
    # plt.plot(rlogbin,bsz,label='From P profile with weighted mean')
    # plt.plot(rlogbin,bszc,label='From P profile  with median')
    # plt.plot(rlogbin,bxmed,label='From T and $n_e$ profile with median')
    # plt.plot(rlogbin,bszc,label='From P profile  with median')
    # plt.axvline(x=1087,color='grey',label='$R_{500}$')
    # plt.axvline(x=2024,color='grey',ls='dashed',label='$R_{Vir}$')

    # print("bx",bx)
    # print("bxc",bxc)

    # plt.fill_between(rlogbin,bxm,bxp,color='blue',alpha=0.5)
    # plt.fill_between(rlogbin,bszm,bszp,color='orange',alpha=0.5)
    #np.savetxt('bsz_virgo_3D_h.txt',bsz)

    fit_plot()


    plt.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
    plt.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)

    un=[1 for i in range(40)]
    zerohuit=[0.8 for i in range(40)]
    un=np.array(un)
    zerohuit=np.array(zerohuit)
    plt.plot(rlogbin,zerohuit,color='grey')
    plt.plot(rlogbin,un,color='black')
    plt.xscale('log')
    plt.xlabel('$\mathrm{R/R_{500}}$')
    plt.legend()
    plt.ylabel("$(1-\mathrm{b})=rac{M_{\mathrm{HE}}}{M_{\mathrm{tot}}}$")
    #plt.ylabel('$(1-b)=\frac{M_{HE}}{M_{tot}}$')

    #plt.axvline(x=1087, color='grey')# label='$R_{500}$')
    #plt.axvline(x=2147, color='grey', ls='dashed')# label='$R_{Vir}$')
    #plt.text(1087, 1.1, "$R_{500}$", rotation=90, size=16)
    #plt.text(2147, 1.1, "$R_{vir}$", rotation=90, size=16)

    plt.axvline(x=1, color='grey')  # label='$R_{500}$')
    plt.axvline(x=2147 / 1087, color='grey', ls='dashed')  # label='$R_{Vir}$')
    plt.text(1, 1.1, "$R_{500}$", rotation=90, size=16)
    plt.text(2147 / 1087, 1.1, "$R_{vir}$", rotation=90, size=16)
    plt.plot(rlogbin/r500, zerohuit, color='grey')
    plt.plot(rlogbin/r500, un, color='black')

    #plt.axvline(x=950, color='grey', ls='dotted')# label='$R_{bump}$')
    #plt.text(950, 1.1, "$R_{bump}$", rotation=90, size=16)
    plt.xlim(200/r500, 4500/r500)
    plt.ylim(0.75,2.1)
    #plt.ylim(0.3, 2.5)
    #plt.grid()
    plt.show()
    sys.exit()

#bias()

def bias_low_res():
    # den,den_norm,mtot2,mdmtot2=f.density(mlog2,mdm,1,1024)
    # den,den_norm,mtot3,mdmtot3=f.density(mlog3,mdm,1,1024)
    # bx7,bsz7=f.ratios(P3,Plog3,T3,Tlog3,n_e3,nelog3,mtot3)
    #Plogp = Plog + Pstdlog
    #Plogm = Plog - Pstdlog
    #Tlogp = Tlog + Tstdlog
    #Tlogm = Tlog - Tstdlog
    #nelogp = nelog + nestdlog
    #nelogm = nelog - nestdlog

    #bxc, bszc = f.ratios(Plog, Tlog, nelog, mtot2)
    #bxcp, bszcp = f.ratios(Plog, Tlogp, nelogp, mtot2)
    #bxcm, bszcm = f.ratios(Plogm, Tlogm, nelogm, mtot2)

    #Plogpm = Plogmean + Pstdlogmean
    #Plogmm = Plogmean - Pstdlogmean
    #Tlogpm = Tlogmean + Tstdlogmean
    #Tlogmm = Tlogmean - Tstdlogmean
    #nelogpm = nelogmean + nestdlogmean
    #nelogmm = nelogmean - nestdlogmean

    r500=1087
    #print("rlogbin_cen",rlogbin_cen[25])
    Plogmean_17_unlog=10**Plogmean_19
    P500=Plogmean_17_unlog[25]
    pratio=np.log10(Plogmean_17_unlog/P500)



    #sys.exit()
    popt,pcov=curve_fit(f.upplog,rlogbin_cen/r500,pratio,p0=[6.41,1.81,1.33,4.13,0.31],bounds=([0,0,0,0,-10],[600,20,10,50,10]),method='trf')
    #popt, pcov = curve_fit(f.upplog, rlogbin_cen, Plogmean_c8, p0=[6.41, 1.81, 1.33, 4.13, 0.31],bounds=([0, 0, 0, 0, -10], [600, 20, 10, 15, 10]), method='trf')
    print("popt p",popt)
    print("sigma p", np.sqrt(np.diagonal(pcov)))

    popt_ne,pcov_ne=curve_fit(f.log_beta_model,rlogbin_cen,nelogmean_19)
    print("popt ne",popt_ne)
    print("sigma ne", np.sqrt(np.diagonal(pcov_ne)))

    #sys.exit()

    alpha = np.load('alpha.npy')

    popt_alpha,pcov_alpha=curve_fit(f.alpha_model,rlogbin_cen,alpha)

    #plt.plot(rlogbin_cen,alpha)
    #plt.plot(rlogbin_cen,f.alpha_model(rlogbin_cen,*popt_alpha))
    #plt.xscale("log")
    #plt.legend()
    #plt.show()
    #sys.exit()




    def plot_pressure_fit():
        plt.plot(rlogbin_cen/r500, pratio, label="3D, $8192^3$", ls='dashed', marker='.',color='orange')
        plt.plot(rlogbin_cen/r500,f.upplog(rlogbin_cen/r500,*popt),label='gNFW fit',color='black')
        plt.xscale('log')
        plt.xlabel('$R$/$R_{500}$', size=16)
        plt.ylabel("lo$\mathrm{g_{10}}$($P$[keV.c$\mathrm{m^{-3}}$])", size=16)
        plt.ylabel('$P$/$P_{500}$')
        plt.axvline(x=1, color='grey')  # label='$R_{500}$')
        plt.axvline(x=2147/1087, color='grey', ls='dashed')  # label='$R_{Vir}$')
        plt.text(1, 1, "$R_{500}$", rotation=90, size=16)
        plt.text(2147/1087, 1, "$R_{vir}$", rotation=90, size=16)
        plt.legend()
        plt.show()
        sys.exit()

    def plot_ne_fit():
        plt.plot(rlogbin_cen, nelogmean_c8, label="3D, $8192^3$", ls='dashed', marker='.',color='orange')
        plt.plot(rlogbin_cen,f.log_beta_model(rlogbin_cen,*popt_ne),label='$beta$-model fit',color='black')
        plt.xscale('log')
        plt.xlabel('R (kpc)', size=16)
        #plt.ylabel("lo$\mathrm{g_{10}}$($P$[keV.c$\mathrm{m^{-3}}$])", size=16)
        #plt.ylabel('$P$/$P_{500}$')
        plt.ylabel("lo$\mathrm{g_{10}}$($n_{\mathrm{e}}$[c$\mathrm{m^{-3}}$])", size=16)
        plt.axvline(x=1087, color='grey')  # label='$R_{500}$')
        plt.axvline(x=2147, color='grey', ls='dashed')  # label='$R_{Vir}$')
        plt.text(1087, -3.5, "$R_{500}$", rotation=90, size=16)
        plt.text(2147, -3.5, "$R_{vir}$", rotation=90, size=16)
        plt.legend()
        plt.show()
        sys.exit()

    #plot_pressure_fit()
    #plot_ne_fit()

    Plogmean_17_fit=np.log10(10**(f.upplog(rlogbin_cen/r500,*popt))*P500)
    nelogmean_17_fit=f.log_beta_model(rlogbin_cen,*popt_ne)
    alpha_fit=f.alpha_model(rlogbin_cen,*popt_alpha)
    #print(Plogmean_c8_fit)
    #sys.exit()
    #plt.plot(rlogbin_cen/r500,Plogmean_c8,label="3D")
    #plt.plot(rlogbin_cen/r500,Plogmean_c8_fit,label='fit')
    #plt.plot(rlogbin_cen,nelogmean_c8,label="3D")
    #plt.plot(rlogbin_cen,nelogmean_c8_fit,label='fit')
    #plt.legend()
    #plt.xscale('log')
    #plt.show()
    #sys.exit()

    m_sum = np.array([np.sum(mdm_c8[0:i + 1]) + np.sum(mlog_c8[0:i + 1]) for i in range(40)])
    #print("msum",m_sum[25],np.log10(m_sum[25]))
    #print("rlogbincen",rlogbin_cen[25])
    #sys.exit()
    #plt.plot(rlogbin_cen,m_sum)
    #plt.show()
    #sys.exit()
    #m_suml = np.array([np.sum(mdm_c8l[0:i + 1]) + np.sum(mlog_c8l[0:i + 1]) for i in range(40)])
    #m_sumr = np.array([np.sum(mdm_c8r[0:i + 1]) + np.sum(mlog_c8r[0:i + 1]) for i in range(40)])
    #m_sum17 = np.array([np.sum(mdm_17[0:i + 1]) + np.sum(mlog_17[0:i + 1]) for i in range(40)])
    #m_sum17o = np.array([np.sum(mdm_17o[0:i + 1]) + np.sum(mlog_17o[0:i + 1]) for i in range(40)])
    print("no fit")
    bx, bsz, smx, smsz, bcorr = f.ratios(Plogmean_17, Tlogmean_17, nelogmean_17, m_sum,1,Pstdlogmean_17,nestdlogmean_17,Tstdlogmean_17,alpha)
    print("fit p")
    bxfit, bszfit, smxfit, smszfit, bcorrfit = f.ratios(Plogmean_17_fit, Tlogmean_17, nelogmean_17, m_sum, 1, Pstdlogmean_17, nestdlogmean_17,Tstdlogmean_17,alpha)
    print("fit p + ne")
    bxfit2, bszfit2, smxfit2, smszfit2, bcorrfit2 = f.ratios(Plogmean_17_fit, Tlogmean_17, nelogmean_17_fit, m_sum, 1, Pstdlogmean_17, nestdlogmean_17, Tstdlogmean_17,alpha_fit)
    #bxl, bszl, smxl, smszl = f.ratios(Plogmean_c8l, Tlogmean_c8l, nelogmean_c8l, m_sum, 1,Pstdlogmean_c8l,nestdlogmean_c8l,Tstdlogmean_c8l)
    #bxr, bszr, smxr, smszr = f.ratios(Plogmean_c8r, Tlogmean_c8r, nelogmean_c8r, m_sum, 1,Pstdlogmean_c8r,nestdlogmean_c8r,Tstdlogmean_c8r)
    #plt.show()
    #bxt, bszt = f.ratios(Plogmean_c8_t, Tlogmean_c8, nelogmean_c8, m_sum, 1)
    #bx17, bsz17 = f.ratios(Plogmean_17, Tlogmean_17, nelogmean_17, m_sum17,1)
    #bx17med, bsz17med = f.ratios(Plogmed_17, Tlogmed_17, nelogmed_17, m_sum17, 1)
    #bx17o, bsz17o = f.ratios(Plogmean_17o, Tlogmean_17o, nelogmean_17o, m_sum17o, 0)
    #bxp, bszp = f.ratios(Plogp, Tlogp, nelogp, mtot2)
    #bxm, bszm = f.ratios(Plogm, Tlogm, nelogm, mtot2)
    #print("plogmean",Plogmean_c8)
    #print("nelogmean*tlogmean",Tlogmean_c8*nelogmean_c8)
    #print("ratio",Plogmean_c8/(Tlogmean_c8*nelogmean_c8))

    #sys.exit()

    #bxmed, bszmed = f.ratios(Plogmed, Tlogmed, nelogmed, mtot2)

    #np.savetxt("fitted_3D_bias",bszfit2)
    #sys.exit()

    end = time.time()

    #print("time=", end - start)

    # r=[100*(i+1) for i in range(39)]
    # r=np.array(r)
    # plt.plot(rlogbin,bx,label='from T, all cells, lvl 17')
    # plt.plot(rlogbin,bsz,label='from P, all cells,lvl 13 ')
    # plt.plot(rlogbin,bx5,label='from T, $T>10^5$, lvl 17')
    # plt.plot(rlogbin,bsz5,label='from P,$T>10^5$, lvl 13')
    #print("plogmed17",Plogmed_17)
    #print("plogmed17o",Plogmed_17o)
    #print(bx)
    #print(bsz)
    #print(bx17)
    #print(bx17o)
    #print("ratio",bszt/bxt)
    #plt.plot(rlogbin, bxt, label='From T and $n_e$ profile with weighted mean')
    #plt.plot(rlogbin, bszt, label='From P profile with weighted mean')

    #plt.plot(rlogbin_cen, bsz,label="From P profile with weighted mean",color='blue')
    def fit_plot():
        r500=1087
        un = [1 for i in range(40)]
        zerohuit = [0.8 for i in range(40)]
        un = np.array(un)
        zerohuit = np.array(zerohuit)

        f, (a, b) = plt.subplots(2, 1)

        a.plot(rlogbin_cen, bsz, label="3D, no fit", color='orange',marker='.')
        a.plot(rlogbin_cen, bszfit, label="3D, fitted pressure profile (gNFW)", color='gold',marker='.',ls='dashed')
        a.plot(rlogbin_cen, bszfit2, label="3D, fitted pressure (gNFW) + electron density(beta-model)", color='gold',marker='.',ls='dotted')
        a.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        a.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)
        #a.plot(rlogbin, zerohuit, color='grey')
        #a.plot(rlogbin, un, color='black')
        a.set_xscale('log')
        a.legend(prop={'size': 16})
        a.set_ylabel("$(1-\mathrm{b})=rac{M_{\mathrm{HE}}}{M_{\mathrm{tot}}}$")

        a.axvline(x=1087, color='grey')  # label='$R_{500}$')
        a.axvline(x=2147, color='grey', ls='dashed')  # label='$R_{Vir}$')
        a.text(1087, 0.83, "$R_{500}$", rotation=90, size=16)
        a.text(2147, 0.83, "$R_{vir}$", rotation=90, size=16)
        a.plot(rlogbin, zerohuit, color='grey')
        a.plot(rlogbin, un, color='black')
        a.set_xlim(200, 4500)
        a.set_ylim(0.75, 2.1)



        b.plot(rlogbin_cen, bsz, label="3D, no fit", color='orange', marker='.')
        b.plot(rlogbin_cen, bcorr, label="Non-thermal pressure correction, no fit",color='firebrick',marker='.',ls='dashed')
        b.plot(rlogbin_cen, bcorrfit2, label="Non-thermal pressure correction, fit ", color='firebrick',ls='dotted',marker='.')
        b.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        b.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)
        #b.plot(rlogbin, zerohuit, color='grey')
        #b.plot(rlogbin, un, color='black')
        b.set_xscale('log')
        b.legend(prop={'size': 16})
        b.set_ylabel("$(1-\mathrm{b})=rac{M_{\mathrm{HE}}}{M_{\mathrm{tot}}}$")
        b.set_xlabel("$R [\mathrm{kpc}]$")

        b.axvline(x=1087, color='grey')  # label='$R_{500}$')
        b.axvline(x=2147, color='grey', ls='dashed')  # label='$R_{Vir}$')
        b.text(1087, 0.83, "$R_{500}$", rotation=90, size=16)
        b.text(2147, 0.83, "$R_{vir}$", rotation=90, size=16)
        b.plot(rlogbin, zerohuit, color='grey')
        b.plot(rlogbin, un, color='black')
        b.set_xlim(200, 4500)
        b.set_ylim(0.75, 2.1)

        # plt.axvline(x=950, color='grey', ls='dotted')# label='$R_{bump}$')
        # plt.text(950, 1.1, "$R_{bump}$", rotation=90, size=16)


        plt.subplots_adjust(wspace=0, hspace=0)
        plt.show()
        print("end")
        sys.exit()
    #plt.plot(rlogbin_cen, bszl, label="Only cells with x<x_Virgo", color='orange')
    #plt.plot(rlogbin_cen, bszr, label="Only cells with x>x_Virgo", color='green')
    #plt.errorbar(rlogbin_cen, bsz, yerr=smsz, label="From P profile with weighted mean")
    #plt.fill_between(rlogbin_cen,bsz-smsz,bsz+smsz,color="blue",alpha=0.1)
    #plt.plot(rlogbin_cen, bx,label="From T and ne profiles with weighted mean",alpha=0.5,marker='+')
    #plt.fill_between(rlogbin_cen,bx-smx,bx+smx,color='orange',alpha=0.2)
    #plt.errorbar(rlogbin_cen, bx, yerr=smx, label="From T and ne profiles with weighted mean")

    #plt.plot(rlogbin, bxg, label='From T and $n_e$ profile with weighted mean and median respectively')
    #lt.plot(rlogbin, bszg / 1.136818, label='From P profile with median')
    #plt.plot(rlogbin, bx17, label='From T and $n_e$ profiles with weighted mean, low res',c="blue",alpha=0.8)

    #plt.plot(rlogbin, bxl, label='From T and $n_e$ profile with weighted mean, only cells with x<x_virgo')
    #plt.plot(rlogbin, bxr, label='From T and $n_e$ profile with weighted mean, only cells with x>x_virgo')
    #plt.plot(rlogbin, bx17o, label='From T and $n_e$ profiles with weighted mean, low res K Pa m-3 units',ls="dashed",alpha=0.8,c="orange")
    #plt.plot(rlogbin, bsz17/1.136818, label='From P profile with weighted mean, low res', c="purple", alpha=0.8)
    #plt.plot(rlogbin, bsz17o/1.136818, label='From P profile with weighted mean, low res K Pa m-3 units', ls="dashed",alpha=0.8, c="red")
    #plt.plot(rlogbin, bsz17, label='From P profile with weighted mean, low res')


    #m_dmsum = np.array([np.sum(mdm_c8[0:i]) for i in range(40)])
    #plt.plot(rlogbin,m_sum,label="total mass: DM+baryons")
    #plt.plot(rlogbin, m_sum17, label="total mass: DM+baryons low res")
    #plt.plot(rlogbin, m_dmsum, label="total mass: DM only")
    # plt.plot(rlogbin,bxc,label='From T and $n_e$ profiles, median for ne, mean for T')
    # plt.plot(rlogbin,bsz,label='From P profile with weighted mean')
    # plt.plot(rlogbin,bszc,label='From P profile  with median')
    # plt.plot(rlogbin,bxmed,label='From T and $n_e$ profile with median')
    # plt.plot(rlogbin,bszc,label='From P profile  with median')
    # plt.axvline(x=1087,color='grey',label='$R_{500}$')
    # plt.axvline(x=2024,color='grey',ls='dashed',label='$R_{Vir}$')

    # print("bx",bx)
    # print("bxc",bxc)

    # plt.fill_between(rlogbin,bxm,bxp,color='blue',alpha=0.5)
    # plt.fill_between(rlogbin,bszm,bszp,color='orange',alpha=0.5)
    #np.savetxt('bsz_virgo_3D_h.txt',bsz)

    fit_plot()


    plt.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
    plt.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)

    un=[1 for i in range(40)]
    zerohuit=[0.8 for i in range(40)]
    un=np.array(un)
    zerohuit=np.array(zerohuit)
    plt.plot(rlogbin,zerohuit,color='grey')
    plt.plot(rlogbin,un,color='black')
    plt.xscale('log')
    plt.xlabel('$\mathrm{R/R_{500}}$')
    plt.legend()
    plt.ylabel("$(1-\mathrm{b})=rac{M_{\mathrm{HE}}}{M_{\mathrm{tot}}}$")
    #plt.ylabel('$(1-b)=\frac{M_{HE}}{M_{tot}}$')

    #plt.axvline(x=1087, color='grey')# label='$R_{500}$')
    #plt.axvline(x=2147, color='grey', ls='dashed')# label='$R_{Vir}$')
    #plt.text(1087, 1.1, "$R_{500}$", rotation=90, size=16)
    #plt.text(2147, 1.1, "$R_{vir}$", rotation=90, size=16)

    plt.axvline(x=1, color='grey')  # label='$R_{500}$')
    plt.axvline(x=2147 / 1087, color='grey', ls='dashed')  # label='$R_{Vir}$')
    plt.text(1, 1.1, "$R_{500}$", rotation=90, size=16)
    plt.text(2147 / 1087, 1.1, "$R_{vir}$", rotation=90, size=16)
    plt.plot(rlogbin/r500, zerohuit, color='grey')
    plt.plot(rlogbin/r500, un, color='black')

    #plt.axvline(x=950, color='grey', ls='dotted')# label='$R_{bump}$')
    #plt.text(950, 1.1, "$R_{bump}$", rotation=90, size=16)
    plt.xlim(200/r500, 4500/r500)
    plt.ylim(0.75,2.1)
    #plt.ylim(0.3, 2.5)
    #plt.grid()
    plt.show()
    sys.exit()

#bias_low_res()

#print(Pstdlogmean_c8)

#Pstdlogmean_c8 /= np.sqrt(nlog_c8)
#nestdlogmean_c8 /= np.sqrt(nlog_c8)
#Tstdlogmean_c8 /= np.sqrt(nlog_c8)

#print("Pstdlogmean_c8",Pstdlogmean_c8)
#sys.exit()
#bias()

#sys.exit()

#Plogmed_c8_2=Plogmean_c8+np.log10(0.76/0.864)
def profiles_plots():
    # plt.plot(rlogbin,nelog4,'-',label="All cells")
    # plt.plot(rlogbin,nelog2,'.',label="Only T>10^7")
    # plt.plot(rlogbin,nelog,'.',label="Only T>10^5")
    # print(nestdlog)

    #plt.scatter(rlogbin, nelogmed, s=6, c='red')
    #plt.errorbar(rlogbin, nelogmed, yerr=nestdlogmed, ls='dotted', label='High Res, med, all cells', c='red',alpha=0.9)
    rlogbin60 = [10 ** ((20 + i) * 0.05) for i in range(60)]
    #plt.scatter(rlogbin60, nelogmed_c8, s=6)
    #plt.errorbar(rlogbin60, nelogmed_c8, yerr=nestdlogmed_c8, ls='dotted', label='High res, med, 1e8.5Msun gal cut', alpha=0.7)


    plt.scatter(rlogbin, Plogmed_c8, s=6)
    plt.errorbar(rlogbin, Plogmed_c8, yerr=Pstdlogmed_c8, ls='dotted', label='High res, med, 1e8.5Msun gal cut', alpha=0.7)
    plt.scatter(rlogbin, Plogmed_c8_t, s=6)
    plt.errorbar(rlogbin, Plogmed_c8_t, yerr=Pstdlogmed_c8_t, ls='dotted', label='High res, med, 1e8.5Msun gal cut, test',alpha=0.7)

    #plt.scatter(rlogbin_test, Plogmed_c8, s=6)
    #plt.errorbar(rlogbin_test, Plogmed_c8, yerr=Pstdlogmed_c8, ls='dotted', label='rlogbintest', alpha=0.7)
    #plt.scatter(rlogbin, Plogmed_c8l, s=6)
    #plt.errorbar(rlogbin, Plogmed_c8l, yerr=Pstdlogmed_c8l, ls='dotted', label='High res, med, 1e8.5Msun gal cut, only cells with x<x_virgo',alpha=0.7)
    #plt.scatter(rlogbin, Plogmed_c8r, s=6)
    #plt.errorbar(rlogbin, Plogmed_c8r, yerr=Pstdlogmed_c8r, ls='dotted', label='High res, med, 1e8.5Msun gal cut, only cells with x>x_virgo',alpha=0.7)

    #plt.scatter(rlogbin, Plogmed_cf, s=6)
    #plt.errorbar(rlogbin, Plogmed_cf, yerr=Pstdlogmed_cf, ls='dotted', label='High res, med, all galaxiees cut', alpha=0.7)


    #plt.plot(P_urban[:, 0], P_urban[:, 1]), '.',label='Observation datapoints from Urban et al.',c='black',ms=7)

    plt.scatter(P_urban[:, 0], P_urban[:, 1],s=7)
    plt.errorbar(P_urban[:, 0], P_urban[:, 1], xerr=[P_urban[:, 2], P_urban[:, 3]],label='Observation datapoints, Urban+ 2011')
    plt.scatter(P_planck[:, 0], P_planck[:, 1], s=7)
    plt.errorbar(P_planck[:, 0], P_planck[:, 1], xerr=[P_planck[:, 2], P_planck[:, 3]],label='Observation datapoints, Planck Collab. 2018')
    #plt.scatter(r_deg16,np.log10(P_planck_deg[:,1]),s=7)
    #plt.errorbar(r_deg16,np.log10(P_planck_deg[:,1]),xerr=[P_planck[:, 2], P_planck[:, 3]],label='Obs data from Planck (degree to kpc, d=16Mpc)')
    #plt.scatter(r_deg21, np.log10(P_planck_deg[:, 1]), s=7)
    #plt.errorbar(r_deg21, np.log10(P_planck_deg[:, 1]), xerr=[P_planck[:, 2], P_planck[:, 3]], label='Obs data from Planck (degree to kpc, d=21Mpc)')

    #plt.scatter(ne_urban[:, 0], ne_urban[:, 1], s=7)
    #plt.errorbar(ne_urban[:, 0], ne_urban[:, 1], xerr=[ne_urban[:, 2], ne_urban[:, 3]],label='Observation datapoints, Urban+ 2011')
    #plt.scatter(ne_planck[:, 0], ne_planck[:, 1], s=7)
    #plt.errorbar(ne_planck[:, 0], ne_planck[:, 1], xerr=[ne_planck[:, 2], ne_planck[:, 3]],label='Observation datapoints, Planck Collab. 2018')
    #plt.plot(ne_ghiz[:,0],np.log10(ne_ghiz[:,1]),'.',label="Observation datapoints, Ghizzardi 2004")


    # plt.scatter(rlogbin,Tlog3,s=6)
    # er2=plt.errorbar(rlogbin,Tlog3,ls='dashed',yerr=Tstdlog3,label='r_vir cut')
    # er2[-1][0].set_linestyle('-.')
    #UPP=np.array([np.log10(f.upp((10 ** ((20 + i) * 0.05))/770,3.5,2.5,1.33,3,0.8)*10**(-3.5371900826446278)) for i in range(60)])
    #UPPlog = np.array([f.upplog((10 ** ((20 + i) * 0.05)) / 770, 3.5, 2.5, 1.33, 3, 0.8)-3.5371900826446278 for i in range(60)])
    #print(UPP)
    #print(rlogbin)

    #plt.plot(rlogbinupp,UPP,label='A10')
    #plt.plot(rlogbinupp, UPPlog, label='A10')

    plt.axvline(x=1087, color='grey', label='$R_{500}$')
    plt.axvline(x=2024, color='grey', ls='dashed', label='$R_{Vir}$')
    plt.axvline(x=850, color='grey', ls='dotted', label='$R_{bump}$')
    plt.legend()


    plt.xscale('log')
    plt.xlabel('R (kpc)')
    #plt.ylabel('$log_{10}(T[keV])$')
    plt.ylabel('$log_{10}(P[keV/cm^3])$')
    #plt.ylabel('$log_{10}(n_e[cm^{-3}])$')
    plt.show()


#profiles_plots()

#sys.exit()

# P500=np.nansum(P3[0:108]*n3[0:108])/np.nansum(n3[0:108])
# print("P500",P500)
# n_e500=np.nansum(n_e[0:108]*n[0:108])/np.sum(n[0:108])
# print("P500",P500)
# print("n_e500",n_e500)
# Plog3-=np.log10(P500)

# f.press(Plog,n,1)


# print("r500")
# compar_bias(1085)
# print("r200")
# compar_bias(1705)
# print("rvir")
# compar_bias(2025)

# plt.plot(rbin,P,'.',label='raw',markersize=3)
# plt.plot(rbin,Ps,'.',label='smoothed',markersize=3,alpha=0.5)
# plt.plot(rbin,P,'.')

def gal_dist():
    gal = np.loadtxt('list_gal_251.dat_js_nocontam')
    # plt.plot(gal[2:,2],gal[2:,23]*(unit_l/3.08567758128E21),'.')
    plt.hist(np.log10(gal[:, 2]), bins=200, color='black')
    plt.axvline(x=9, label='Resolution limit', color='grey', ls='dashed')
    plt.legend()
    plt.show()

def three_d_plots():
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    ax1.scatter(rlogbin_cen, Plogmean_17, s=6, c='blue')
    ax1.errorbar(rlogbin_cen, Plogmean_17, yerr=Pstdlogmean_17, ls='dotted', label='$2048^3$', alpha=0.7,c='blue')  # 'Smallest grid size: 5.62kpc'
    ax1.scatter(rlogbin_cen, Plogmean_c8, s=6, c='orange')
    ax1.errorbar(rlogbin_cen, Plogmean_c8, yerr=Pstdlogmean_c8, ls='dotted', label='$8192^3$', alpha=0.7,c='orange')  # '''Smallest grid size: 0.35kpc
    ax1.set_xlabel("$R$[kpc]", size=16)
    ax1.set_ylabel("lo$\mathrm{g_{10}}$($P$[keV.c$\mathrm{m^{-3}}$])", size=16)
    ax1.set_xscale('log')
    ax1.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
    ax1.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)
    ax1.axvline(x=1087, color='grey')
    ax1.text(1087, -7, "$\mathrm{R_{500}}$", rotation=90, size=16)
    ax1.axvline(x=2147, color='grey', ls='dashed')
    ax1.text(2147, -7, "$\mathrm{R_{vir}}$", rotation=90, size=16)
    ax1.legend(prop={'size': 16})

    ax2.scatter(rlogbin_cen, nelogmean_17, s=6, c='blue')
    ax2.errorbar(rlogbin_cen, nelogmean_17, yerr=nestdlogmean_17, ls='dotted', label='$2048^3$', alpha=0.7,c='blue')  # 'Smallest grid size: 5.62kpc'
    ax2.scatter(rlogbin_cen, nelogmean_c8, s=6, c='orange')
    ax2.errorbar(rlogbin_cen, nelogmean_c8, yerr=nestdlogmean_c8, ls='dotted', label='$8192^3$', alpha=0.7,c='orange')  # '''Smallest grid size: 0.35kpc
    ax2.set_ylabel("lo$\mathrm{g_{10}}$($n_{\mathrm{e}}$[c$\mathrm{m^{-3}}$])", size=16)
    ax2.set_xlabel("$R$[kpc]", size=16)
    ax2.set_xscale('log')
    ax2.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
    ax2.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)
    ax2.axvline(x=1087, color='grey')
    ax2.text(1087, -6.2, "$\mathrm{R_{500}}$", rotation=90, size=16)
    ax2.axvline(x=2147, color='grey', ls='dashed')
    ax2.text(2147, -6.2, "$\mathrm{R_{vir}}$", rotation=90, size=16)
    ax2.legend(prop={'size': 16})

    ax3.scatter(rlogbin_cen, 10 ** Tlogmean_17, s=6, c='blue')
    ax3.errorbar(rlogbin_cen, 10 ** Tlogmean_17, yerr=10 ** Tlogmean_17 * Tstdlogmean_17, ls='dotted', label='$2048^3$',alpha=0.7, c='blue')
    ax3.scatter(rlogbin_cen, 10 ** Tlogmean_c8, s=6, c='orange')
    ax3.errorbar(rlogbin_cen, 10 ** Tlogmean_c8, yerr=10 ** Tlogmean_c8 * Tstdlogmean_c8, ls='dotted', label='$8192^3$',alpha=0.7, c='orange')
    ax3.set_ylabel("$T$[keV]", size=16)
    ax3.set_xlabel("$R$[kpc]", size=16)
    ax3.set_xscale('log')
    ax3.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
    ax3.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)
    ax3.axvline(x=1087, color='grey')
    ax3.text(1087, 0.2, "$\mathrm{R_{500}}$", rotation=90, size=16)
    ax3.axvline(x=2147, color='grey', ls='dashed')
    ax3.text(2147, 0.2, "$\mathrm{R_{vir}}$", rotation=90, size=16)
    ax3.legend(prop={'size': 16})

    ax4.scatter(rlogbin_cen, np.log10(k_l), s=6, c='blue')
    ax4.errorbar(rlogbin_cen, np.log10(k_l), yerr=k_l_std / k_l, ls='dotted', label='$2048^3$', alpha=0.7, c='blue')
    ax4.scatter(rlogbin_cen, np.log10(k_h), s=6, c='orange')
    ax4.errorbar(rlogbin_cen, np.log10(k_h), yerr=k_h_std / k_h, ls='dotted', label='$8192^3$', alpha=0.7, c='orange')
    ax4.set_ylabel('lo$\mathrm{g_{10}}$($K$[keV.c$\mathrm{m^{2}}$])')
    ax4.set_xlabel("$R$[kpc]", size=16)
    ax4.set_xscale('log')
    ax4.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
    ax4.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)
    ax4.axvline(x=1087, color='grey')
    ax4.text(1087, 2.6, "$\mathrm{R_{500}}$", rotation=90, size=16)
    ax4.axvline(x=2147, color='grey', ls='dashed')
    ax4.text(2147, 2.6, "$\mathrm{R_{vir}}$", rotation=90, size=16)
    ax4.legend(prop={'size': 16})

    fig.suptitle('3D profiles')

    plt.show()
    sys.exit()


def three_d_plots_2():

    #print('Plogmean_c8',Plogmean_c8)
    #print('nelogmean_c8',nelogmean_c8)
    #sys.exit()

    #print("rlogbincen",rlogbin_cen)
    #print("rlogbincen 1:16", rlogbin_cen[:16])



    #relat_dif_p=(Plogmean_c8[8:]-Plogmean_17[8:])/Plogmean_c8[8:]
    #print("relat dif p",np.mean(np.absolute(relat_dif_p)))

    #relat_dif_ne = (nelogmean_c8[8:] - nelogmean_17[8:]) / nelogmean_c8[8:]
    #print("relat dif ne", np.mean(np.absolute(relat_dif_ne)))

    #relat_dif_t = (10**Tlogmean_c8[:3] - 10**Tlogmean_19[:3]) / 10**Tlogmean_c8[:3]
    #relat_dif_t = (10 ** Tlogmean_c8[16:37] - 10 ** Tlogmean_19[16:37]) / 10 ** Tlogmean_c8[16:37]
    #relat_dif_t = (10**Tlogmean_c8 - 10**Tlogmean_17) / 10**Tlogmean_c8
    #print('mean relat diff',relat_dif_t)
    #print("mean abs relat dif T", np.mean(np.absolute(relat_dif_t)))

    #sys.exit()
    fig, ((ax1), (ax2), (ax3)) = plt.subplots(ncols=1,nrows=3)

    #ax1.scatter(rlogbin_cen, Plogmean_19, s=6, c='blue')
    #ax1.errorbar(rlogbin_cen, Plogmean_19, yerr=Pstdlogmean_19, ls='dotted', label='$2048^3$', alpha=0.7, c='blue')
    #ax1.scatter(rlogbin_cen, Plogmean_17, s=6, c='blue')
    #ax1.errorbar(rlogbin_cen, Plogmean_17, yerr=Pstdlogmean_17, ls='dotted', label='$2048^3$', alpha=0.7,c='blue')  # 'Smallest grid size: 5.62kpc'
    ax1.scatter(rlogbin_cen, Plogmean_c8, s=6, c='orange')
    ax1.errorbar(rlogbin_cen, Plogmean_c8, yerr=Pstdlogmean_c8, ls='dotted', label='$8192^3$', alpha=0.7,c='orange')  # '''Smallest grid size: 0.35kpc
    #ax1.set_xlabel("$R$[kpc]", size=16)
    ax1.set_ylabel("lo$\mathrm{g_{10}}$($P$[keV.c$\mathrm{m^{-3}}$])", size=16)
    ax1.set_xscale('log')
    ax1.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
    ax1.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)
    ax1.axvline(x=1087, color='grey')
    ax1.text(1087, -7, "$\mathrm{R_{500}}$", rotation=90, size=16)
    ax1.axvline(x=2147, color='grey', ls='dashed')
    ax1.text(2147, -7, "$\mathrm{R_{vir}}$", rotation=90, size=16)
    #ax1.legend(prop={'size': 16})
    ax1.set_xticks([])

    #ax2.scatter(rlogbin_cen, nelogmean_19, s=6, c='blue')
    #ax2.errorbar(rlogbin_cen, nelogmean_19, yerr=nestdlogmean_19, ls='dotted', label='$2048^3$', alpha=0.7, c='blue')
    #ax2.scatter(rlogbin_cen, nelogmean_17, s=6, c='blue')
    #ax2.errorbar(rlogbin_cen, nelogmean_17, yerr=nestdlogmean_17, ls='dotted', label='$2048^3$', alpha=0.7,c='blue')  # 'Smallest grid size: 5.62kpc'
    ax2.scatter(rlogbin_cen, nelogmean_c8, s=6, c='orange')
    ax2.errorbar(rlogbin_cen, nelogmean_c8, yerr=nestdlogmean_c8, ls='dotted', label='$8192^3$', alpha=0.7,c='orange')  # '''Smallest grid size: 0.35kpc
    ax2.set_ylabel("lo$\mathrm{g_{10}}$($n_{\mathrm{e}}$[c$\mathrm{m^{-3}}$])", size=16)
    #ax2.set_xlabel("$R$[kpc]", size=16)
    ax2.set_xscale('log')
    ax2.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
    ax2.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)
    ax2.axvline(x=1087, color='grey')
    #ax2.text(1087, -6.2, "$\mathrm{R_{500}}$", rotation=90, size=16)
    ax2.axvline(x=2147, color='grey', ls='dashed')
    #ax2.text(2147, -6.2, "$\mathrm{R_{vir}}$", rotation=90, size=16)
    #ax2.legend(prop={'size': 16})
    ax2.set_xticks([])

    #ax3.scatter(rlogbin_cen, 10**Tlogmean_19, s=6, c='blue')
    #ax3.errorbar(rlogbin_cen, 10**Tlogmean_19, yerr=10**Tlogmean_19 * Tstdlogmean_19, ls='dotted', label='$2048^3$', alpha=0.7, c='blue')
    #ax3.scatter(rlogbin_cen, 10 ** Tlogmean_17, s=6, c='blue')
    #ax3.errorbar(rlogbin_cen, 10 ** Tlogmean_17, yerr=10 ** Tlogmean_17 * Tstdlogmean_17, ls='dotted', label='$2048^3$',alpha=0.7, c='blue')
    ax3.scatter(rlogbin_cen, 10 ** Tlogmean_c8, s=6, c='orange')
    ax3.errorbar(rlogbin_cen, 10 ** Tlogmean_c8, yerr=10 ** Tlogmean_c8 * Tstdlogmean_c8, ls='dotted', label='$8192^3$',alpha=0.7, c='orange')

    #ax3.scatter(rlogbin_cen, Tlogmean_17, s=6, c='blue')
    #ax3.errorbar(rlogbin_cen, Tlogmean_17, yerr=Tstdlogmean_17, ls='dotted', label='$2048^3$',alpha=0.7, c='blue')
    #ax3.scatter(rlogbin_cen, Tlogmean_c8, s=6, c='orange')
    #ax3.errorbar(rlogbin_cen, Tlogmean_c8, yerr=Tstdlogmean_c8, ls='dotted', label='$8192^3$',alpha=0.7, c='orange')

    ax3.set_ylabel("$T$[keV]", size=16)
    ax3.set_xlabel("$R$[kpc]", size=16)
    ax3.set_xscale('log')
    ax3.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
    ax3.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)
    ax3.axvline(x=1087, color='grey')
    #ax3.text(1087, 0.2, "$\mathrm{R_{500}}$", rotation=90, size=16)
    ax3.axvline(x=2147, color='grey', ls='dashed')
    #ax3.text(2147, 0.2, "$\mathrm{R_{vir}}$", rotation=90, size=16)
    #ax3.legend(prop={'size': 16})

    #ax4.scatter(rlogbin_cen, np.log10(k_l), s=6, c='blue')
    #ax4.errorbar(rlogbin_cen, np.log10(k_l), yerr=k_l_std / k_l, ls='dotted', label='$2048^3$', alpha=0.7, c='blue')
    #ax4.scatter(rlogbin_cen, np.log10(k_h), s=6, c='orange')
    #ax4.errorbar(rlogbin_cen, np.log10(k_h), yerr=k_h_std / k_h, ls='dotted', label='$8192^3$', alpha=0.7, c='orange')
    #ax4.set_ylabel('lo$\mathrm{g_{10}}$($K$[keV.c$\mathrm{m^{2}}$])')
    #ax4.set_xlabel("$R$[kpc]", size=16)
    #ax4.set_xscale('log')
    #ax4.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
    #ax4.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)
    #ax4.axvline(x=1087, color='grey')
    #ax4.text(1087, 2.6, "$\mathrm{R_{500}}$", rotation=90, size=16)
    #ax4.axvline(x=2147, color='grey', ls='dashed')
    #ax4.text(2147, 2.6, "$\mathrm{R_{vir}}$", rotation=90, size=16)
    #ax4.legend(prop={'size': 16})

    plt.subplots_adjust(wspace=0, hspace=0)

    fig.suptitle('3D profiles')

    plt.show()
    sys.exit()

#three_d_plots_2()

def projection():
    #plt.scatter(rlogbin_cen, Plogmed_c8, s=6)
    #plt.errorbar(rlogbin_cen, Plogmed_c8, yerr=Pstdlogmed_c8, ls='dotted', label='3D high res, med, 1e8.5Msun gal cut',alpha=0.7)
    #plt.scatter(rlogbin_cen, nelogmed_c8, s=6)
    #plt.errorbar(rlogbin_cen, nelogmed_c8, yerr=nestdlogmed_c8, ls='dotted', label='3D high res, med, 1e8.5Msun gal cut',alpha=0.7)

    #plt.scatter(rlogbin_cen, 10**Tlogmean_c8, s=6,c='orange')  #T in Kev (normal scale)
    #plt.errorbar(rlogbin_cen, 10**Tlogmean_c8, yerr=Tstdlogmean_c8*np.log(10)*10**Tlogmean_c8, ls='dotted', label='3D, 8192^3',alpha=0.7, c='orange' )

    #plt.scatter(rlogbin_cen, Plogmean_17, s=6, c='blue')
    #plt.errorbar(rlogbin_cen, Plogmean_17, yerr=Pstdlogmean_17, ls='dotted', label='3D, $2048^3$', alpha=0.7,c='blue')  # 'Smallest grid size: 5.62kpc'
    #plt.scatter(rlogbin_cen, Plogmean_c8, s=6, c='orange')
    #plt.errorbar(rlogbin_cen, Plogmean_c8, yerr=Pstdlogmean_c8, ls='dotted', label='3D, $8192^3$', alpha=0.7,c='orange')  # '''Smallest grid size: 0.35kpc
    #plt.xlabel("R(kpc)", size=16)
    #plt.ylabel("$log_{10}(P[keV/cm^3])$", size=16)
    #plt.xscale('log')
    #plt.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
    #plt.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)
    #plt.axvline(x=1087, color='grey')
    #plt.text(1087, -6, "$R_{500}$", rotation=90, size=16)
    #plt.axvline(x=2147, color='grey', ls='dashed')
    #plt.text(2147, -6, "$R_{vir}$", rotation=90, size=16)
    #plt.legend(prop={'size': 16})
    #f.los_distrib("virgo_xyz_hydro_l19_MW_los.dat")
    #f.los_distrib("virgo_xyz_hydro_l19_gal_clean_m1e8.5.dat")

    #f.projection("./maps/low_res/map_low_17_fil_P_los.bin", "P", 0, nelogmean_c8, nestdlogmean_c8, 0, 0, 1)

    #f.projection("./maps/high_res/P_maps/map_high_19_fil_P_los.bin", "P", 0, Plogmean_c8, Pstdlogmean_c8, 1, 0, 0)

    #f.projection("./maps/high_res/map_high_21_fil_P_los_nline.bin", "P", 0, Plogmean_c8, Pstdlogmean_c8, 0, 1,0)

    #f.projection("./toy_models/sphere_toy_map_column_density_1e5-1e-5.npy", "ne", 0, Plogmean_c8, Pstdlogmean_c8, 0, 0, 1)

    #plt.show()
    #sys.exit()



    #plt.show()
    #sys.exit()

    #f.projection("./maps/high_res/P_maps/map_10000px_fily_P_los_cleanm1e8.5+T1e5_highr.txt", "P",10000, Plogmean_c8, Pstdlogmean_c8, 0, 1, 0)
    #f.projection("./maps/map_z_P_rescale_galrmv_low_test.txt", "P", 0, Plogmean_c8,Pstdlogmean_c8, 0, 0, 1)

    f.projection("./maps/low_res/map_low_19_fil_ne_los.bin", "ne",0,nelogmean_c8,nestdlogmean_c8, 1, 0)
    #f.projection("./maps/high_res/ne_maps/map_4000px_z_ne_los_cleanm1e8.5+T1e5_highr.txt","ne",4000, nelogmean_c8, nestdlogmean_c8, 0, 0, 1)

    #f.projection("./maps/low_res/map_2000px_z_P_los_l14_lowr.txt", "P", 2000, Plogmean_17, Pstdlogmean_17, 0,1, 0)
#plt.scatter(rlogbin_cen, Plogmean_c8, s=6, c='orange')
#plt.errorbar(rlogbin_cen, Plogmean_c8, yerr=Pstdlogmean_c8, ls='dotted', label='3D, 8192^3', alpha=0.7, c='orange')
#f.toy_model()
#sys.exit()
#f.gianfagna_overplot()
#sys.exit()
projection()

sys.exit()

# f.deprojection("map_2000px_P_Tsup7.txt","z")

def multiplot():
    n = 1
    for i in range(0, 8):
        print(i)
        plt.subplot(2, 4, n)
        m = 0
        # for i in range(10):
        # f.temp_prof_gal(i)
        f.local_env(i)
        # m+=1
        # titre="gal "+str(i)#str(10*(n-1))+"to "+str(10*n)
        # plt.title(titre)
        # plt.legend()
        n += 1

    plt.show()

# bias()
# deproj()

#g = np.loadtxt('/data/cluster/tlebeau/virgo/HighRes/list_gal_251.dat_js_nocontam')
#print("g")
#plt.hist(np.log10(g[:,2]),bins=200, color='black')
#plt.xscale('log')
#plt.show()
# print((np.max(gals[:,3])-np.min(gals[:,3]))*737.441)
# print((np.max(gals[:,4])-np.min(gals[:,4]))*737.441)
# print((np.max(gals[:,5])-np.min(gals[:,5]))*737.441)

xv = 0.48461068  # Virgo
yv = 0.50809848
zv = 0.49687076

def plot_jenny():
    plt.plot(g[:, 19][np.logical_and(g[:, 18] > 0, np.sqrt(
        (xv - g[:, 3]) ** 2 + (yv - g[:, 4]) ** 2 + (zv - g[:, 5]) ** 2) < 0.016272488)], g[:, 2][
                 np.logical_and(g[:, 18] > 0,
                                np.sqrt(
                                    (xv - g[:, 3]) ** 2 + (yv - g[:, 4]) ** 2 + (zv - g[:, 5]) ** 2) < 0.016272488)],
             '.', color="black", label="Satellite Galaxies", alpha=0.8)

    plt.plot(g[:, 19][np.logical_and(g[:, 18] == 1, np.sqrt(
        (xv - g[:, 3]) ** 2 + (yv - g[:, 4]) ** 2 + (zv - g[:, 5]) ** 2) < 0.016272488)], g[:, 2][
                 np.logical_and(g[:, 18] == 1,
                                np.sqrt(
                                    (xv - g[:, 3]) ** 2 + (yv - g[:, 4]) ** 2 + (zv - g[:, 5]) ** 2) < 0.016272488)],
             '.', color="red", label="Central Galaxies", alpha=0.8)

    plt.axhline(y=1e9)
    plt.axvline(x=2e10)


#p = np.loadtxt("/data/cluster/tlebeau/P_projected_profile_along_z")
#print(p)
#plt.plot(p,'-')
# plt.legend()
#plt.show()
def slow_virgo_prof():

    #data_p_3000 = np.loadtxt("./slow_virgo_prof/Virgo_pressure_3000_log.dat.txt")
    # print(data_p_3000[:,0])
    #p_mean_3000 = data_p_3000[:, 0]
    #print(len(data_p_3000))
    #r_low_3000 = data_p_3000[:, 1]
    #r_up_3000 = data_p_3000[:, 2]
    #rbin_3000 = r_low_3000 + 0.5 * (r_up_3000 - r_low_3000)
    #p_median_3000 = data_p_3000[:, 8]
    #p_median_16_3000 = data_p_3000[:, 6]
    #p_median_84_3000 = data_p_3000[:, 7]

    rlogbin_cen = np.array([10 ** ((35.5 + i) * 0.05) for i in range(40)])

    def p_plot():
        #data_p_1500 = np.loadtxt("./slow_virgo_prof/Virgo_pressure_1500_log.dat.txt")
        data_p_1500 = np.loadtxt("./slow_virgo_prof/pressure_1500_virgo_cutout.txt")
        p_mean_1500 = data_p_1500[:, 0]
        print(len(data_p_1500))
        r_low_1500 = data_p_1500[:, 1]
        r_up_1500 = data_p_1500[:, 2]
        rbin_1500 = r_low_1500 + 0.5 * (r_up_1500 - r_low_1500)
        p_median_1500 = data_p_1500[:, 8]
        p_median_16_1500 = data_p_1500[:, 6]
        p_median_84_1500 = data_p_1500[:, 7]

        plt.scatter(rlogbin_cen, Plogmean_17, s=6, c='blue')
        plt.errorbar(rlogbin_cen, Plogmean_17, yerr=Pstdlogmean_17, ls='dotted', label='CLONE, 3D $2048^3$', alpha=0.7,c='blue')  # 'Smallest grid size: 5.62kpc'

        plt.scatter(rlogbin_cen, Plogmean_c8, s=6, c='orange')
        plt.errorbar(rlogbin_cen, Plogmean_c8, yerr=Pstdlogmean_c8, ls='dotted', label='CLONE, 3D $8192^3$', alpha=0.7,c='orange')  # '''Smallest grid size: 0.35kpc

        plt.scatter(rbin_1500, np.log10(p_mean_1500),c='forestgreen',marker='.',ls='dotted')
        plt.errorbar(rbin_1500, np.log10(p_mean_1500), label="SLOW, 3D $1576^3$", c='tomato', marker='.', ls='dotted')

        plt.ylabel("$log_{10}(P[keV/cm^3])$", size=16)

    def ne_plot():
        data_ne_1500 = np.loadtxt("./slow_virgo_prof/Virgo_edensity_1500_log.dat.txt")
        ne_mean_1500 = data_ne_1500[:, 0]
        print(len(data_ne_1500))
        r_low_1500_ne = data_ne_1500[:, 1]
        r_up_1500_ne = data_ne_1500[:, 2]
        rbin_1500_ne = r_low_1500_ne + 0.5 * (r_up_1500_ne - r_low_1500_ne)
        ne_median_1500 = data_ne_1500[:, 8]
        ne_median_16_1500 = data_ne_1500[:, 6]
        ne_median_84_1500 = data_ne_1500[:, 7]

        plt.scatter(rlogbin_cen, nelogmean_17, s=6, c='blue')
        plt.errorbar(rlogbin_cen, nelogmean_17, yerr=nestdlogmean_17, ls='dotted', label='CLONE, 3D $2048^3$', alpha=0.7, c='blue')  # 'Smallest grid size: 5.62kpc'

        plt.scatter(rlogbin_cen, nelogmean_c8, s=6, c='orange')
        plt.errorbar(rlogbin_cen, nelogmean_c8, yerr=nestdlogmean_c8, ls='dotted', label='CLONE, 3D $8192^3$', alpha=0.7, c='orange')  # '''Smallest grid size: 0.35kpc

        plt.scatter(rbin_1500_ne, np.log10(ne_mean_1500), c='forestgreen', marker='.', ls='dotted')
        plt.errorbar(rbin_1500_ne, np.log10(ne_mean_1500), label="SLOW, 3D $1576^3$", c='tomato', marker='.', ls='dotted')

        plt.ylabel("lo$\mathrm{g_{10}}$($n_{\mathrm{e}}$[c$\mathrm{m^{-3}}$])")

    def t_plot():
        data_t_1500 = np.loadtxt("./slow_virgo_prof/Virgo_temp_1500_log.dat.txt")
        t_mean_1500 = data_t_1500[:, 0]
        print(len(data_t_1500))
        r_low_1500_t = data_t_1500[:, 1]
        r_up_1500_t = data_t_1500[:, 2]
        rbin_1500_t = r_low_1500_t + 0.5 * (r_up_1500_t - r_low_1500_t)
        t_median_1500 = data_t_1500[:, 8]
        t_median_16_1500 = data_t_1500[:, 6]
        t_median_84_1500 = data_t_1500[:, 7]

        plt.scatter(rlogbin_cen, 10 ** Tlogmean_17, s=6, c='blue')
        plt.errorbar(rlogbin_cen, 10 ** Tlogmean_17, yerr=10 ** Tlogmean_17 * Tstdlogmean_17, ls='dotted', label='$2048^3$', alpha=0.7, c='blue')
        plt.scatter(rlogbin_cen, 10 ** Tlogmean_c8, s=6, c='orange')
        plt.errorbar(rlogbin_cen, 10 ** Tlogmean_c8, yerr=10 ** Tlogmean_c8 * Tstdlogmean_c8, ls='dotted', label='$8192^3$', alpha=0.7, c='orange')

        plt.scatter(rbin_1500_t, t_mean_1500, c='forestgreen', marker='.', ls='dotted')
        plt.errorbar(rbin_1500_t, t_mean_1500, label="SLOW, 3D $1576^3$", c='tomato', marker='.',ls='dotted')

        plt.ylabel("$T$[keV]")


    p_plot()

    #plt.plot(rbin,p_median)
    #plt.plot(rbin,p_mean)
    plt.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
    plt.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)
    plt.axvline(x=1087, color='grey')
    plt.text(1087, -6, "$R_{500}$", rotation=90, size=16)
    plt.axvline(x=2147, color='grey', ls='dashed')
    plt.text(2147, -6, "$R_{vir}$", rotation=90, size=16)
    plt.legend(prop={'size': 16})
    plt.xscale('log')
    plt.xlabel("R(kpc)", size=16)
    plt.show()

def find_cluster_Cen_Los():
    d=np.loadtxt("./HighRes/list_halo_251.dat_js_nocontam_high_res")
    cond=np.logical_and(d[:,3]>0.494,np.logical_and(d[:,3]<0.496,np.logical_and(d[:,4]>0.502,np.logical_and(d[:,4]<0.504,np.logical_and(d[:,5]>0.498,d[:,5]<0.5)))))
    print(len(d[:,0][cond]))
    print(d[:,2:8][cond])


#slow_virgo_prof()

#find_cluster_Cen_Los()


##### Project 2 on Virgo ######

def load_sectors_data(file):

    rad_data = np.load(file)
    p_rad_log = rad_data[0, :]
    err_p = rad_data[1, :]
    ne_rad_log = rad_data[2, :]
    err_ne = rad_data[3, :]
    t_rad_log = rad_data[4, :]
    err_t = rad_data[5, :]
    k_rad_log = rad_data[6, :]
    err_k = rad_data[7, :]
    m_rad_log = rad_data[8, :]
    n_log = rad_data[9, :]

    return p_rad_log, err_p, ne_rad_log, err_ne, t_rad_log, err_t, k_rad_log, err_k, m_rad_log, n_log

def rad_press_prof_sectors():

    def dm_mass_prof(filedm):
        #filedm = "./virgo_xyz_files/virgo_xyz_dm_high_res.dat"
        d = FortranFile(filedm, 'r')

        ncell = d.read_ints()

        xdm, ydm, zdm, mdm, vx, vy, vz = ftp.f90_to_py.read_dm(ncell, filedm)

        rdm = np.sqrt(xdm ** 2 + ydm ** 2 + zdm ** 2)

        rlog = np.log10(rdm)
        rlog = np.array(rlog)
        # print(rlog)

        print('rayon part dm fini')

        # nbin=40
        nbin = 77

        m_rad = np.zeros(nbin)
        n_rad = np.zeros(nbin)

        # rlogbin = np.array([10 ** ((17.5 + i) * 0.1) for i in range(nbin + 1)])
        #rlogbin = np.array([10 ** ((35 + i) * 0.05) for i in range(40)])

        #nbin = 77
        rlogbin_cen = np.array([10 ** ((92.5 + i) * 0.025) for i in range(nbin)])

        m_inside_200kpc= np.sum(mdm[rdm<200])

        for i in range(nbin):
            # cond = np.logical_and(rlog > (i + 35) * 0.1, rlog < (i + 36) * 0.1)
            cond = np.logical_and(rlog > (i + 92) * 0.025, rlog < (i + 93) * 0.025)
            m_rad[i] = np.sum(mdm[cond])
            n_rad[i] = len(mdm[cond])

        return m_rad, n_rad,m_inside_200kpc

    #m_dm,n_dm,m_inside_200kpc = dm_mass_prof("./virgo_xyz_files/virgo_xyz_dm_high_res.dat")


    #np.save("m_dm_77_bins.npy",m_dm)
    #np.save("m_dm_inside_200kpc_77_bins.npy",m_inside_200kpc)

    m_dm = np.load("m_dm_77_bins.npy")
    m_dm_inside_200kpc = np.load("m_dm_inside_200kpc_77_bins.npy")
    m_bar_inside_200kpc = np.load("m_bar_inside_200kpc_77_bins.npy")
    print("mbar(r<200kpc)=", "{:.2e}".format(m_bar_inside_200kpc))
    print("mdm(r<200kpc)=","{:.2e}".format(m_dm_inside_200kpc))
    #Sectors definition for simplicity of variables :
    # 1 : x<x_cen, y<y_cen, z<z_cen
    # 2 : x<x_cen, y<y_cen, z>z_cen
    # 3 : x<x_cen, y>y_cen, z<z_cen
    # 4 : x<x_cen, y>y_cen, z>z_cen
    # 5 : x>x_cen, y<y_cen, z<z_cen
    # 6 : x>x_cen, y<y_cen, z>z_cen
    # 7 : x>x_cen, y>y_cen, z<z_cen
    # 8 : x>x_cen, y>y_cen, z>z_cen


    #p1, err_p1, ne1, err_ne1, t1, err_t1, k1, err_k1, m1, n1 = load_sectors_data("./sectors_study/lograd_profs_h_xi_yi_zi.npy")
    #p2, err_p2, ne2, err_ne2, t2, err_t2, k2, err_k2, m2, n2 = load_sectors_data("./sectors_study/lograd_profs_h_xi_yi_zs.npy")
    #p3, err_p3, ne3, err_ne3, t3, err_t3, k3, err_k3, m3, n3 = load_sectors_data("./sectors_study/lograd_profs_h_xi_ys_zi.npy")
    #p4, err_p4, ne4, err_ne4, t4, err_t4, k4, err_k4, m4, n4 = load_sectors_data("./sectors_study/lograd_profs_h_xi_ys_zs.npy")
    #p5, err_p5, ne5, err_ne5, t5, err_t5, k5, err_k5, m5, n5 = load_sectors_data("./sectors_study/lograd_profs_h_xs_yi_zi.npy")
    #p6, err_p6, ne6, err_ne6, t6, err_t6, k6, err_k6, m6, n6 = load_sectors_data("./sectors_study/lograd_profs_h_xs_yi_zs.npy")
    #p7, err_p7, ne7, err_ne7, t7, err_t7, k7, err_k7, m7, n7 = load_sectors_data("./sectors_study/lograd_profs_h_xs_ys_zi.npy")
    #p8, err_p8, ne8, err_ne8, t8, err_t8, k8, err_k8, m8, n8 = load_sectors_data("./sectors_study/lograd_profs_h_xs_ys_zs.npy")
    #p_all, err_p_all, ne_all, err_ne_all, t_all, err_t_all, k_all, err_k_all, m_all, n_all = load_sectors_data("./sectors_study/lograd_profs_h_all.npy")


    ##extension of range for splashback radius study
    #p1, err_p1, ne1, err_ne1, t1, err_t1, k1, err_k1, m1, n1 = load_sectors_data("./sectors_study/splashback/lograd_profs_h_sect1_77.npy")
    #p5, err_p5, ne5, err_ne5, t5, err_t5, k5, err_k5, m5, n5 = load_sectors_data("./sectors_study/splashback/lograd_profs_h_sect5_77.npy")
    #p8, err_p8, ne8, err_ne8, t8, err_t8, k8, err_k8, m8, n8 = load_sectors_data("./sectors_study/splashback/lograd_profs_h_sect8_77.npy")
    #p_all, err_p_all, ne_all, err_ne_all, t_all, err_t_all, k_all, err_k_all, m_all, n_all = load_sectors_data("./sectors_study/splashback/lograd_profs_h_all_77.npy")

    ##Method modif : log after mean:
    p1, err_p1, ne1, err_ne1, t1, err_t1, k1, err_k1, m1, n1 = load_sectors_data("./sectors_study/splashback/lograd_profs_h_sect1_77_mean_before_log.npy")
    p5, err_p5, ne5, err_ne5, t5, err_t5, k5, err_k5, m5, n5 = load_sectors_data("./sectors_study/splashback/lograd_profs_h_sect5_77_mean_before_log.npy")
    p8, err_p8, ne8, err_ne8, t8, err_t8, k8, err_k8, m8, n8 = load_sectors_data("./sectors_study/splashback/lograd_profs_h_sect8_77_mean_before_log.npy")
    p_all, err_p_all, ne_all, err_ne_all, t_all, err_t_all, k_all, err_k_all, m_all, n_all = load_sectors_data("./sectors_study/splashback/lograd_profs_h_all_77_mean_before_log.npy")

    p1_med, err_p1_med, ne1_med, err_ne1_med, t1_med, err_t1_med, k1_med, err_k1_med, m1_med, n1_med = load_sectors_data("./sectors_study/splashback/lograd_profs_h_sect1_77_med_before_log.npy")
    p5_med, err_p5_med, ne5_med, err_ne5_med, t5_med, err_t5_med, k5_med, err_k5_med, m5_med, n5_med = load_sectors_data("./sectors_study/splashback/lograd_profs_h_sect5_77_med_before_log.npy")
    p8_med, err_p8_med, ne8_med, err_ne8_med, t8_med, err_t8_med, k8_med, err_k8_med, m8_med, n8_med = load_sectors_data("./sectors_study/splashback/lograd_profs_h_sect8_77_med_before_log.npy")
    p_all_med, err_p_all_med, ne_all_med, err_ne_all_med, t_all_med, err_t_all_med, k_all_med, err_k_all_med, m_all_med, n_all_med = load_sectors_data("./sectors_study/splashback/lograd_profs_h_all_77_med_before_log.npy")

    m_cumul_sum = np.cumsum(m_all)+np.cumsum(m_dm)+m_dm_inside_200kpc+m_bar_inside_200kpc #cumulative sum of mass
    #print("m_cumul_all_bar", np.cumsum(m_all))
    #print("m_cumul_dm", np.cumsum(m_dm))
    #print("m_cumul_sum",m_cumul_sum)
    nbin=77

    rlogbin_cen = np.array([10 ** ((92.5 + i) * 0.025) for i in range(nbin)])

    #n=90
    #rlogbin_cen_bar = np.array([10 ** ((92.5 + i) * 0.025) for i in range(n)])

    def plot_pressure():

        plt.scatter(rlogbin_cen, p_all, s=6, c='black')
        plt.errorbar(rlogbin_cen, p_all, yerr=err_p_all, ls='dashed', label='Full box', alpha=0.7, c='black')

        plt.scatter(rlogbin_cen, p1, s=6,c='green')#, c='blue')
        plt.errorbar(rlogbin_cen, p1, yerr=err_p1, ls='dotted', label='Isolated sector', alpha=0.7,c='green')#, c='blue')

        #plt.scatter(rlogbin_cen, p2, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, p2, yerr=err_p2, ls='dotted', label='sector 2', alpha=0.7)  # , c='blue')

        #plt.scatter(rlogbin_cen, p3, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, p3, yerr=err_p3, ls='dotted', label='sector 3', alpha=0.7)  # , c='blue')

        #plt.scatter(rlogbin_cen, p4, s=6, c='red')
        #plt.errorbar(rlogbin_cen, p4, yerr=err_p4, ls='dotted', label='sector 4', alpha=0.7, c='red')

        plt.scatter(rlogbin_cen, p5, s=6, c='red')
        plt.errorbar(rlogbin_cen, p5, yerr=err_p5, ls='dotted', label='Matter infall sector', alpha=0.7, c='red')

        #plt.scatter(rlogbin_cen, p6, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, p6, yerr=err_p6, ls='dotted', label='sector 6', alpha=0.7)  # , c='blue')

        #plt.scatter(rlogbin_cen, p7, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, p7, yerr=err_p7, ls='dotted', label='sector 7', alpha=0.7)  # , c='blue')

        plt.scatter(rlogbin_cen, p8, s=6, c='blue')
        plt.errorbar(rlogbin_cen, p8, yerr=err_p8, ls='dotted', label='In-between sector', alpha=0.7, c='blue')


        plt.ylabel("lo$\mathrm{g_{10}}$($P$[keV.c$\mathrm{m^{-3}}$])", size=16)
        plt.xscale('log')
        plt.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        plt.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)
        plt.axvline(x=1087, color='grey')
        plt.text(1087, -7, "$\mathrm{R_{500}}$", rotation=90, size=16)
        plt.axvline(x=2147, color='grey', ls='dashed')
        plt.text(2147, -7, "$\mathrm{R_{vir}}$", rotation=90, size=16)

        #plt.axvline(x=350, color='grey', ls='dotted')
        #plt.text(350, -7, "Shock front (AGN)", rotation=90, size=16)

        #plt.axvline(x=850, color='grey', ls='dotted')
        #plt.text(850, -7, "Shock front (Matter infall)", rotation=90, size=16)

        plt.legend(prop={'size': 16})
        plt.xlabel("$R$[kpc]", size=16)
        #plt.set_xticks([])
        plt.xlim(50,4000)

        plt.show()

    def plot_temperature():

        plt.scatter(rlogbin_cen, 10 ** t_all, s=6, c='black')
        plt.errorbar(rlogbin_cen, 10 ** t_all, yerr=10 ** t_all * err_t_all, ls='dashed', label='Full box', alpha=0.7, c='black')

        plt.scatter(rlogbin_cen, 10**t1, s=6,c='green')  # , c='blue')
        plt.errorbar(rlogbin_cen, 10**t1, yerr=10**t1*err_t1, ls='dotted', label='Isolated sector', alpha=0.7, c='green')  # , c='blue')

        #plt.scatter(rlogbin_cen, 10**t2, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, 10**t2, yerr=10**t2*err_t2, ls='dotted', label='sector 2', alpha=0.7)  # , c='blue')

        #plt.scatter(rlogbin_cen, 10**t3, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, 10**t3, yerr=10**t3*err_t3, ls='dotted', label='sector 3', alpha=0.7)  # , c='blue')

        #plt.scatter(rlogbin_cen, 10**t4, s=6, c='red')
        #plt.errorbar(rlogbin_cen, 10**t4, yerr=10**t4*err_t4, ls='dotted', label='sector 4', alpha=0.7, c='red')

        plt.scatter(rlogbin_cen, 10**t5, s=6, c='red')
        plt.errorbar(rlogbin_cen, 10**t5, yerr=10**t5*err_t5, ls='dotted', label='Matter infall sector', alpha=0.7, c='red')

        #plt.scatter(rlogbin_cen, 10**t6, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, 10**t6, yerr=10**t6*err_t6, ls='dotted', label='sector 6', alpha=0.7)  # , c='blue')

        #plt.scatter(rlogbin_cen, 10**t7, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, 10**t7, yerr=10**t7*err_t7, ls='dotted', label='sector 7', alpha=0.7)  # , c='blue')

        plt.scatter(rlogbin_cen, 10**t8, s=6, c='blue')
        plt.errorbar(rlogbin_cen, 10**t8, yerr=10**t8*err_t8, ls='dotted', label='In-between sector', alpha=0.7, c='blue')

        plt.ylabel("$T$[keV]", size=16)
        plt.xscale('log')
        plt.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        plt.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,
                        labelleft=True)
        plt.axvline(x=1087, color='grey')
        plt.text(1087, 0, "$\mathrm{R_{500}}$", rotation=90, size=16)
        plt.axvline(x=2147, color='grey', ls='dashed')
        plt.text(2147, 0, "$\mathrm{R_{vir}}$", rotation=90, size=16)

        #plt.axvline(x=350, color='grey', ls='dotted')
        #plt.text(350, 0, "AGN feedback", rotation=90, size=16)

        plt.axvline(x=850, color='grey', ls='dotted')
        plt.text(850, 0, "Matter infall", rotation=90, size=16)

        plt.legend(prop={'size': 16})
        plt.xlabel("$R$[kpc]", size=16)
        # plt.set_xticks([])

        plt.show()

    def plot_density():

        plt.scatter(rlogbin_cen, ne_all, s=6, c='black')
        plt.errorbar(rlogbin_cen, ne_all, yerr=err_ne_all, ls='dashed', label='All sectors', alpha=0.7, c='black')

        plt.scatter(rlogbin_cen, ne1, s=6, c='green')  # , c='blue')
        plt.errorbar(rlogbin_cen, ne1, yerr=err_ne1, ls='dotted', label='Matter infall sector', alpha=0.7,c='green')  # , c='blue')

        #plt.scatter(rlogbin_cen, ne2, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, ne2, yerr=err_ne2, ls='dotted', label='sector 2', alpha=0.7)  # , c='blue')

        #plt.scatter(rlogbin_cen, ne3, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, ne3, yerr=err_ne3, ls='dotted', label='sector 3', alpha=0.7)  # , c='blue')

        #plt.scatter(rlogbin_cen, ne4, s=6 , c='red')
        #plt.errorbar(rlogbin_cen, ne4, yerr=err_ne4, ls='dotted', label='sector 4', alpha=0.7 , c='red')

        plt.scatter(rlogbin_cen, ne5, s=6, c='red')
        plt.errorbar(rlogbin_cen, ne5, yerr=err_ne5, ls='dotted', label='Isolated sector', alpha=0.7, c='red')

        #plt.scatter(rlogbin_cen, ne6, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, ne6, yerr=err_ne6, ls='dotted', label='sector 6', alpha=0.7)  # , c='blue')

        #plt.scatter(rlogbin_cen, ne7, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, ne7, yerr=err_ne7, ls='dotted', label='sector 7', alpha=0.7)  # , c='blue')

        plt.scatter(rlogbin_cen, ne8, s=6, c='blue')
        plt.errorbar(rlogbin_cen, ne8, yerr=err_ne8, ls='dotted', label='In-between sector', alpha=0.7, c='blue')


        plt.ylabel("lo$\mathrm{g_{10}}$($n_{\mathrm{e}}$[c$\mathrm{m^{-3}}$])", size=16)
        plt.xscale('log')
        plt.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        plt.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,
                        labelleft=True)
        plt.axvline(x=1087, color='grey')
        plt.text(1087, -6, "$\mathrm{R_{500}}$", rotation=90, size=16)
        plt.axvline(x=2147, color='grey', ls='dashed')
        plt.text(2147, -6, "$\mathrm{R_{vir}}$", rotation=90, size=16)

        #plt.axvline(x=350, color='grey', ls='dotted')
        #plt.text(350, -6, "Shock front (AGN)", rotation=90, size=16)

        #plt.axvline(x=850, color='grey', ls='dotted')
        #plt.text(850, -6, "Shock front (Matter infall)", rotation=90, size=16)

        plt.legend(prop={'size': 16})
        plt.xlabel("$R$[kpc]", size=16)
        # plt.set_xticks([])

        plt.show()

    def plot_entropy():
        plt.scatter(rlogbin_cen, k_all, s=6, c='black')
        plt.errorbar(rlogbin_cen, k_all, yerr=err_k_all, ls='dotted', label='All sectors', alpha=0.7, c='black')

        #plt.scatter(rlogbin_cen, k1, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, k1, yerr=err_k1, ls='dotted', label='sector 1', alpha=0.7)  # , c='blue')

        #plt.scatter(rlogbin_cen, k2, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, k2, yerr=err_k2, ls='dotted', label='sector 2', alpha=0.7)  # , c='blue')

        #plt.scatter(rlogbin_cen, k3, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, k3, yerr=err_k3, ls='dotted', label='sector 3', alpha=0.7)  # , c='blue')

        plt.scatter(rlogbin_cen, k4, s=6, c='red')
        plt.errorbar(rlogbin_cen, k4, yerr=err_k4, ls='dotted', label='sector 4', alpha=0.7, c='red')

        plt.scatter(rlogbin_cen, k5, s=6, c='purple')
        plt.errorbar(rlogbin_cen, k5, yerr=err_k5, ls='dotted', label='sector 5', alpha=0.7, c='purple')

        #plt.scatter(rlogbin_cen, k6, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, k6, yerr=err_k6, ls='dotted', label='sector 6', alpha=0.7)  # , c='blue')

        #plt.scatter(rlogbin_cen, k7, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, k7, yerr=err_k7, ls='dotted', label='sector 7', alpha=0.7)  # , c='blue')

        #plt.scatter(rlogbin_cen, k8, s=6)  # , c='blue')
        #plt.errorbar(rlogbin_cen, k8, yerr=err_k8, ls='dotted', label='sector 8', alpha=0.7)  # , c='blue')

        #plt.ylabel("lo$\mathrm{g_{10}}$($k$)", size=16)
        plt.ylabel('lo$\mathrm{g_{10}}$($K$[keV.c$\mathrm{m^{2}}$])',size=16)
        plt.xscale('log')
        plt.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        plt.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,
                        labelleft=True)
        plt.axvline(x=1087, color='grey')
        plt.text(1087, 4.5, "$\mathrm{R_{500}}$", rotation=90, size=16)
        plt.axvline(x=2147, color='grey', ls='dashed')
        plt.text(2147, 4.5, "$\mathrm{R_{vir}}$", rotation=90, size=16)

        plt.axvline(x=350, color='grey', ls='dotted')
        plt.text(350, 3.8, "Shock front (AGN)", rotation=90, size=16)

        plt.axvline(x=850, color='grey', ls='dotted')
        plt.text(850, 3.8, "Shock front (Matter infall)", rotation=90, size=16)

        plt.legend(prop={'size': 16})
        plt.xlabel("$R$[kpc]", size=16)
        plt.ylim(2.3,5.5)
        # plt.set_xticks([])

        plt.show()

    def plot_bias():
        m_cumul_dm = np.load("m_cumul_dm.npy")
        m_cumul_ba = np.load("m_cumul_ba.npy")
        m_sum = m_cumul_dm + m_cumul_ba

        alpha = np.load('alpha.npy') ##not good data, just to make the following function work, need to add turbulent pressure for the sectors

        #bx_all, bsz_all, smx_all, smsz_all, bcorr_all = f.ratios(p_all, t_all, ne_all, m_sum, 1, err_p_all, err_ne_all, err_t_all, alpha)
        bx1, bsz1, smx1, smsz1, bcorr1 = f.ratios(p1, t1, ne1, m_sum, 1, err_p1,err_ne1, err_t1, alpha)
        bx2, bsz2, smx2, smsz2, bcorr2 = f.ratios(p2, t2, ne2, m_sum, 1, err_p2, err_ne2, err_t2, alpha)
        bx3, bsz3, smx3, smsz3, bcorr3 = f.ratios(p3, t3, ne3, m_sum, 1, err_p3, err_ne3, err_t3, alpha)
        bx4, bsz4, smx4, smsz4, bcorr4 = f.ratios(p4, t4, ne4, m_sum, 1, err_p4, err_ne4, err_t4, alpha)
        bx5, bsz5, smx5, smsz5, bcorr5 = f.ratios(p5, t5, ne5, m_sum, 1, err_p5, err_ne5, err_t5, alpha)
        bx6, bsz6, smx6, smsz6, bcorr6 = f.ratios(p6, t6, ne6, m_sum, 1, err_p6, err_ne6, err_t6, alpha)
        bx7, bsz7, smx7, smsz7, bcorr7 = f.ratios(p7, t7, ne7, m_sum, 1, err_p7, err_ne7, err_t7, alpha)
        bx8, bsz8, smx8, smsz8, bcorr8 = f.ratios(p8, t8, ne8, m_sum, 1, err_p8, err_ne8, err_t8, alpha)

        bx_all, bsz_all, smx_all, smsz_all, bcorr_all = f.ratios(p_all, t_all, ne_all, m_sum, 1, err_p_all, err_ne_all, err_t_all, alpha)

        plt.plot(rlogbin_cen, bsz_all, label="All sectors", ls='solid', color='black',alpha=1)
        plt.scatter(rlogbin_cen, bsz_all, s=6, c='black', alpha=1)

        plt.plot(rlogbin_cen, bsz1, label="sector 1",ls='dashed', color='green',alpha=0.7)  # , color='orange', marker='.')
        #plt.plot(rlogbin_cen, bsz2, label="sector 2",ls='dashed', alpha=0.7)  # , color='orange', marker='.')
        #plt.plot(rlogbin_cen, bsz3, label="sector 3",ls='dashed', alpha=0.7)  # , color='orange', marker='.')
        #plt.plot(rlogbin_cen, bsz4, label="sector 4",ls='dashed',color='red', alpha=0.7)  # , color='orange', marker='.')
        plt.plot(rlogbin_cen, bsz5, label="sector 5",ls='dashed',color='purple', alpha=0.7)  # , color='orange', marker='.')
        #plt.plot(rlogbin_cen, bsz6, label="sector 6",ls='dashed', alpha=0.7)  # , color='orange', marker='.')
        #plt.plot(rlogbin_cen, bsz7, label="sector 7",ls='dashed', alpha=0.7)  # , color='orange', marker='.')
        #plt.plot(rlogbin_cen, bsz8, label="sector 8",ls='dashed', alpha=0.7)  # , color='orange', marker='.')

        plt.scatter(rlogbin_cen, bsz1, s=6, alpha=0.7, c='green')  # , color='orange', marker='.')
        #plt.scatter(rlogbin_cen, bsz2, s=6, alpha=0.7)  # , color='orange', marker='.')
        #plt.scatter(rlogbin_cen, bsz3, s=6, alpha=0.7)  # , color='orange', marker='.')
        #plt.scatter(rlogbin_cen, bsz4, s=6, c='red', alpha=0.7)  # , color='orange', marker='.')
        plt.scatter(rlogbin_cen, bsz5, s=6, c='purple', alpha=0.7)  # , color='orange', marker='.')
        #plt.scatter(rlogbin_cen, bsz6, s=6, alpha=0.7)  # , color='orange', marker='.')
        #plt.scatter(rlogbin_cen, bsz7, s=6, alpha=0.7)  # , color='orange', marker='.')
        #plt.scatter(rlogbin_cen, bsz8, s=6, alpha=0.7)  # , color='orange', marker='.')



        plt.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        plt.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,labelleft=True)
        # a.plot(rlogbin, zerohuit, color='grey')
        # a.plot(rlogbin, un, color='black')
        plt.xscale('log')
        plt.legend(prop={'size': 16})
        plt.ylabel("$(1-\mathrm{b})=rac{M_{\mathrm{HE}}}{M_{\mathrm{tot}}}$")

        un = [1 for i in range(40)]
        zerohuit = [0.8 for i in range(40)]
        un = np.array(un)
        zerohuit = np.array(zerohuit)

        plt.axvline(x=1087, color='grey')  # label='$R_{500}$')
        plt.axvline(x=2147, color='grey', ls='dashed')  # label='$R_{Vir}$')
        plt.text(1087, 0.83, "$R_{500}$", rotation=90, size=16)
        plt.text(2147, 0.83, "$R_{vir}$", rotation=90, size=16)
        plt.axvline(x=350, color='grey', ls='dotted')
        plt.text(350, 0.55, "AGN feedback", rotation=90, size=16)
        plt.axvline(x=850, color='grey', ls='dotted')
        plt.text(850, 0.55, "Matter infall", rotation=90, size=16)
        plt.plot(rlogbin, zerohuit, color='grey')
        plt.plot(rlogbin, un, color='black')
        plt.xlim(200, 4500)
        #plt.ylim(0.75, 2.45)
        plt.ylim(0.5,3.5)
        plt.xlabel("$R [\mathrm{kpc}]$")

        plt.show()

    def double_plot():

        r200m = 2895
        rcar = r200m

        fig, ((ax1), (ax2)) = plt.subplots(ncols=1, nrows=2)

        ax1.scatter(rlogbin_cen/rcar, 10 ** t_all, s=6, c='black')
        ax1.errorbar(rlogbin_cen/rcar, 10 ** t_all, yerr=10 ** t_all * err_t_all, ls='dashed', label='Full box', alpha=0.7,
                     c='black')


        # plt.scatter(rlogbin_cen, 10**t2, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, 10**t2, yerr=10**t2*err_t2, ls='dotted', label='sector 2', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, 10**t3, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, 10**t3, yerr=10**t3*err_t3, ls='dotted', label='sector 3', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, 10**t4, s=6, c='red')
        # plt.errorbar(rlogbin_cen, 10**t4, yerr=10**t4*err_t4, ls='dotted', label='sector 4', alpha=0.7, c='red')

        ax1.scatter(rlogbin_cen/rcar, 10 ** t5, s=6, c='red')
        ax1.errorbar(rlogbin_cen/rcar, 10 ** t5, yerr=10 ** t5 * err_t5, ls='dotted', label='Filament',
                     alpha=0.7, c='red')

        # plt.scatter(rlogbin_cen, 10**t6, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, 10**t6, yerr=10**t6*err_t6, ls='dotted', label='sector 6', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, 10**t7, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, 10**t7, yerr=10**t7*err_t7, ls='dotted', label='sector 7', alpha=0.7)  # , c='blue')

        ax1.scatter(rlogbin_cen/rcar, 10 ** t8, s=6, c='blue')
        ax1.errorbar(rlogbin_cen/rcar, 10 ** t8, yerr=10 ** t8 * err_t8, ls='dotted', label='Spherical Collapse', alpha=0.7,
                     c='blue')

        ax1.scatter(rlogbin_cen / rcar, 10 ** t1, s=6, c='green')  # , c='blue')
        ax1.errorbar(rlogbin_cen / rcar, 10 ** t1, yerr=10 ** t1 * err_t1, ls='dotted', label='Relaxed', alpha=0.7,
                     c='green')  # , c='blue')

        ax1.set_ylabel("$T$[keV]", size=16)
        ax1.set_xscale('log')
        ax1.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        ax1.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,
                        labelleft=True)
        #ax1.axvline(x=1087, color='grey')
        #ax1.text(1087, 0, "$\mathrm{R_{500}}$", rotation=90, size=16)
        #ax1.axvline(x=2147, color='grey', ls='dashed')
        #ax1.text(2147, 0, "$\mathrm{R_{vir}}$", rotation=90, size=16)

        # plt.axvline(x=350, color='grey', ls='dotted')
        # plt.text(350, 0, "AGN feedback", rotation=90, size=16)

        #plt.axvline(x=850, color='grey', ls='dotted')
        #plt.text(850, 0, "Matter infall", rotation=90, size=16)

        ax1.axvline(x=1.689, color='green', ls='dashed')
        ax1.axvline(x=1.3456, color='blue', ls='dashed')
        #ax1.axvline(x=0.23, color='red', ls = 'dashed')

        ax1.legend(prop={'size': 16})
        #plt.xlabel("$R$[kpc]", size=16)

        ax1.set_xticks([])

        ax1.set_title('Temperature (top) and Electron density (bottom) radial profiles', size=16)

        ax2.scatter(rlogbin_cen/rcar, ne_all, s=6, c='black')
        ax2.errorbar(rlogbin_cen/rcar, ne_all, yerr=err_ne_all, ls='dashed', alpha=0.7, c='black')

        # plt.scatter(rlogbin_cen, p2, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p2, yerr=err_p2, ls='dotted', label='sector 2', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, p3, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p3, yerr=err_p3, ls='dotted', label='sector 3', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, p4, s=6, c='red')
        # plt.errorbar(rlogbin_cen, p4, yerr=err_p4, ls='dotted', label='sector 4', alpha=0.7, c='red')

        ax2.scatter(rlogbin_cen/rcar, ne5, s=6, c='red')
        ax2.errorbar(rlogbin_cen/rcar, ne5, yerr=err_ne5, ls='dotted', alpha=0.7, c='red')

        # plt.scatter(rlogbin_cen, p6, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p6, yerr=err_p6, ls='dotted', label='sector 6', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, p7, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p7, yerr=err_p7, ls='dotted', label='sector 7', alpha=0.7)  # , c='blue')

        ax2.scatter(rlogbin_cen/rcar, ne8, s=6, c='blue')
        ax2.errorbar(rlogbin_cen/rcar, ne8, yerr=err_ne8, ls='dotted', alpha=0.7, c='blue')

        ax2.scatter(rlogbin_cen / rcar, ne1, s=6, c='green')  # , c='blue')
        ax2.errorbar(rlogbin_cen / rcar, ne1, yerr=err_ne1, ls='dotted', alpha=0.7,
                     c='green')  # , c='blue')

        ax2.set_ylabel("lo$\mathrm{g_{10}}$($n_e$[$\mathrm{cm^{-3}}$])", size=16)
        ax2.set_xscale('log')
        ax2.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        ax2.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,
                        labelleft=True)
        ax2.axvline(x=1.3456, color='blue', label="$R_{sp,spherical\,collapse}$ from pressure",ls='dashed')
        ax2.axvline(x=1.689, color='green', label="$R_{sp,relaxed}$ from pressure", ls='dashed')
        #ax2.axvline(x=0.23, color='red', label="Accretion shock in filament sector",ls='dashed')
        #plt.text(1087, -7, "$\mathrm{R_{500}}$", rotation=90, size=16)
        #plt.axvline(x=2147, color='grey', ls='dashed')
        #plt.text(2147, -7, "$\mathrm{R_{vir}}$", rotation=90, size=16)

        # plt.axvline(x=350, color='grey', ls='dotted')
        # plt.text(350, -7, "Shock front (AGN)", rotation=90, size=16)

        # plt.axvline(x=850, color='grey', ls='dotted')
        # plt.text(850, -7, "Shock front (Matter infall)", rotation=90, size=16)

        ax2.legend(prop={'size': 16})
        ax2.set_xlabel("$R/R_{200m}$", size=16)

        plt.subplots_adjust(wspace=0, hspace=0)



        plt.show()
        sys.exit()

    def triple_plot():

        r200m = 2895
        rcar = r200m

        #print("err_t5",err_t5)
        #sys.exit()

        fig, ((ax1), (ax3), (ax2)) = plt.subplots(ncols=1, nrows=3)

        ax1.scatter(rlogbin_cen / rcar, 10 ** t8, s=6, c='dodgerblue')
        ax1.errorbar(rlogbin_cen / rcar, 10 ** t8, yerr=10 ** t8 * err_t8, ls='solid', label='Collapsing Material', alpha=0.7, c='dodgerblue')

        ax1.scatter(rlogbin_cen / rcar, 10 ** t1, s=6, c='orange')  # , c='blue')
        ax1.errorbar(rlogbin_cen / rcar, 10 ** t1, yerr=10 ** t1 * err_t1, ls='solid', label='Outflowing Material',alpha=0.7,c='orange')  # , c='blue')

        ax1.scatter(rlogbin_cen / rcar, 10 ** t5, s=6, c='mediumvioletred')
        ax1.errorbar(rlogbin_cen / rcar, 10 ** t5, yerr=10 ** t5 * err_t5, ls='solid', label='Filament Material', alpha=0.7, c='mediumvioletred')

        ax1.scatter(rlogbin_cen / rcar, 10 ** t_all, s=6, c='black')
        ax1.errorbar(rlogbin_cen / rcar, 10 ** t_all, yerr=10 ** t_all * err_t_all, ls='solid', label='Full Material', alpha=0.7, c='black')

        # plt.scatter(rlogbin_cen, 10**t2, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, 10**t2, yerr=10**t2*err_t2, ls='dotted', label='sector 2', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, 10**t3, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, 10**t3, yerr=10**t3*err_t3, ls='dotted', label='sector 3', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, 10**t4, s=6, c='red')
        # plt.errorbar(rlogbin_cen, 10**t4, yerr=10**t4*err_t4, ls='dotted', label='sector 4', alpha=0.7, c='red')


        # plt.scatter(rlogbin_cen, 10**t6, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, 10**t6, yerr=10**t6*err_t6, ls='dotted', label='sector 6', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, 10**t7, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, 10**t7, yerr=10**t7*err_t7, ls='dotted', label='sector 7', alpha=0.7)  # , c='blue')


        ax1.set_ylabel("$T~[keV]$", size=22)
        ax1.set_xscale('log')
        ax1.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        ax1.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=18, labelbottom=True,
                        labelleft=True)
        # ax1.axvline(x=1087, color='grey')
        # ax1.text(1087, 0, "$\mathrm{R_{500}}$", rotation=90, size=16)
        # ax1.axvline(x=2147, color='grey', ls='dashed')
        # ax1.text(2147, 0, "$\mathrm{R_{vir}}$", rotation=90, size=16)

        # plt.axvline(x=350, color='grey', ls='dotted')
        # plt.text(350, 0, "AGN feedback", rotation=90, size=16)

        # plt.axvline(x=850, color='grey', ls='dotted')
        # plt.text(850, 0, "Matter infall", rotation=90, size=16)

        ax1.axvline(x=1, color='grey',ls='dashed')
        ax1.axvline(x=1.894, color='orange', ls='dashed')
        ax1.axvline(x=1.502, color='dodgerblue', ls='dashed')
        #ax1.axvline(x=0.23, color='red', ls='dashed')

        ax1.set_xlim(0.3, 4.9)

        ax1.legend(prop={'size': 20},loc="upper right")
        # plt.xlabel("$R$[kpc]", size=16)

        ax1.set_xticks([])

        ##adding baryons den rsp for comparison

        #ax1.axvline(x=1.689, color='green', ls='dotted')
        #ax1.axvline(x=1.3456, color='blue', ls='dotted')

        ax2.scatter(rlogbin_cen / rcar, k8, s=6, c='dodgerblue')
        ax2.errorbar(rlogbin_cen / rcar, k8, yerr=err_k8, ls='solid', alpha=0.7, c='dodgerblue')

        ax2.scatter(rlogbin_cen / rcar, k1, s=6, c='orange')  # , c='blue')
        ax2.errorbar(rlogbin_cen / rcar, k1, yerr=err_k1, ls='solid', alpha=0.7,c='orange')  # , c='blue')

        # plt.scatter(rlogbin_cen, p2, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p2, yerr=err_p2, ls='dotted', label='sector 2', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, p3, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p3, yerr=err_p3, ls='dotted', label='sector 3', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, p4, s=6, c='red')
        # plt.errorbar(rlogbin_cen, p4, yerr=err_p4, ls='dotted', label='sector 4', alpha=0.7, c='red')

        ax2.scatter(rlogbin_cen / rcar, k5, s=6, c='mediumvioletred')
        ax2.errorbar(rlogbin_cen / rcar, k5, yerr=err_k5, ls='solid', alpha=0.7, c='mediumvioletred')

        ax2.scatter(rlogbin_cen / rcar, k_all, s=6, c='black')
        ax2.errorbar(rlogbin_cen / rcar, k_all, yerr=err_k_all, ls='solid', alpha=0.7, c='black')

        # plt.scatter(rlogbin_cen, p6, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p6, yerr=err_p6, ls='dotted', label='sector 6', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, p7, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p7, yerr=err_p7, ls='dotted', label='sector 7', alpha=0.7)  # , c='blue')

        ax2.set_ylabel("lo$\mathrm{g_{10}}$($K~[\mathrm{keV~cm^2}$])", size=22)
        ax2.set_xscale('log')
        ax2.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        ax2.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=18, labelbottom=True,
                        labelleft=True)

        ax2.axvline(x=1, color='grey', ls='dashed') #label='$R_{200m}$'
        ax2.axvline(x=1.894, color='orange', ls='dashed') #label="$R_{sb}$ in relaxed sector",
        ax2.axvline(x=1.502, color='dodgerblue', ls='dashed') #label="$R_{sb}$ in spherical collapse sector",
        #ax2.axvline(x=0.23, color='red', label="Accretion shock in filament sector", ls='dashed')
        # plt.text(1087, -7, "$\mathrm{R_{500}}$", rotation=90, size=16)
        # plt.axvline(x=2147, color='grey', ls='dashed')
        # plt.text(2147, -7, "$\mathrm{R_{vir}}$", rotation=90, size=16)

        # plt.axvline(x=350, color='grey', ls='dotted')
        # plt.text(350, -7, "Shock front (AGN)", rotation=90, size=16)

        # plt.axvline(x=850, color='grey', ls='dotted')
        # plt.text(850, -7, "Shock front (Matter infall)", rotation=90, size=16)

        #ax2.legend(prop={'size': 16})

        ##adding baryons den rsp for comparison

        #ax2.axvline(x=1.689, color='green', ls='dotted', label="$R_{sp,relaxed}$ from bar density")
        #ax2.axvline(x=1.3456, color='blue', ls='dotted', label="$R_{sp,spherical\,collapse}$ from bar density")

        ax2.set_xlim(0.3, 4.9)


        ax3.scatter(rlogbin_cen / rcar, ne8, s=6, c='dodgerblue')
        ax3.errorbar(rlogbin_cen / rcar, ne8, yerr=err_ne8, ls='solid', alpha=0.7, c='dodgerblue')

        ax3.scatter(rlogbin_cen/ rcar, ne1, s=6, c='orange')  # , c='blue')
        ax3.errorbar(rlogbin_cen/ rcar, ne1, yerr=err_ne1, ls='solid', alpha=0.7,c='orange')  # , c='blue')

        # plt.scatter(rlogbin_cen, ne2, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, ne2, yerr=err_ne2, ls='dotted', label='sector 2', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, ne3, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, ne3, yerr=err_ne3, ls='dotted', label='sector 3', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, ne4, s=6 , c='red')
        # plt.errorbar(rlogbin_cen, ne4, yerr=err_ne4, ls='dotted', label='sector 4', alpha=0.7 , c='red')

        ax3.scatter(rlogbin_cen/ rcar, ne5, s=6, c='mediumvioletred')
        ax3.errorbar(rlogbin_cen/ rcar, ne5, yerr=err_ne5, ls='solid', alpha=0.7, c='mediumvioletred')

        ax3.scatter(rlogbin_cen / rcar, ne_all, s=6, c='black')
        ax3.errorbar(rlogbin_cen / rcar, ne_all, yerr=err_ne_all, ls='solid', alpha=0.7, c='black')

        # plt.scatter(rlogbin_cen, ne6, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, ne6, yerr=err_ne6, ls='dotted', label='sector 6', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, ne7, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, ne7, yerr=err_ne7, ls='dotted', label='sector 7', alpha=0.7)  # , c='blue')



        #ax3.axvline(x=1.689, color='green', label="$R_{sb}$ in relaxed sector", ls='dashed')
        #ax3.axvline(x=1.3456, color='blue', label="$R_{sb}$ in spherical collapse sector", ls='dashed')
        #ax3.axvline(x=0.23, color='red', label="Accretion shock in filament sector", ls='dashed')

        ax3.axvline(x=1, color='grey', label='$R_{200m}$', ls='dashed')
        ax3.axvline(x=1.502, color='dodgerblue', label="$R_{sp,pres,Collapsing}$", ls='dashed')
        ax3.axvline(x=1.894, color='orange', label="$R_{sp,pres,Outflowing}$", ls='dashed')

        ax3.set_ylabel("lo$\mathrm{g_{10}}$($n_{\mathrm{e}}~$[c$\mathrm{m^{-3}}$])", size=22)
        ax3.set_xscale('log')
        ax3.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        ax3.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=18, labelbottom=True,
                        labelleft=True)
        ax3.set_ylim(-7.5, -2.5)
        ax3.set_xlim(0.3, 4.9)

        ##adding baryons den rsp for comparison

        #ax3.axvline(x=1.689, color='green', ls='dotted')
        #ax3.axvline(x=1.3456, color='blue', ls='dotted')

        #ax3.axvline(x=1087, color='grey')
        #ax3.text(1087, -6, "$\mathrm{R_{500}}$", rotation=90, size=16)
        #ax3.axvline(x=2147, color='grey', ls='dashed')
        #ax3.text(2147, -6, "$\mathrm{R_{vir}}$", rotation=90, size=16)

        # plt.axvline(x=350, color='grey', ls='dotted')
        # plt.text(350, -6, "Shock front (AGN)", rotation=90, size=16)

        # plt.axvline(x=850, color='grey', ls='dotted')
        # plt.text(850, -6, "Shock front (Matter infall)", rotation=90, size=16)

        #plt.legend(prop={'size': 16})
        #plt.xlabel("$R$[kpc]", size=16)
        # plt.set_xticks([])

        ax2.set_xlabel("$R/R_{200m}$", size=22)

        ax3.legend(prop={'size': 20},loc="lower left")
        #ax2.legend(prop={'size': 16})


        plt.subplots_adjust(wspace=0, hspace=0)

        #fig.suptitle('Radial profiles, temperature (top), pressure (middle) and electron density (bottom)', size=16)

        plt.show()
        sys.exit()

    def plot_pressure_grad():

        r200m = 2895

        r_car = r200m

        r_car = 1

        b = 0.025

        grad_p_all = np.gradient(p_all, b)
        grad_p1 = np.gradient(p1, b)
        grad_p5 = np.gradient(p5, b)
        grad_p8 = np.gradient(p8, b)

        grad_p_all_med = np.gradient(p_all_med, b)
        grad_p1_med = np.gradient(p1_med, b)
        grad_p5_med = np.gradient(p5_med, b)
        grad_p8_med = np.gradient(p8_med, b)

        f, (c, d) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [2, 2]},sharex=True)

        c.scatter(rlogbin_cen / r_car, p8, s=6, c='dodgerblue')
        c.errorbar(rlogbin_cen / r_car, p8, yerr=err_p8, ls='solid', label='Collapsing Material', alpha=0.7, c='dodgerblue')

        c.scatter(rlogbin_cen / r_car, p1, s=6, c='orange')  # , c='blue')
        c.errorbar(rlogbin_cen / r_car, p1, yerr=err_p1, ls='solid', label='Outflowing Material', alpha=0.7, c='orange')  # , c='blue')

        c.scatter(rlogbin_cen / r_car, p5, s=6, c='mediumvioletred')
        c.errorbar(rlogbin_cen / r_car, p5, yerr=err_p5, ls='solid', label='Filament Material', alpha=0.7, c='mediumvioletred')

        c.scatter(rlogbin_cen / r_car, p_all, s=6, c='black')
        c.errorbar(rlogbin_cen / r_car, p_all, yerr=err_p_all, ls='solid', label='Full Material', alpha=0.7, c='black')

        #c.scatter(rlogbin_cen / r_car, p8_med, s=6, c='darkblue')
        #c.errorbar(rlogbin_cen / r_car, p8_med, yerr=err_p8_med, ls='dashed', label='Collapsing Material (median)', alpha=0.7, c='darkblue')

        #c.scatter(rlogbin_cen / r_car, p1_med, s=6, c='orangered')  # , c='blue')
        #c.errorbar(rlogbin_cen / r_car, p1_med, yerr=err_p1_med, ls='dashed', label='Outflowing Material (median)', alpha=0.7, c='orangered')  # , c='blue')

        #c.scatter(rlogbin_cen / r_car, p5_med, s=6, c='red')
        #c.errorbar(rlogbin_cen / r_car, p5_med, yerr=err_p5_med, ls='dashed', label='Filament Material (median)', alpha=0.7, c='red')

        #c.scatter(rlogbin_cen / r_car, p_all_med, s=6, c='black')
        #c.errorbar(rlogbin_cen / r_car, p_all_med, yerr=err_p_all_med, ls='dashed', label='Full Material (median)', alpha=0.7, c='black')

        # plt.scatter(rlogbin_cen, p2, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p2, yerr=err_p2, ls='dotted', label='sector 2', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, p3, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p3, yerr=err_p3, ls='dotted', label='sector 3', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, p4, s=6, c='red')
        # plt.errorbar(rlogbin_cen, p4, yerr=err_p4, ls='dotted', label='sector 4', alpha=0.7, c='red')



        # plt.scatter(rlogbin_cen, p6, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p6, yerr=err_p6, ls='dotted', label='sector 6', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, p7, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p7, yerr=err_p7, ls='dotted', label='sector 7', alpha=0.7)  # , c='blue')

        c.set_ylabel("lo$\mathrm{g_{10}}$($P~[\mathrm{keV~cm^{-3}}$])", size=22)
        c.set_xscale('log')
        c.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        c.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=18, labelbottom=True,
                        labelleft=True)

        #c.axvline(x=1, color='grey', ls='dashed')
        # d.axvline(x=1.605, color='black', label="All: $R_{sb}\sim$3450kpc", ls='dashed')
        #c.axvline(x=1.502, color='dodgerblue', ls='dashed')
        #c.axvline(x=1.894, color='orange', ls='dashed')

        c.axvline(x=1087, color='grey', ls='dashed')
        c.axvline(x=2147, color='grey', ls='dotted')
        c.set_ylim(-7,-1.5)

        #c.axvline(x=1087, color='grey')
        #c.text(1087, -7, "$\mathrm{R_{500}}$", rotation=90, size=16)
        #plt.axvline(x=2147, color='grey', ls='dashed')
        #plt.text(2147, -7, "$\mathrm{R_{vir}}$", rotation=90, size=16)

        # plt.axvline(x=350, color='grey', ls='dotted')
        # plt.text(350, -7, "Shock front (AGN)", rotation=90, size=16)

        # plt.axvline(x=850, color='grey', ls='dotted')
        # plt.text(850, -7, "Shock front (Matter infall)", rotation=90, size=16)

        c.legend(prop={'size': 12})
        #plt.xlabel("$R$[kpc]", size=16)
        # plt.set_xticks([])

        c.set_xscale('log')
        # c.set_title("Baryons density profiles (including galaxies)")
        c.xaxis.set_ticklabels([])
        #c.set_xlim(0.3, 4.9)
        # plt.xticks(fontsize=14)
        # a.set_yticks()

        # d.plot(rlogbin_cen/rvir, grad_log_den_bar, ls="dotted", marker='.',c="black")
        d.plot(rlogbin_cen / r_car, grad_p8, ls="solid", marker='.', c="dodgerblue")
        d.plot(rlogbin_cen / r_car, grad_p1, ls="solid", marker='.', c="orange")
        d.plot(rlogbin_cen / r_car, grad_p5, ls="solid", marker='.', c="mediumvioletred")
        d.plot(rlogbin_cen / r_car, grad_p_all, ls="solid", marker='.', c="black")

        #d.plot(rlogbin_cen / r_car, grad_p8_med, ls="dashed", marker='.', c="darkblue")
        #d.plot(rlogbin_cen / r_car, grad_p1_med, ls="dashed", marker='.', c="orangered")
        #d.plot(rlogbin_cen / r_car, grad_p5_med, ls="dashed", marker='.', c="red")
        #d.plot(rlogbin_cen / r_car, grad_p_all_med, ls="dashed", marker='.', c="black")

        d.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        d.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=18, labelbottom=True,
                      labelleft=True)

        #d.axvline(x=1, color='grey', label='$R_{200m}$', ls='dashed')
        # d.axvline(x=1.605, color='black', label="All: $R_{sb}\sim$3450kpc", ls='dashed')
        #d.axvline(x=1.502, color='dodgerblue', label="$R_{sp,Collapsing\,Material}$", ls='dashed') #$R_{sp,spherical\,collapse}\sim$4.1Mpc"
        #d.axvline(x=1.894, color='orange', label="$R_{sp,Outflowing\,Material}$", ls='dashed') #"$R_{sp,relaxed}\sim$4.9Mpc"
        d.axvline(x=1087, color='grey', label='$R_{500c}$', ls='dashed')
        d.axvline(x=2147, color='grey', label='$R_{vir}$', ls='dotted')

        #d.axvline(x=1087, color='grey')
        # d.text(1087, -7, "$R_{500}$", rotation=90, size=16)
        # d.axvline(x=2147, color='grey', ls='dashed')
        # d.text(2147, -7, "$R_{vir}$", rotation=90, size=16)
        # d.set_xlabel("$R$ [kpc]", size=16)
        #d.set_xlabel("$R/R_{200m}$", size=22)
        d.set_xlabel("R [kpc]", size=22)
        # d.yaxis.set_ticklabels([])
        # plt.ylabel("$log_{10}(P[keV/cm^3])$", size=16)
        # plt.ylabel("$log_{10}(n_e[cm^{-3}])$",size=14)
        # plt.ylabel("$log_{10}(T[keV])$", size=14)
        d.set_ylabel("$\\frac{dlog_{10}(P)}{dlog_{10}(r)}$",size=22)  # code à copier ds plot : $\frac{dlog_{10}(P)}{dlog_{10}(r)}$
        # plt.legend(prop={'size': 16})

        #d.axvline(x=1, color='grey', label='$R_{200m}$=2895kpc', ls='dashed')
        # d.axvline(x=1.605, color='black', label="All: $R_{sb}\sim$3450kpc", ls='dashed')
        #d.axvline(x=1.689, color='green', label="Relaxed: $R_{sb}\sim$4890kpc", ls='dashed')
        #d.axvline(x=1.3456, color='blue', label="Spherical collapse: $R_{sb}\sim$3890kpc", ls='dashed')

        d.set_xscale('log')
        #d.set_xlim(0.3, 4.9)
        d.set_xlim(200, 4000)
        #d.set_ylim(-10, 2.5)
        d.set_ylim(-11,1)
        zerox = np.linspace(50, 8000, 500)
        zeroy = np.zeros(500)
        d.legend(prop={'size': 12})
        # d.plot(zerox, zeroy, c='orange', ls='dotted')
        # b.set_xticks()
        # b.set_yticks()
        # plt.title("Projected temperature profiles, high res, MW los, 100 MC loops, \n comparison between 3 weighting methods")
        # plt.xlim(50,7000)
        plt.subplots_adjust(wspace=0, hspace=0)
        plt.show()
        print("end")
        sys.exit()

    def plot_density_grad():

        r200m = 2895

        r_car = r200m

        b = 0.025

        grad_ne_all = np.gradient(ne_all, b)
        grad_ne1 = np.gradient(ne1, b)
        grad_ne5 = np.gradient(ne5, b)
        grad_ne8 = np.gradient(ne8, b)

        f, (c, d) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [2, 2]})

        c.scatter(rlogbin_cen / r_car, ne_all, s=6, c='black')
        c.errorbar(rlogbin_cen / r_car, ne_all, yerr=err_ne_all, ls='dashed', label='Full Material', alpha=0.7, c='black')

        c.scatter(rlogbin_cen / r_car, ne5, s=6, c='red')
        c.errorbar(rlogbin_cen / r_car, ne5, yerr=err_ne5, ls='dotted', label='Filament Material', alpha=0.7,c='red')

        c.scatter(rlogbin_cen / r_car, ne8, s=6, c='blue')
        c.errorbar(rlogbin_cen / r_car, ne8, yerr=err_ne8, ls='dotted', label='Collapsing Material', alpha=0.7, c='blue')

        c.scatter(rlogbin_cen / r_car, ne1, s=6, c='green')  # , c='blue')
        c.errorbar(rlogbin_cen / r_car, ne1, yerr=err_ne1, ls='dotted', label='Outflowing Material', alpha=0.7, c='green')  # , c='blue')

        # plt.scatter(rlogbin_cen, p2, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p2, yerr=err_p2, ls='dotted', label='sector 2', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, p3, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p3, yerr=err_p3, ls='dotted', label='sector 3', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, p4, s=6, c='red')
        # plt.errorbar(rlogbin_cen, p4, yerr=err_p4, ls='dotted', label='sector 4', alpha=0.7, c='red')



        # plt.scatter(rlogbin_cen, p6, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p6, yerr=err_p6, ls='dotted', label='sector 6', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, p7, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p7, yerr=err_p7, ls='dotted', label='sector 7', alpha=0.7)  # , c='blue')

        c.set_ylabel("lo$\mathrm{g_{10}}$($n_e$[$\mathrm{cm^{-3}}$])", size=18)
        c.set_xscale('log')
        c.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        c.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=18, labelbottom=True,
                        labelleft=True)

        c.axvline(x=1, color='grey', ls='dashed')
        # d.axvline(x=1.605, color='black', label="All: $R_{sb}\sim$3450kpc", ls='dashed')
        #c.axvline(x=1.418, color='blue', ls='dashed')
        #c.axvline(x=1.812, color='green', ls='dashed')

        #c.axvline(x=1087, color='grey')
        #c.text(1087, -7, "$\mathrm{R_{500}}$", rotation=90, size=16)
        #plt.axvline(x=2147, color='grey', ls='dashed')
        #plt.text(2147, -7, "$\mathrm{R_{vir}}$", rotation=90, size=16)

        # plt.axvline(x=350, color='grey', ls='dotted')
        # plt.text(350, -7, "Shock front (AGN)", rotation=90, size=16)

        # plt.axvline(x=850, color='grey', ls='dotted')
        # plt.text(850, -7, "Shock front (Matter infall)", rotation=90, size=16)

        c.legend(prop={'size': 16})
        #plt.xlabel("$R$[kpc]", size=16)
        # plt.set_xticks([])

        c.set_xscale('log')
        # c.set_title("Baryons density profiles (including galaxies)")
        c.xaxis.set_ticklabels([])
        c.set_xlim(0.3, 10)
        c.set_ylim(-8.5,-2)
        # plt.xticks(fontsize=14)
        # a.set_yticks()

        # d.plot(rlogbin_cen/rvir, grad_log_den_bar, ls="dotted", marker='.',c="black")
        d.plot(rlogbin_cen / r_car, grad_ne_all, ls="dotted", marker='.', c="black")
        d.plot(rlogbin_cen / r_car, grad_ne5, ls="dotted", marker='.', c="red")
        d.plot(rlogbin_cen / r_car, grad_ne8, ls="dotted", marker='.', c="blue")
        d.plot(rlogbin_cen / r_car, grad_ne1, ls="dotted", marker='.', c="green")
        d.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        d.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=18, labelbottom=True,
                      labelleft=True)

        d.axvline(x=1, color='grey', label='$R_{200m}$', ls='dashed')
        # d.axvline(x=1.605, color='black', label="All: $R_{sb}\sim$3450kpc", ls='dashed')
        #d.axvline(x=1.418, color='blue', label="$R_{sp,spherical\,collapse}\sim$4.1Mpc", ls='dashed')
        #d.axvline(x=1.812, color='green', label="$R_{sp,relaxed}\sim$5.2Mpc", ls='dashed')

        ##adding baryons den rsp for comparison

        #d.axvline(x=1, color='grey', ls='dashed')
        d.axvline(x=1.894, color='green', ls='dashed', label="from pressure in relaxed")
        d.axvline(x=1.502, color='blue', ls='dashed', label="from pressure in spher coll")
        d.axvline(x=1.689, color='green', ls='solid', label="from bar density in relaxed")
        d.axvline(x=1.3456, color='blue', ls='solid', label="from bar density in spher coll")

        #d.axvline(x=1087, color='grey')
        # d.text(1087, -7, "$R_{500}$", rotation=90, size=16)
        # d.axvline(x=2147, color='grey', ls='dashed')
        # d.text(2147, -7, "$R_{vir}$", rotation=90, size=16)
        # d.set_xlabel("$R$ [kpc]", size=16)
        d.set_xlabel("$R/R_{200m}$", size=18)
        # d.yaxis.set_ticklabels([])
        # plt.ylabel("$log_{10}(P[keV/cm^3])$", size=16)
        # plt.ylabel("$log_{10}(n_e[cm^{-3}])$",size=14)
        # plt.ylabel("$log_{10}(T[keV])$", size=14)
        d.set_ylabel("$dlog_{10}(n_e)/log_{10}(r)$",size=18)  # code à copier ds plot : $\frac{dlog_{10}(n_e)}{dlog_{10}(r)}$
        # plt.legend(prop={'size': 16})

        #d.axvline(x=1, color='grey', label='$R_{200m}$=2895kpc', ls='dashed')
        # d.axvline(x=1.605, color='black', label="All: $R_{sb}\sim$3450kpc", ls='dashed')
        #d.axvline(x=1.689, color='green', label="Relaxed: $R_{sb}\sim$4890kpc", ls='dashed')
        #d.axvline(x=1.3456, color='blue', label="Spherical collapse: $R_{sb}\sim$3890kpc", ls='dashed')

        d.set_xscale('log')
        d.set_xlim(0.3, 10)
        d.set_ylim(-8, 8)
        zerox = np.linspace(50, 8000, 500)
        zeroy = np.zeros(500)
        d.legend()
        # d.plot(zerox, zeroy, c='orange', ls='dotted')
        # b.set_xticks()
        # b.set_yticks()
        # plt.title("Projected temperature profiles, high res, MW los, 100 MC loops, \n comparison between 3 weighting methods")
        # plt.xlim(50,7000)
        plt.subplots_adjust(wspace=0, hspace=0)
        plt.show()
        print("end")
        sys.exit()

    def plot_entropy_grad():
        r200m = 2895

        r_car = r200m

        b = 0.025

        grad_k_all = np.gradient(k_all, b)
        grad_k1 = np.gradient(k1, b)
        grad_k5 = np.gradient(k5, b)
        grad_k8 = np.gradient(k8, b)

        f, (c, d) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [2, 2]})

        c.scatter(rlogbin_cen / r_car, k_all, s=6, c='black')
        c.errorbar(rlogbin_cen / r_car, k_all, yerr=err_k_all, ls='dashed', label='Full box', alpha=0.7, c='black')

        c.scatter(rlogbin_cen / r_car, k1, s=6, c='green')  # , c='blue')
        c.errorbar(rlogbin_cen / r_car, k1, yerr=err_k1, ls='dotted', label='Relaxed', alpha=0.7,
                   c='green')  # , c='blue')

        # plt.scatter(rlogbin_cen, p2, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p2, yerr=err_p2, ls='dotted', label='sector 2', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, p3, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p3, yerr=err_p3, ls='dotted', label='sector 3', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, p4, s=6, c='red')
        # plt.errorbar(rlogbin_cen, p4, yerr=err_p4, ls='dotted', label='sector 4', alpha=0.7, c='red')

        c.scatter(rlogbin_cen / r_car, k5, s=6, c='red')
        c.errorbar(rlogbin_cen / r_car, k5, yerr=err_k5, ls='dotted', label='Matter infall from filament', alpha=0.7,
                   c='red')

        # plt.scatter(rlogbin_cen, p6, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p6, yerr=err_p6, ls='dotted', label='sector 6', alpha=0.7)  # , c='blue')

        # plt.scatter(rlogbin_cen, p7, s=6)  # , c='blue')
        # plt.errorbar(rlogbin_cen, p7, yerr=err_p7, ls='dotted', label='sector 7', alpha=0.7)  # , c='blue')

        c.scatter(rlogbin_cen / r_car, k8, s=6, c='blue')
        c.errorbar(rlogbin_cen / r_car, k8, yerr=err_k8, ls='dotted', label='Spherical collapse', alpha=0.7, c='blue')

        c.set_ylabel("lo$\mathrm{g_{10}}$($K$[keV.c$\mathrm{m^{2}}$])", size=16)
        c.set_xscale('log')
        c.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        c.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,
                      labelleft=True)
        # c.axvline(x=1087, color='grey')
        # c.text(1087, -7, "$\mathrm{R_{500}}$", rotation=90, size=16)
        # plt.axvline(x=2147, color='grey', ls='dashed')
        # plt.text(2147, -7, "$\mathrm{R_{vir}}$", rotation=90, size=16)

        # plt.axvline(x=350, color='grey', ls='dotted')
        # plt.text(350, -7, "Shock front (AGN)", rotation=90, size=16)

        # plt.axvline(x=850, color='grey', ls='dotted')
        # plt.text(850, -7, "Shock front (Matter infall)", rotation=90, size=16)

        c.legend(prop={'size': 16})
        # plt.xlabel("$R$[kpc]", size=16)
        # plt.set_xticks([])

        c.set_xscale('log')
        # c.set_title("Baryons density profiles (including galaxies)")
        c.xaxis.set_ticklabels([])
        c.set_xlim(0.3, 10)
        # plt.xticks(fontsize=14)
        # a.set_yticks()

        # d.plot(rlogbin_cen/rvir, grad_log_den_bar, ls="dotted", marker='.',c="black")
        d.plot(rlogbin_cen / r_car, grad_k_all, ls="dotted", marker='.', c="black")
        d.plot(rlogbin_cen / r_car, grad_k1, ls="dotted", marker='.', c="green")
        d.plot(rlogbin_cen / r_car, grad_k5, ls="dotted", marker='.', c="red")
        d.plot(rlogbin_cen / r_car, grad_k8, ls="dotted", marker='.', c="blue")
        d.tick_params(which="minor", right=True, top=True, direction="in", length=3, labelbottom=False)
        d.tick_params(which="major", right=True, top=True, direction="in", length=5, labelsize=16, labelbottom=True,
                      labelleft=True)
        # d.axvline(x=1087, color='grey')
        # d.text(1087, -7, "$R_{500}$", rotation=90, size=16)
        # d.axvline(x=2147, color='grey', ls='dashed')
        # d.text(2147, -7, "$R_{vir}$", rotation=90, size=16)
        # d.set_xlabel("$R$ [kpc]", size=16)
        d.set_xlabel("$R/R_{200m}$", size=16)
        # d.yaxis.set_ticklabels([])
        # plt.ylabel("$log_{10}(P[keV/cm^3])$", size=16)
        # plt.ylabel("$log_{10}(n_e[cm^{-3}])$",size=14)
        # plt.ylabel("$log_{10}(T[keV])$", size=14)
        d.set_ylabel("$dlog_{10}(K)/log_{10}(r)$",
                     size=16)  # code à copier ds plot : $\frac{dlog_{10}(\rho_{bar})}{dlog_{10}(r)}$
        # plt.legend(prop={'size': 16})

        # d.axvline(x=1, color='grey', label='$R_{200m}$=2895kpc', ls='dashed')
        # d.axvline(x=1.605, color='black', label="All: $R_{sb}\sim$3450kpc", ls='dashed')
        # d.axvline(x=1.689, color='green', label="Relaxed: $R_{sb}\sim$4890kpc", ls='dashed')
        # d.axvline(x=1.3456, color='blue', label="Spherical collapse: $R_{sb}\sim$3890kpc", ls='dashed')

        #test
        d.axvline(x=1.894, color='green', ls='dashed', label="from pressure in relaxed")
        d.axvline(x=1.502, color='blue', ls='dashed', label="from pressure in spher coll")
        d.axvline(x=1.689, color='green', ls='solid', label="from bar density in relaxed")
        d.axvline(x=1.3456, color='blue', ls='solid', label="from bar density in spher coll")

        d.set_xscale('log')
        d.set_xlim(0.3, 10)
        # d.set_ylim(-10, 2.5)
        zerox = np.linspace(50, 8000, 500)
        zeroy = np.zeros(500)
        d.legend()
        # d.plot(zerox, zeroy, c='orange', ls='dotted')
        # b.set_xticks()
        # b.set_yticks()
        # plt.title("Projected temperature profiles, high res, MW los, 100 MC loops, \n comparison between 3 weighting methods")
        # plt.xlim(50,7000)
        plt.subplots_adjust(wspace=0, hspace=0)
        plt.show()
        print("end")
        sys.exit()

    def compare_gradient():
        #log_den_bar_all = np.load('/data/cluster/tlebeau/virgo/sectors_study/splashback/log_den_bar_all.npy')
        #grad_log_den_bar_all = np.load('/data/cluster/tlebeau/virgo/sectors_study/splashback/grad_log_den_bar_all.npy')

        #log_den_bar_sec5 = np.load('/data/cluster/tlebeau/virgo/sectors_study/splashback/log_den_bar_sec5.npy')
        #grad_log_den_bar_sec5 = np.load('/data/cluster/tlebeau/virgo/sectors_study/splashback/grad_log_den_bar_sec5.npy')

        #log_den_bar_sec1 = np.load('/data/cluster/tlebeau/virgo/sectors_study/splashback/log_den_bar_sec1.npy')
        grad_log_den_bar_sec1 = np.load('/data/cluster/tlebeau/virgo/sectors_study/splashback/grad_log_den_bar_sec1.npy')

        #log_den_bar_sec8 = np.load('/data/cluster/tlebeau/virgo/sectors_study/splashback/log_den_bar_sec8.npy')
        grad_log_den_bar_sec8 = np.load('/data/cluster/tlebeau/virgo/sectors_study/splashback/grad_log_den_bar_sec8.npy')

        grad_log_den_dm_sec1 = np.load('/data/cluster/tlebeau/virgo/sectors_study/splashback/grad_log_den_dm_sec1.npy')
        grad_log_den_dm_sec8 = np.load('/data/cluster/tlebeau/virgo/sectors_study/splashback/grad_log_den_dm_sec8.npy')

        r200m = 2895

        r_car = r200m

        b = 0.025

        #grad_p_all = np.gradient(p_all, b)
        grad_p1 = np.gradient(p1, b)
        #grad_p5 = np.gradient(p5, b)
        grad_p8 = np.gradient(p8, b)

        #grad_k_all = np.gradient(k_all, b)
        grad_k1 = np.gradient(k1, b)
        #grad_k5 = np.gradient(k5, b)
        grad_k8 = np.gradient(k8, b)

        grad_t1 = np.gradient(t1, b)
        grad_t8 = np.gradient(t8, b)

        plt.plot(rlogbin_cen / r_car, grad_k1, ls="solid", marker='.', c="cyan", label="$K$")
        plt.plot(rlogbin_cen_bar / r_car, grad_log_den_dm_sec1, ls="solid", marker='.', c="black",label=r"$\rho_{DM}$")
        plt.plot(rlogbin_cen_bar / r_car, grad_log_den_bar_sec1, ls="solid", marker='.', c="red", label=r"$\rho_{bar}$")
        plt.plot(rlogbin_cen / r_car, grad_p1, ls="solid", marker='.', c="orange", label="$P$")
        plt.plot(rlogbin_cen / r_car, grad_t1, ls="solid", marker='.', c="pink", label="$T$")

        #plt.plot(rlogbin_cen / r_car, grad_k8, ls="solid", marker='.', c="cyan", label="$K$")
        #plt.plot(rlogbin_cen_bar / r_car, grad_log_den_dm_sec8, ls="solid", marker='.', c="black", label=r"$\rho_{DM}$")
        #plt.plot(rlogbin_cen_bar / r_car, grad_log_den_bar_sec8, ls="solid", marker='.', c="red",label=r"$\rho_{bar}$")
        #plt.plot(rlogbin_cen / r_car, grad_p8, ls="solid", marker='.', c="orange", label="$P$")
        #plt.plot(rlogbin_cen / r_car, grad_t8, ls="solid", marker='.', c="pink", label="$T$")

        plt.xlim(0.3, 4.9)
        plt.ylim(-30,10)

        #for sector 1

        plt.axvline(x=1.505, color='black', label="$R_{sp,DM}$", ls='dashed')
        plt.axvline(x=1.683, color='red', ls='dashed', label="$R_{sp,bar}$")
        plt.axvline(x=1.886, color='orange', ls='dashed', label="$R_{sp,P}$")
        plt.axvline(x=1, color='grey', label='$R_{200m}$', ls='dashed')
        #plt.axvline(x=1.990, color='cyan', ls='dashed', label="$R_{sp,K}$")
        #plt.axvline(x=1.886, color='pink', ls='dashed', label="$R_{sp,T}$")

        #for sector 8
        #plt.axvline(x=1.19, color='black', label="$R_{sp,DM}$", ls='dashed')
        #plt.axvline(x=1.335, color='red', ls='dashed', label="$R_{sp,bar}$")
        #plt.axvline(x=1.503, color='orange', ls='dashed', label="$R_{sp,P}$")
        #plt.axvline(x=1, color='grey', label='$R_{200m}$', ls='dashed')
        #plt.axvline(x=1.503, color='cyan', ls='dashed', label="$R_{sp,K}$")

        plt.xlabel("$R/R_{200m}$", size=16)
        plt.xscale('log')
        plt.ylabel(r"$\frac{dlog(x)}{dlog(r)}$", size=16)
        plt.title("Outflowing material")
        #plt.title("Collapsing material")
        plt.legend(ncol=2)
        plt.show()

    def bias_w_prop_from_sectors_77():

        def mhex(i, tlog, ne_log, unit):
            r_bin = i
            b = 0.025
            dt = np.gradient(tlog, b)
            dne = np.gradient(ne_log, b)
            dt = dt[r_bin]
            dne = dne[r_bin]
            if unit == 0:
                # m = (((-10 ** ((i + 35.5) * 0.05) * kb * 10 ** tlog[r_bin]) / (g * mu * mp)) * (dne + dt)) * (
                #        (1E3 * pc) / m_sun)

                m = (((-10 ** ((18 + i) * 0.1) * kb * 10 ** tlog[r_bin]) / (g * mu * mp)) * (dne + dt)) * (
                        (1E3 * pc) / m_sun)

                sm = 0
                dne = 0
                dt = 0

            if unit == 1:
                m = (((-10 ** ((i + 35) * 0.05) * 10 ** tlog[r_bin]) / (g * mu * mp)) * (dne + dt)) * (
                        (1E3 * pc * kev) / m_sun)

                invt = np.array([1 / (10 ** tlog[i]) for i in range(len(tlog))])
                # print("invt", invt)
                dinvt = np.gradient(invt, b)
                dinvt = dinvt[r_bin]

                invne = np.array([1 / (10 ** ne_log[i]) for i in range(len(ne_log))])
                # print("invp", invne)
                dinvne = np.gradient(invne, b)
                dinvne = dinvne[r_bin]

                dmdne = -((10 ** ((i + 35) * 0.05) * dinvne * 10 ** tlog[i]) / (g * mu * mp)) * (
                        (1E3 * pc * kev) / m_sun)
                dmdt = m / (10 ** tlog[r_bin]) - ((10 ** ((i + 35) * 0.05) * dinvt * 10 ** tlog[i]) / (g * mu * mp)) * (
                        (1E3 * pc * kev) / m_sun)
                sm = np.sqrt(
                    dmdne ** 2 * (10 ** sne[i] * 10 ** ne_log[i]) ** 2 + dmdt ** 2 * (10 ** st[i] * 10 ** tlog[i]) ** 2)
                # print("sm", sm)
                # print("m", m)

            if unit == 2:
                # m = (((-10 ** ((i + 35.5) * 0.05) * kb * 10 ** tlog[r_bin]) / (g * mu * mp)) * (dne + dt)) * (
                #        (1E3 * pc) / m_sun)

                m = (((-10 ** ((92.5 + i) * 0.025) * kb * 10 ** tlog[r_bin]) / (g * mu * mp)) * (dne + dt)) * (
                        (1E3 * pc) / m_sun)

                sm = 0
                dne = 0
                dt = 0

            return m, sm, dne, dt

        def mhesz(i, plog, ne_log, unit):
            b = 0.025
            r_bin = i
            dp = np.gradient(plog, b)
            dp = dp[r_bin]
            if unit == 0:
                # m = ((-10 ** ((i + 35.5) * 0.05) * 10 ** plog[r_bin] * dp) / (g * mu * mp * 10 ** ne_log[r_bin])) * (
                #        (1E3 * pc * kev) / m_sun)

                m = ((-10 ** ((18 + i) * 0.1) * 10 ** plog[r_bin] * dp) / (g * mu * mp * 10 ** ne_log[r_bin])) * (
                        (1E3 * pc * kev) / m_sun)

                sm = 0
                dp = 0

            if unit == 1:
                m = ((-10 ** ((i + 35) * 0.05) * 10 ** plog[r_bin] * dp) / (g * mu * mp * 10 ** ne_log[r_bin])) * (
                        (1E3 * pc * kev) / m_sun)

                invp = np.array([1 / (10 ** plog[i]) for i in range(len(plog))])
                # print("invp",invp)
                dinvp = np.gradient(invp, b)
                dinvp = dinvp[r_bin]
                # print("dinvp",dinvp)
                dmdne = -m / (10 ** ne_log[r_bin])
                dmdp = m / (10 ** plog[r_bin]) - (
                        (10 ** ((i + 35) * 0.05) * dinvp * 10 ** plog[i]) / (g * mu * mp * 10 ** ne_log[r_bin])) * (
                               (1E3 * pc * kev) / m_sun)
                # dmdp=0
                sm = np.sqrt(
                    dmdne ** 2 * (10 ** sne[i] * 10 ** ne_log[i]) ** 2 + dmdp ** 2 * (10 ** sp[i] * 10 ** plog[i]) ** 2)
                # print("sne",sne[i])
                # print("sp",sp[i])
                # print("dmdne",dmdne)
                # print("dmdp",dmdp)
                # print("sm",sm)
                # print("m ",m)
                # sys.exit()

            if unit == 2:
                # m = ((-10 ** ((i + 35.5) * 0.05) * 10 ** plog[r_bin] * dp) / (g * mu * mp * 10 ** ne_log[r_bin])) * (
                #        (1E3 * pc * kev) / m_sun)

                m = ((-10 ** ((92.5 + i) * 0.025) * 10 ** plog[r_bin] * dp) / (g * mu * mp * 10 ** ne_log[r_bin])) * (
                        (1E3 * pc * kev) / m_sun)

                sm = 0
                dp = 0
            return m, sm, dp

        def compar_bias(i, plog, tlog, ne_log):
            mxlog, smx, dne, dt = mhex(i, tlog, ne_log, 2)
            # print("mhexlog", f"{mxlog:.3e}", "M_sun")
            mszlog, smsz, dp = mhesz(i, plog, ne_log, 2)

            # print("i",i,"msz", f"{mszlog:.3e}"," mcorr", f"{mcorrlog:.3e}", " mrand", f"{mrandlog:.3e}", " mrot", f"{mrotlog:.3e}", " mall", f"{malllog:.3e}"," mtot[i]",f"{m_tot[i]:.3e}", "M_sun")

            # sys.exit()

            # print("mrand", f"{mrandlog:.3e}", "M_sun")
            # return mxlog / m_tot[i], mszlog / m_tot[i], dne, dt, dp, smx / m_tot[i], smsz / m_tot[i], mcorr / m_tot[i], mall / m_tot[i]
            return mxlog, mszlog
            # return mxlog,mszlog,dne,dt,dp,smx,smsz

        bx = np.zeros(nbin)
        bsz = np.zeros(nbin)
        dne = np.zeros(nbin)
        dt = np.zeros(nbin)
        dp = np.zeros(nbin)
        smx = np.zeros(nbin)
        smsz = np.zeros(nbin)
        bcorr = np.zeros(nbin)
        ball = np.zeros(nbin)

        mx_out = np.zeros(nbin) #Outflowing material == sector 1
        msz_out = np.zeros(nbin)

        mx_col = np.zeros(nbin) #Collapsing material == sector 8
        msz_col = np.zeros(nbin)

        mx_full = np.zeros(nbin)
        msz_full = np.zeros(nbin)

        mx_fil = np.zeros(nbin) #Filament material == sector 5
        msz_fil = np.zeros(nbin)

        alpha = np.load('alpha.npy')

        print("Mass calculation")
        for i in range(nbin):
            #print("i",i)
            # bx[i], bsz[i], dne[i], dt[i], dp[i], smx[i], smsz[i], bcorr[i], ball[i] = compar_bias(i, plog, tlog, ne_log, m_tot,alpha, sigma_r2_log, sigma_t2_log, vt_log)
            mx_out[i], msz_out[i] = compar_bias(i, p1, t1, ne1)
            mx_col[i], msz_col[i] = compar_bias(i, p8, t8, ne8)
            mx_full[i], msz_full[i] = compar_bias(i, p_all, t_all, ne_all)
            mx_fil[i], msz_fil[i] = compar_bias(i, p5, t5, ne5)
        #print("msz_out",msz_out)
        #print("msz_col",msz_col)
        #print("msz_full",msz_full)
        #print("msz_fil",msz_fil)
        #sys.exit()
        plt.plot(rlogbin_cen,msz_out/m_cumul_sum, label="Outflowing material",color='orange',marker='.')
        plt.plot(rlogbin_cen,msz_col/m_cumul_sum, label="Collapsing material",color='dodgerblue',marker='.')
        plt.plot(rlogbin_cen,msz_full/m_cumul_sum, label="Full material",color='black',marker='.')
        plt.plot(rlogbin_cen,msz_fil/m_cumul_sum, label="Filament material",color='mediumvioletred',marker='.')
        plt.xlim(200,4000)
        plt.xscale('log')
        plt.ylim(0.5,3.5)
        plt.axvline(x=1087, color='grey',ls='dashed')
        plt.text(1087, 0.6, "$R_{500}$", rotation=90, size=16)
        plt.axvline(x=2147, color='grey', ls='dashed')
        plt.text(2147, 0.6, "$R_{vir}$", rotation=90, size=16)
        plt.axhline(1,color='black',ls='dashed')
        plt.xlabel('R [kpc]',size=16)
        plt.ylabel('$M_{HE}/M_{tot}$',size=16)
        #plt.yscale('log')
        plt.legend()
        plt.show()

    #plot_pressure()
    #plot_density()
    #plot_bias()
    #plot_temperature()

    #plot_entropy()

    #triple_plot()
    #double_plot()
    #triple_plot()
    #plot_pressure_grad()
    #plot_entropy_grad()

    #plot_density_grad()
    #compare_gradient()

    bias_w_prop_from_sectors_77()




rad_press_prof_sectors()

