import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
import matplotlib.lines as mlines
#import proplot as pplt
from scipy import signal
from astropy.cosmology import Planck18 as cosmo
#import pickle
from mpl_toolkits import mplot3d
#from ipywidgets import interact
from matplotlib.colors import Normalize
import sys


from astropy import constants as const
import matplotlib.colors as colors
kb = const.k_B.to_value('J/K')

from matplotlib.patches import Rectangle
from matplotlib.ticker import LogLocator, LogFormatter

import matplotlib.gridspec as gridspec

#import seaborn as sns

from scipy.stats import skew
from scipy.stats import kurtosis
from scipy.stats import norm

from lmfit.models import GaussianModel

import f_to_py as ftp

unit_kpc = 2.2755E27 / 3.08567758128E21
unit_l = 0.227550518265197E+28
unit_d = 0.264088847033921E-29
unit_t = 0.455724529603677E+18
unit_p = unit_d * (unit_l / unit_t) ** 2

kb = const.k_B.to_value('J/K')
mp = const.m_p.to_value('kg')
me = const.m_e.to_value('kg')
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
m_sun = const.M_sun.to_value('kg')
kev = 1.602176634e-16

def mapf90(type, lvl, proj,file,titre,filament,savecond,savefile):

    #plt.figure(facecolor='white')

    def map_bin_to_npy():
        #ncell = 888 ** 3
        ncell = 983 ** 3
        #ncell = 1966 ** 3

        #file = "/data/cluster/tlebeau/virgo/maps/3D_pyvista_lvl16_with_gals_P_fil.bin"
        #file = "/data/cluster/tlebeau/virgo/maps/3D_pyvista_lvl15_with_gals_2.bin"
        file = "/data/cluster/tlebeau/virgo/maps/3D_pyvista_lvl15_without_gals_P_fil_los.bin"

        #map = np.zeros(1966**3)

        map = ftp.f90_to_py.read_map_file(ncell,file)

        print("map",map)
        print("successfull data reading")
        #sys.exit()

        np.save("/data/cluster/tlebeau/virgo/maps/3D_pyvista_lvl15_without_gals_P_fil_los.npy",map)

        print("file saved in npy format")

        sys.exit()

    #map_bin_to_npy()

    #map = np.loadtxt(file)
    h = FortranFile(file, 'r')

    if filament==1:
        nx,ny,nz = h.read_ints()
        cen_x,cen_y,cen_z = h.read_reals()
        print("nx,ny,nz",nx,ny,nz)
        #print("cen_x,cen_y,cen_z",cen_x,cen_y,cen_z)
        #sys.exit()
        #ncell = 0
        #print(ncell.dtype)
        #sys.exit()
        #ncell = np
        #sys.exit()
        if proj == "x":
            ncell = nz * ny
        elif proj == "y":
            ncell = nx * nz
        elif proj == "z":
            #ncell = np.uint64(nx * ny)
            ncell = nx * ny
        print("ncell",ncell)

        map = np.zeros(ncell)

        map = ftp.f90_to_py.read_map_file(ncell, file,0)
        #map=h.read_reals()

        if type=="vnon":
            #f = FortranFile('./maps/high_res/filament/map_high_19_xyz_left_map_vy_z_d0.08.bin', 'r')

            #f = FortranFile('./maps/high_res/filament/map_high_19_xyz_full_map_vy_z_d0.02.bin', 'r')

            #filetwo = './maps/high_res/filament/map_high_19_xyz_left_minus_0Mpc_map_vy_x_d0.02.bin'

            filetwo = './maps/high_res/velocity/map_high_21_z_pm_500_kpc_map_vy_d0.02.bin'

            f = FortranFile(filetwo, 'r')

            nx, ny, nz = f.read_ints()

            if proj == "x":
                ncell = nz * ny
            elif proj == "y":
                ncell = nx * nz
            elif proj == "z":
                ncell = nx * ny
            print("ncell",ncell)

            map2 = np.zeros(ncell)

            #map2 = ftp.f90_to_py.read_map_file(ncell, './maps/high_res/filament/map_high_19_xyz_full_map_vy_z_d0.02.bin',0)

            map2 = ftp.f90_to_py.read_map_file(ncell,filetwo, 0)



            if proj == "x":
                map2 = np.reshape(map2, (nz, ny))
            elif proj == "y":
                map2 = np.reshape(map2, (nx, nz))
            elif proj == "z":
                map2 = np.reshape(map2, (ny, nx))


            #g = FortranFile( './maps/high_res/filament/map_high_19_xyz_full_map_vz_z_d0.08.bin', 'r')
            #nx, ny, nz = g.read_ints()

            #if proj == "x":
            #    ncell = nz * ny
            #elif proj == "y":
            #    ncell = nx * nz
            #elif proj == "z":
            #    ncell = nx * ny
            #print("ncell", ncell)

            #map3 = np.zeros(ncell)

            #map3 = ftp.f90_to_py.read_map_file(ncell,'./maps/high_res/filament/map_high_19_xyz_full_map_vz_z_d0.08.bin',0)

            #map = np.abs(map-map2)

            #map[np.isnan(map) == True] = 1e15
            #map[map==0] = 1e-15

            #print("max diff",np.max(map),"min diff",np.min(map))

        if type == "f":
            # f = FortranFile('./maps/high_res/filament/map_high_19_xyz_left_map_vy_z_d0.08.bin', 'r')

            # f = FortranFile('./maps/high_res/filament/map_high_21_xyz_left_map_vx_z_d0.02.bin', 'r')

            f = FortranFile('./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_fy_x_d0.02.bin', 'r')

            nx, ny, nz = f.read_ints()

            if proj == "x":
                ncell = nz * ny
            elif proj == "y":
                ncell = nx * nz
            elif proj == "z":
                ncell = nx * ny
            print("ncell", ncell)

            map2 = np.zeros(ncell)

                # map2 = ftp.f90_to_py.read_map_file(ncell, './maps/high_res/filament/map_high_19_xyz_left_map_vy_z_d0.08.bin',0)

            map2 = ftp.f90_to_py.read_map_file(ncell,'./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_fy_x_d0.02.bin',0)

            if proj == "x":
                map2 = np.reshape(map2, (nz, ny))
            elif proj == "y":
                map2 = np.reshape(map2, (nx, nz))
            elif proj == "z":
                map2 = np.reshape(map2, (ny, nx))


        print("len nan",len(map[np.isnan(map)==True]),'ratio, ',len(map[np.isnan(map)==True])/len(map))
        #sys.exit()

        if proj == "x":
            map = np.reshape(map, (nz, ny))
            #map2 = np.reshape(map2, (nz, ny))
            #map3 = np.reshape(map3, (nz, ny))
            nl = nx
        elif proj == "y":
            map = np.reshape(map, (nx, nz))
            #map2 = np.reshape(map2, (nx, nz))
            #map3 = np.reshape(map3, (nx, nz))
            nl = ny
        elif proj == "z":
            map = np.reshape(map, (ny, nx))
            #map2 = np.reshape(map2, (ny, nx))
            #map3 = np.reshape(map3, (ny, nx))
            nl = nz

    elif (filament==2):

        nx, ny, nz = h.read_ints()
        cen_x, cen_y, cen_z = h.read_reals()

        nl = nx

        #nx = 15
        #ny = 15728
        #nz = 15728

        file = './maps/high_res/filament/map_high_19_xyz_left_minus_0Mpc_map_vzcomp_x_d0.02.bin.npy'

        map = np.load(file)

        #print('shape map',map.shape)



        filetwo = './maps/high_res/filament/map_high_19_xyz_left_minus_0Mpc_map_vycomp_x_d0.02.bin.npy'

        map2 = np.load(filetwo)

        #print('shape map2',map2.shape)

        #sys.exit()

    else :
        map = h.read_reals()
    #map=[]
    #nline=h.read_ints()
    #print("nline",nline)
    #nline = 16
    #sys.exit()
    #for i in range(nline):
        #print(i)
        #m=h.read_reals()
        #print("m",m)
        #sys.exit()
        #map=np.concatenate((map,m))
    #map=np.array(map)
        print("len map",len(map))
        print("map",map)
    #sys.exit()
        px=np.int_(np.sqrt(len(map)))
        print("nbr px per line",px)
        map = np.reshape(map, (px, px))

        nx=px




    #sys.exit()
    #plt.hist(map,bins=100)
    #plt.show()

    print('file loaded')
    #sys.exit()

    #for i in range(len(map[0, :])):
    #    for j in range(len(map[:, 0])):
    #        if (map[i, j] == 0):
    #            map[i, j] = np.nan

    #map2 = np.loadtxt("map_2000px_z_ne_los_cleanall+T1e5_highr.txt")
    #for i in range(len(map2[0, :])):
    #    for j in range(len(map2[:, 0])):
            #if (map2[i, j] == 0):
            #    map2[i, j] = 1e-50

    #map-=map2
    print(map)
    print("min=",np.min(map))
    print("max=",np.max(map))
    #plt.figure()

    if (type=="SD"):

        file = './maps/high_res/map_high_19_cen_SD_DM_los.bin'
        #file = './maps/high_res/random_proj/sd_dm/map_high_19_rand_sum100_sd_dm_los.bin'
        h2 = FortranFile(file, 'r')
        map2 = h2.read_reals()
        #map=[]
        # nline=h.read_ints()
        # print("nline",nline)
        # nline = 16
        # sys.exit()
        # for i in range(nline):
        # print(i)
        # m=h.read_reals()
        # print("m",m)
        # sys.exit()
        # map=np.concatenate((map,m))
        # map=np.array(map)
        print("len map", len(map2))
        print("map", map2)
        # sys.exit()
        px = np.int_(np.sqrt(len(map2)))
        print("nbr px per line", px)
        # sys.exit()
        # plt.hist(map,bins=100)
        # plt.show()
        map2 = np.reshape(map2, (px, px))
        print('file 2 loaded')

        map = map + map2

    def dessins_gal():
        cen_x = 0.48461068
        cen_y = 0.50809848
        cen_z = 0.49687076

        cen_x = (cen_x - 0.5) * (unit_l / 3.08567758128E21)
        cen_y = (cen_y - 0.5) * (unit_l / 3.08567758128E21)
        cen_z = (cen_z - 0.5) * (unit_l / 3.08567758128E21)

        width = 5000
        minx = cen_x - width
        maxx = cen_x + width
        miny = cen_y - width
        maxy = cen_y + width
        minz = cen_z - width
        maxz = cen_z + width

        nx = 2000
        ny = 2000
        nz = 2000

        lenx = maxx - minx
        leny = maxy - miny
        lenz = maxz - minz

        dx = lenx / nx
        dy = leny / ny
        dz = lenz / nz

        gal_list = np.loadtxt("/data/cluster/tlebeau/virgo/list_gal_251.dat_js_nocontam")
        cond = np.logical_and(gal_list[:, 20] > 0.47783, np.logical_and(gal_list[:, 20] < 0.49139,
                                                                        np.logical_and(gal_list[:, 21] > 0.50131,
                                                                                       np.logical_and(
                                                                                           gal_list[:, 21] < 0.51487,
                                                                                           np.logical_and(gal_list[:,
                                                                                                          22] > 0.49009,
                                                                                                          gal_list[:,
                                                                                                          22] < 0.50365)))))
        print("cond", cond[0:40])
        x = (gal_list[:, 3] - 0.5) * (unit_l / 3.08567758128E21)
        y = (gal_list[:, 4] - 0.5) * (unit_l / 3.08567758128E21)
        x = (x - cen_x) / 1000
        y = (y - cen_y) / 1000
        r200dm = gal_list[:, 23] * (unit_l / 3.08567758128E21)
        mgal = gal_list[:, 2]

        print("len(x)", len(x))

        # print(x)
        # print(y)
        # plt.scatter(gal_pos[:,0][np.logical_and(r500>0,r500<500)],gal_pos[:,1][np.logical_and(r500>0,r500<500)],alpha=1,s=1,c='black',label='0<r500<500kpc')
        # plt.scatter(gal_pos[:,0][r500>500],gal_pos[:,1][r500>500],alpha=1,s=1,c='grey',label='r500>500kpc')
        # plt.scatter(gal_pos[:,0][r500==0],gal_pos[:,1][r500==0],alpha=1,s=1,c='white',label='r500=0kpc')
        # plt.scatter(gal_pos[:,0][np.logical_and(r200>0,r200<800)],gal_pos[:,1][np.logical_and(r200>0,r200<800)],alpha=1,s=1,c='black',label='0<r200<800kpc')
        # plt.scatter(gal_pos[:,0][r200>800],gal_pos[:,1][r200>800],alpha=1,s=1,c='grey',label='r200>800kpc')
        # plt.scatter(gal_pos[:,0][r200==0],gal_pos[:,1][r200==0],alpha=1,s=1,c='white',label='r200=0kpc')
        # plt.legend()

        def mw():
            mw=np.zeros(3)
            mw[0]=-(0.48461068-0.5)*(unit_l/3.08567758128E21)
            mw[1]=-(0.50809848-0.5)*(unit_l/3.08567758128E21)
            mw[2]=-(0.49687076-0.5)*(unit_l/3.08567758128E21)

            print("mw",mw)

            #mw_x=0.5-0.48461068
            #mw_y=0.5-0.50809848
            #mw_z=0.5-0.49687076

            phi=np.arctan(mw_x/mw_z)
            theta=np.arctan(mw_x/mw_y)

            print(phi)
            print(theta)

            # theta=0

            # c,s=np.cos(theta),np.sin(theta)
            # print("cos",c,"sin",s)
            # if(dir=="z")
            # rx=np.array(((1,0,0),(0,c,-s),(0,s,c)))
            # ry=np.array(((c,0,s),(0,1,0),(-s,0,c)))
            # rz=np.array(((c,-s,0),(s,c,0),(0,0,1)))
            # print("rz",rz)
            # print(rx)

            # coord=np.array((x,y,z))
            # print(coord)

            # coord2=np.matmul(mw.T,rz)
            # coord2=coord2.T

            # print("ccord2",coord2)

            # phi=0

            # c,s=np.cos(phi),np.sin(phi)
            # ry=np.array(((c,0,s),(0,1,0),(-s,0,c)))

            # coord3=np.matmul(coord2,ry)
            # coord3=coord3.T
            # print("coord3",coord3)

            # mw[0]=coord3[0]+(0.48461068-0.5)*(unit_l/3.08567758128E21)
            # mw[1]=coord3[1]+(0.50809848-0.5)*(unit_l/3.08567758128E21)
            # mw[2]=coord3[2]+(0.49687076-0.5)*(unit_l/3.08567758128E21)

            # print("mw",mw)

            # mw[0]=(mw[0]-cen_x)/1000
            # mw[1]=(mw[1]-cen_y)/1000
            # mw[2]=(mw[2]-cen_z)/1000

            # mw[:]=[11.34871563702,-5.97215124594,2.30762989658]
            # mw[:]=[-9.0949470177292824E-013,-12.824196543046277,2.3076327260135868]
            # mw[:]=[-9.0949470177292824E-013,-4.5474735088646412E-013,13.030164456861203]

            # print("mw",mw)

            # plt.scatter(x[9],y[9],alpha=0.9,s=1,c='white')
            # plt.scatter(mw[0],mw[1],alpha=0.9,s=10,c='white')

        xpx = np.linspace(-5, 5, 2000)
        ypx = np.linspace(-5, 5, 2000)

        Xpx, Ypx = np.meshgrid(xpx, ypx)
        for i in range(len(mgal)):

            if (mgal[i] > 1e9 and mgal[i] < 1e13):
                print(i)
                # F = (Xpx-gal_pos[i,0])**2 + (Ypx-gal_pos[i,1])**2 - (r200dm[i]*(10/1e4))**2
                F = (Xpx - x[i]) ** 2 + (Ypx - y[i]) ** 2 - (r200dm[i] * (10 / 1e4)) ** 2
                # F = Xpx**2 + Ypx**2 - (10**((35+i)*0.05)*(10/1e4))**2
                plt.contour(Xpx, Ypx, F, [0], colors='white', linewidths=0.6, alpha=0.7)
                plt.scatter(x[i], y[i], alpha=0.9, s=1, c='white')

    #dessins_gal()

    def show_radii():
        #xpx = np.linspace(-10,10,4000)
        #ypx = np.linspace(-10,10,4000)

        #xpx = np.linspace(-5,5,4000)
        #ypx = np.linspace(-5,5,4000)

        xpx = np.linspace(-2.5,2.5,4000)
        ypx = np.linspace(-2.5,2.5,4000)

        Xpx, Ypx = np.meshgrid(xpx, ypx)

        # F = (Xpx-gal_pos[i,0])**2 + (Ypx-gal_pos[i,1])**2 - (r200dm[i]*(10/1e4))**2
        #F = (Xpx - x[i]) ** 2 + (Ypx - y[i]) ** 2 - (r200dm[i] * (10 / 1e4)) ** 2

        if(type=="ne" or type=="P" or type=="dm" or type=="v"):

            F = Xpx**2 + Ypx**2 - (2024*(10/1e4))**2  #rvir
            rvir=plt.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9)
            plt.clabel(rvir,rvir.levels,inline=True,fmt="$R_{vir}$",fontsize=20)

            F = Xpx ** 2 + Ypx ** 2 - (1087 * (10 / 1e4)) ** 2  # r500
            r500 = plt.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9)#, linestyles='dashed')
            plt.clabel(r500, r500.levels, inline=True, fmt="$R_{500}$", fontsize=20)#, manual=pos)

            #F = Xpx ** 2 + Ypx ** 2 - (5031 * (10 / 1e4)) ** 2  #rzv
            #rzv=plt.contour(Xpx, Ypx, F, [0], colors='white',ls='dashed', linewidths=2, alpha=0.9)
            #plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{zv}$", fontsize=20)

            #F = Xpx ** 2 + Ypx ** 2 - (3516 * (10 / 1e4)) ** 2  # 3.516 Mpc
            #rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', ls='dashed', linewidths=2, alpha=0.9)
            #plt.clabel(rzv, rzv.levels, inline=True, fmt="$3.516 Mpc$", fontsize=20)

            #F = Xpx ** 2 + Ypx ** 2 - (3162 * (10 / 1e4)) ** 2  # 3.162 Mpc
            #rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', ls='dashed', linewidths=2, alpha=0.9)
            #plt.clabel(rzv, rzv.levels, inline=True, fmt="$3.162 Mpc$", fontsize=20)

            #F = Xpx ** 2 + Ypx ** 2 - (5000 * (10 / 1e4)) ** 2  # 5 Mpc
            #rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', ls='dashed', linewidths=2, alpha=0.9)
            #plt.clabel(rzv, rzv.levels, inline=True, fmt="$5 Mpc$", fontsize=20)

        if (type=="y" or type=="SD" or type=="EM"):

            F = Xpx ** 2 + Ypx ** 2 - (2895 * (10 / 1e4)) ** 2  # r200m
            rvir = plt.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9)
            plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{200m}$", fontsize=20)

        ##Rsp from baryons 3D rad prof

        #F = Xpx ** 2 + Ypx ** 2 - (4890 * (10 / 1e4)) ** 2  # rsb_relaxed bar
        #rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dotted', linewidths=2, alpha=0.9)
        #plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,rel}$", fontsize=20)

        #F = Xpx ** 2 + Ypx ** 2 - (3890 * (10 / 1e4)) ** 2  # rsb_spherical_collapse bar
        #rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
        #plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,sph}$", fontsize=20)

        ##Rsp from DM 3D rad prof

        #F = Xpx ** 2 + Ypx ** 2 - (4300 * (10 / 1e4)) ** 2  # rsb_relaxed DM
        #rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', ls='dashed', linewidths=2, alpha=0.9)
        #plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,rel}$", fontsize=20)

        #F = Xpx ** 2 + Ypx ** 2 - (3400 * (10 / 1e4)) ** 2  # rsb_spherical_collapse DM
        #rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', ls='dotted', linewidths=2, alpha=0.9)
        #plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,sph}$", fontsize=20)

        if (type=="EM"): ##Rsp from EM

            F = Xpx ** 2 + Ypx ** 2 - (3900 * (10 / 1e4)) ** 2
            rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,EM}$", fontsize=20)

            #F = Xpx ** 2 + Ypx ** 2 - (3600 * (10 / 1e4)) ** 2
            #rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            #plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,EM}$", fontsize=20)

            #F = Xpx ** 2 + Ypx ** 2 - (1098 * (10 / 1e4)) ** 2
            #rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            #plt.clabel(rzv, rzv.levels, inline=True, fmt="$test$", fontsize=20)

        if (type=="y"): ##Rsp from y_c

            F = Xpx ** 2 + Ypx ** 2 - (4900 * (10 / 1e4)) ** 2
            rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,y_c}$", fontsize=20)

        if (type=="SD"): ##Rsp from Surface density

            F = Xpx ** 2 + Ypx ** 2 - (3300 * (10 / 1e4)) ** 2
            rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,SD}$", fontsize=20)#,manual=True)

        if(type=="vproj" or type=="d"):
            #xpx = np.linspace(-2.5, 2.5, 4000)
            #ypx = np.linspace(-2.5, 2.5, 4000)

            xpx = np.linspace(-2.5, 2.5, 4000)
            ypx = np.linspace(-2.5, 2.5, 4000)

            Xpx, Ypx = np.meshgrid(xpx, ypx)

            F = Xpx**2 + Ypx**2 - (1087*(10/1e4))**2  #r500
            r500=plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=2, alpha=0.9,linestyles='dashed')
            pos=[(-0.75,2)]
            plt.clabel(r500,r500.levels,inline=True,fmt="$R_{500}$",fontsize=20,manual=pos)
            F = Xpx ** 2 + Ypx ** 2 - (2147 * (10 / 1e4)) ** 2  # rvir
            rvir = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=2, alpha=0.9,linestyles='solid')
            plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{vir}$", fontsize=20,manual=pos)

            F = Xpx ** 2 + Ypx ** 2 - (515 * (10 / 1e4)) ** 2  # rvir
            r2500 = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=2, alpha=0.9,linestyles='dotted')
            plt.clabel(r2500, r2500.levels, inline=True, fmt="$R_{2500}$", fontsize=20,manual=pos)

            #F = Xpx ** 2 + Ypx ** 2 - (100 * (10 / 1e4)) ** 2  # rvir
            #rvir = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=2, alpha=0.9,linestyles='solid')
            #plt.clabel(rvir, rvir.levels, inline=True, fmt="$100\,kpc$", fontsize=20)

        if (type == "mach" or type=='T'):
            # xpx = np.linspace(-2.5, 2.5, 4000)
            # ypx = np.linspace(-2.5, 2.5, 4000)

            xpx = np.linspace(-2.5, 2.5, 4000)
            ypx = np.linspace(-2.5, 2.5, 4000)

            Xpx, Ypx = np.meshgrid(xpx, ypx)

            F = Xpx ** 2 + Ypx ** 2 - (1087 * (10 / 1e4)) ** 2  # r500
            r500 = plt.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9, linestyles='dashed')
            pos = [(-0.75, 2)]
            plt.clabel(r500, r500.levels, inline=True, fmt="$R_{500}$", fontsize=20, manual=pos)

            F = Xpx ** 2 + Ypx ** 2 - (2147 * (10 / 1e4)) ** 2  # rvir
            rvir = plt.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9, linestyles='solid')
            plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{vir}$", fontsize=20, manual=pos)

            F = Xpx ** 2 + Ypx ** 2 - (515 * (10 / 1e4)) ** 2  # rvir
            r2500 = plt.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9,linestyles='dotted')
            plt.clabel(r2500, r2500.levels, inline=True, fmt="$R_{2500}$", fontsize=20,manual=pos)


            # F = Xpx ** 2 + Ypx ** 2 - (500 * (10 / 1e4)) ** 2  # rvir
            # rvir = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=2, alpha=0.9,linestyles='dotted')
            # plt.clabel(rvir, rvir.levels, inline=True, fmt="$500\,kpc$", fontsize=10)

            # F = Xpx ** 2 + Ypx ** 2 - (100 * (10 / 1e4)) ** 2  # rvir
            # rvir = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=2, alpha=0.9,linestyles='solid')
            # plt.clabel(rvir, rvir.levels, inline=True, fmt="$100\,kpc$", fontsize=20)

        #F = Xpx ** 2 + Ypx ** 2 - (8000 * (10 / 1e4)) ** 2  # background lower limit
        #rbgl = plt.contour(Xpx, Ypx, F, [0], colors='grey', linewidths=0.6, alpha=0.9)
        #plt.clabel(rbgl, rbgl.levels, inline=True, fmt="$Background lower limit}$", fontsize=10)

        #F = Xpx ** 2 + Ypx ** 2 - (10000 * (10 / 1e4)) ** 2  # background upper limit
        #rbgh = plt.contour(Xpx, Ypx, F, [0], colors='grey', linewidths=0.6, alpha=0.9)
        #plt.clabel(rbgh, rbgh.levels, inline=True, fmt="$Background upper limit}$", fontsize=10)

        #plt.fill(rbgl,rbgh)

        #F = Xpx ** 2 + Ypx ** 2 - (1648 * (10 / 1e4)) ** 2 #r200
        #r200=plt.contour(Xpx, Ypx, F, [0], colors='blue', linewidths=0.6, alpha=0.9)
        #plt.clabel(r200, r200.levels, inline=True, fmt="$R_{200}$", fontsize=10)
        #F = Xpx ** 2 + Ypx ** 2 - (1085 * (10 / 1e4)) ** 2 #r500
        #r500=plt.contour(Xpx, Ypx, F, [0], colors='blue',linestyles="dashed", linewidths=0.6, alpha=0.9)
        #plt.clabel(r500, r500.levels, inline=True, fmt="$R_{500}$", fontsize=10)
        #F = Xpx ** 2 + Ypx ** 2 - (900 * (10 / 1e4)) ** 2  # bump radius
        #rbump = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
        #plt.clabel(rbump, rbump.levels, inline=True, fmt="$R_{bump}$", fontsize=10)
        #F = Xpx ** 2 + Ypx ** 2 - (850 * (10 / 1e4)) ** 2  # bump radius
        #rbump = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
        #plt.clabel(rbump, rbump.levels, inline=True, fmt="$R_{bump}$", fontsize=10)

        #F = Xpx ** 2 + Ypx ** 2 - (400 * (10 / 1e4)) ** 2  # low mw deproj radius
        #rdec = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
        #plt.clabel(rdec, rdec.levels, inline=True, fmt="$R_{decrease in deproj prof}$", fontsize=10)
        #plt.scatter(x[i], y[i], alpha=0.9, s=1, c='white')
        #plt.show()
        #sys.exit()

    def show_quadrants():

        if(type=="vproj" or type=="d"):
            #xpx = np.linspace(-2.5, 2.5, 4000)
            #ypx = np.linspace(-2.5, 2.5, 4000)

            xpx = np.linspace(-2.5, 2.5, 4000)
            ypx = np.linspace(-2.5, 2.5, 4000)

            rvir=2.147
            xrvir=rvir*0.5*np.sqrt(2)
            xpx_vir = np.linspace(-xrvir, xrvir, 4000)

            #Xpx, Ypx = np.meshgrid(xpx, ypx)

            #F = Xpx**2 + Ypx**2 - (1087*(10/1e4))**2  #r500
            #r500=plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=2, alpha=0.9,linestyles='dashed')
            #plt.clabel(r500,r500.levels,inline=True,fmt="$R_{500}$",fontsize=20)

            #F = Xpx ** 2 + Ypx ** 2 - (2147 * (10 / 1e4)) ** 2  # rvir
            #rvir = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=2, alpha=0.9,linestyles='dashdot')
            #plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{vir}$", fontsize=20)

            plt.plot(xpx_vir,xpx_vir,color='black',lw=2,alpha=0.9)
            plt.plot(xpx_vir,-xpx_vir,color='black',lw=2,alpha=0.9)
            plt.axvline(x=0,ymin=(2.5-rvir)/5,ymax=(2.5+rvir)/5,color='black',lw=2,alpha=0.9)
            plt.axhline(y=0,xmin=(2.5-rvir)/5,xmax=(2.5+rvir)/5,color='black',lw=2,alpha=0.9)



    #show_radii()

    x=[-10,10]
    y1=[10,10]
    y2=[-10,-10]
    #plt.fill_between(x,y1,y2,color='none',hatch='/',edgecolor="black")

    if filament==1 or filament==2:
        if proj == "x":
            dimx = (ny / 2) * (737.441 / 2 ** lvl)
            dimy = (nz / 2) * (737.441 / 2 ** lvl)
        elif proj == "y":
            dimx = (nz / 2) * (737.441 / 2 ** lvl)
            dimy = (nx / 2) * (737.441 / 2 ** lvl)
        elif proj == "z":
            map = np.reshape(map, (nx, ny))
            dimx = (nx / 2) * (737.441 / 2 ** lvl)
            dimy = (ny / 2) * (737.441 / 2 ** lvl)

    else :
        dimx=(nx/2)*(737.441/2**lvl)
        dimy=dimx


    dim = [-dimx, dimx, -dimy, dimy]
    print("dim",dim)
    #sys.exit()
    if (type == "ne"):
        #map[map == 0] = 1
        #map /= np.mean(map)
        #map = np.log10(map)
        #map[np.isnan(map)==True] = 16


        fig, ax = plt.subplots(figsize=(8, 8),facecolor='white')
        im=ax.imshow(map, cmap="Solar", origin='lower', alpha=1, extent=dim,norm=colors.LogNorm(vmin=10**(16), vmax=10**(20)))#,norm=colors.LogNorm(vmin=10**(-7), vmax=10**(-2.5)))#, vmin=18,vmax=20)#,vmin=14,vmax=19)#,vmin=16)#, vmin=16, vmax=21)  # 3D: 0 6.5 #filament : vmin=14, vmax=19
        show_radii()
        #ax.axis('off')

        #divider = make_axes_locatable(ax)
        #cax = divider.new_vertical(size='5%', pad=0.6, pack_start=True)
        #fig.add_axes(cax)
        #cb = fig.colorbar(im, cax=cax, orientation='horizontal')

        #cb.set_label('$log_{10}(n_e[cm^{-2}])$',size='large')

        ax.set_xlabel("x [Mpc]")
        ax.set_ylabel("y [Mpc]")

        ratio = map.shape[1] / map.shape[0]

        cb = fig.colorbar(im, ax=ax, orientation="horizontal",fraction=0.047 * ratio)  # , label='$T~\mathrm{[K]}$',size="large")#,pad=0.1)
        cb.ax.tick_params(labelsize=14)
        #cb=fig.colorbar(im)
        cb.set_label('$n_e~\mathrm{[cm^{-3}]}$',size='large')
        #cb.set_label('$log_{10}(n_e/\overline{n_e})$', size='large')

    elif (type == "T"):
        fig, ax = plt.subplots(figsize=(8, 8),facecolor='white')
        ##test in K :
        #map /= (kb / 1.602e-16)
        #map=np.log10(map)

        #im=ax.imshow(map, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e6, vmax=1e8))#, vmin=0, vmax=6)  "hot"

        im = ax.imshow(map, cmap="magma", origin='lower', alpha=1, extent=dim)#,norm=colors.Normalize(vmin=1.7, vmax=12))  # , vmin=0, vmax=6)  "hot"
        show_radii()
        #ax.axis('off')


        ratio = map.shape[1]/map.shape[0]

        #

        size=22
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction= 0.047*ratio)#, label='$T~\mathrm{[K]}$',size="large")#,pad=0.1)
        cb.ax.tick_params(labelsize=14)
        cb.set_label('$T~\mathrm{[keV]}$', size=size)
        cb.ax.tick_params(labelsize=size)
        ax.set_xlabel("Mpc", size=size)
        ax.set_ylabel("Mpc", size=size)
        plt.yticks(fontsize=size)
        plt.xticks(fontsize=size)
        #divider = make_axes_locatable(ax)
        #cax = divider.new_vertical(size='5%', pad=0.6, pack_start=True)
        #fig.add_axes(cax)
        #cb = fig.colorbar(im, cax=cax, orientation='horizontal')

        #cb.set_label('$T[keV]$',size='large')
        #cb.set_label('T[K]', size='large')

        #ax.get_xaxis().set_visible(False)
        #ax.get_yaxis().set_visible(False)

    elif (type == "P"):
        #map[np.isnan(map) == True] = 1e-15
        #map = np.log10(map)
        #print("test")

        mid = int(nx / 2)

        px_size = 737.441 / 2 ** lvl

        npx500 = int(0.5 / px_size)

        min = mid - npx500
        max = mid + npx500

        map_small = map[min:max, min:max]
        dim_small = [-0.5, 0.5, -0.5, 0.5]


        fig,ax= plt.subplots(facecolor="white")
        #im=ax.imshow(map_small, cmap="inferno", origin='lower', alpha=1, extent=dim_small, norm=colors.LogNorm(vmin=1e-4,vmax=4e-3) )#,vmin=-8,vmax=0)#,vmin=-10)#, vmin=-8, vmax=-1.5) #inferno #for pm 1Mpc maps : 4e-6
        im = ax.imshow(map, cmap="inferno", origin='lower', alpha=1, extent=dim,norm=colors.LogNorm(vmin=1e-8, vmax=3e-2))
        show_radii()

        size=22

        #if filament==1:
        #    im = ax.imshow(map, cmap="inferno", origin='lower', alpha=1, extent=dim, vmin=-8, vmax=-1.5) #, extent=dim)#, vmin=-8, vmax=-1.5)
            #im = ax.imshow(map, cmap="binary", origin='lower', alpha=1, extent=dim)

        #show_radii()
        #ax.axis('off')

        #divider = make_axes_locatable(ax)
        #cax = divider.new_vertical(size='5%', pad=0.6, pack_start=True)
        #fig.add_axes(cax)
        #cb=fig.colorbar(im, cax=cax, orientation='horizontal')
        #cb.ax.tick_params(labelsize=22)
        cb=fig.colorbar(im,ax=ax)
        ax.set_xlabel("x [Mpc]")
        ax.set_ylabel("y [Mpc]")

        #cb=plt.colorbar(im,ax=[ax], location='bottom')
        cb.set_label('$P~[\mathrm{keV~cm^{-3}}]$',size='large')

    elif (type == "dm"):
        plt.imshow(map, cmap="cividis", origin='lower', alpha=1, extent=dim)

    elif (type == "y"):
        def gaussian(x, mean, sigma):
            return (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-((x - mean) ** 2) / (2 * sigma ** 2))

        #x = np.linspace(0, 15728, 15728)
        #FWHM=14
        #filter = gaussian(x, 7864, FWHM)
        #filter /= np.trapz(filter)
        # print(np.trapz(filter))
        #plt.plot(x,filter)
        #plt.show()

        #kernel = filter[:, np.newaxis] * filter[np.newaxis, :]

        #print("convolution")

        #map2 = signal.fftconvolve(map,kernel,mode='same')
        #map2 *= 1.602e-10

        #print("convulotion done")

        map[map == 0] = 1e-12
        #map = np.log10(map)
        map[np.isnan(map) == True] = 1e-12
        map = np.log10(map)

        #sys.exit()

        #map *=1.602e-10

        #map = np.log10(map)
        fig, ax = plt.subplots()
        im = ax.imshow(map, cmap="viridis", origin='lower', alpha=1, extent=dim, vmin=-10)#, vmax=6e-6)
        #im = ax.imshow(kernel, cmap="Greys", origin='lower')
        show_radii()
        # ax.axis('off')

        # divider = make_axes_locatable(ax)
        # cax = divider.new_vertical(size='5%', pad=0.6, pack_start=True)
        # fig.add_axes(cax)
        # cb=fig.colorbar(im, cax=cax, orientation='horizontal')
        cb = fig.colorbar(im, ax=ax,orientation='horizontal')#,size=20)
        cb.ax.tick_params(labelsize=20)
        cb.set_label('$log_{10}(y)$', size=30)
        ax.axis('off')
        #ax.set_xlabel("x [Mpc]")
        #ax.set_ylabel("y [Mpc]")

        # cb=plt.colorbar(im,ax=[ax], location='bottom')
        #cb.set_label('$log_{10}(P[keV.cm^{-3}])$', size='large')

    elif (type == "v"):

        vx_virgo = -509.1301
        vy_virgo = 228.9488
        vz_virgo = -131.9249

        map-=vx_virgo

        if titre == "vx":
            map -= vx_virgo
            #map2 -= vy_virgo
        elif titre == "vy":
            map -= vy_virgo
        elif titre == "vz":
            map -= vz_virgo

        fig, ax = plt.subplots(facecolor='white', figsize=(8, 8))

        nbin = 10

        def vel_in_polar_coord(map,map2):

            v_x = map
            v_y = map2


            # Define the origin (center of the map)
            x_0 = v_x.shape[1] // 2
            y_0 = v_y.shape[0] // 2

            # Create coordinate grids
            y, x = np.indices(v_x.shape)

            # Calculate relative coordinates from the origin
            x_rel = x - x_0
            y_rel = y - y_0

            # Compute the angle theta for each point
            theta = np.arctan2(y_rel, x_rel)

            # Compute v_r and v_theta
            v_r = v_x * np.cos(theta) + v_y * np.sin(theta)
            v_theta = -v_x * np.sin(theta) + v_y * np.cos(theta)

            # v_r and v_theta are now the radial and angular components of the velocity

            return v_r, v_theta
        #print("converting in polar coord")
        #v_r,v_theta = vel_in_polar_coord(map,map2)
        #print("conversion done")



        #plt.hist(map.flatten(), bins=100, alpha=0.6, color='blue', label='vz')
        #plt.show()
        #sys.exit()

        #x, y = np.meshgrid(np.linspace(-dimx, dimx, int(nx/nbin)), np.linspace(-dimy, dimy, int(ny/nbin)))

        #np.save("mesh_x_pm_500_kpc",x)
        #np.save("mesh_y_pm_500_kpc",y)

        x = np.load("mesh_x.npy")
        y = np.load("mesh_y.npy")

        #x = np.load("mesh_x_pm_500_kpc.npy")
        #y = np.load("mesh_y_pm_500_kpc.npy")

        #print(np.shape(x))
        #print(np.shape(y))

        #sys.exit()

        #x = np.load("mesh_x_left_fil.npy")
        #y = np.load("mesh_y_left_fil.npy")

        print("mean calculation")

        #map -= vz_virgo
        #map2 -= vy_virgo
        #map3 -= vz_virgo

        #map_mean = np.sqrt(map**2 + map2**2)# + map3**2)

        #np.save("mean_v_2D_map.npy",map_mean)

        #map_mean = np.load("mean_v_map.npy")

        #print("max=",np.max(map_mean))
        #print("min=",np.min(map_mean))

        #print(np.shape(map))
        #print(np.shape(map2))

        #mean10_map = np.array([[np.mean(map[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for j in range(int(ny / nbin))] for i in range(int(nx / nbin))])

        #mean10_map = np.array([[np.nanmean(map[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for j in range(int(ny / nbin))] for i in range(int(nz / nbin))])

        #print("mean10_map done, shape: ",np.shape(mean10_map))

        #np.save("stream_vz_x_proj_pm_500_kpc",mean10_map)

        #mean10_map2 = np.array([[np.mean(map2[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for j in range(int(ny / nbin))] for i in range(int(nx / nbin))])

        #mean10_map2 = np.array([[np.nanmean(map2[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for j in range(int(ny / nbin))] for i in range(int(nz / nbin))])

        #print("mean10_map2 done, shape: ",np.shape(mean10_map2))

        #np.save("stream_vy_x_proj_pm_500_kpc",mean10_map2)



        mean10_map = np.load("stream_vx.npy")
        mean10_map2 = np.load("stream_vy.npy")

        #mean10_map = np.load("stream_vx_left_fil.npy")
        #mean10_map2 = np.load("stream_vy_left_fil.npy")

        #mean10_map = np.load("stream_vz_left_fil_minus_2Mpc.npy")
        #mean10_map2 = np.load("stream_vy_left_fil_minus_2Mpc.npy")

        #mean10_map = np.load("stream_vx_full_map.npy")
        #mean10_map2 = np.load("stream_vy_full_map.npy")

        #mean10_map = np.load("stream_vz_right_fil_plus_6Mpc.npy")
        #mean10_map2 = np.load("stream_vy_right_fil_plus_6Mpc.npy")


        #print("max mean map vx =",np.max(mean10_map))
        #print("min mean map vx =",np.min(mean10_map))
        #print("max mean map vy =",np.max(mean10_map2))
        #print("min mean map vy =",np.min(mean10_map2))

        #sys.exit()

        #mean10_map3 = np.array([[np.mean(map3[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for i in range(int(nx / nbin))] for j in range(int(ny / nbin))])

        #d = np.array([[np.log10(np.sqrt((cen_x - (i + 0.5)) ** 2 + (cen_y - (j + 0.5)) ** 2) * 5.6262) for i in range(max_x)] for j in range(max_y)])

        #print("mean calculated")

        print("mean mesh calculation")

        #print(np.shape(mean10_map))
        #print(np.shape(mean10_map2))

        #print("len nan", len(mean10_map[np.isnan(mean10_map) == True]), 'ratio, ', len(mean10_map[np.isnan(mean10_map) == True]) / len(mean10_map))
        #print("len nan", len(mean10_map2[np.isnan(mean10_map2) == True]), 'ratio, ', len(mean10_map2[np.isnan(mean10_map2) == True]) / len(mean10_map2))

        map_mean_mesh = np.sqrt(mean10_map**2 + mean10_map2**2)

        #print("mean mesh done, shape: ",np.shape(map_mean_mesh))

        #print("len nan", len(map_mean_mesh[np.isnan(map_mean_mesh) == True]), 'ratio, ', len(map_mean_mesh[np.isnan(map_mean_mesh) == True]) / len(map_mean_mesh))

        map_mean_mesh[np.isnan(map_mean_mesh) == True] = 0

        #plt.hist(map.flatten(),bins=100,range=(-2000,2000),histtype='step',color='black',linewidth=2)
        #plt.show()
        #sys.exit()

        #plt.hist(map_mean_mesh.flatten(),bins=100,range=(0,2000),histtype='step',color='black',linewidth=2)
        #plt.show()
        #sys.exit()

        #map_mean_mesh = np.sqrt(mean10_map[np.isnan(mean10_map) == False] ** 2 + mean10_map2[np.isnan(mean10_map2) == False] ** 2)

        #print(np.shape(map_mean_mesh))

        #map_diff = np.abs(map-map2)

        #print("max mean map =",np.max(map_mean_mesh))
        #print("min mean map =",np.min(map_mean_mesh))

        #norm = Normalize(vmin = np.min(map_mean_mesh), vmax = np.max(map_mean_mesh))

        #print("max mean map =", np.max(mean10_map2))
        #print("min mean map =", np.min(mean10_map2))

        #norm = Normalize(vmin=0, vmax=700)

        norm = Normalize(vmin=0, vmax=1200) #longitudinal slice
        #norm = Normalize(vmin=0, vmax=800) #trasnverse slice

        #plt.hist(map_mean_mesh.flatten(), bins=100, range=(0, 2000), histtype='step', color='black', linewidth=2)
        #plt.show()

        #sys.exit()

        #lw = 2 * map_mean_mesh / np.max(map_mean_mesh)

        #print("diff",map-map2)

        #print("lw",np.shape(lw))
        #sys.exit()
        #plt.quiver(x, y, mean10_map2, mean10_map3)
        stream = ax.streamplot(x, y, mean10_map, mean10_map2, linewidth=1, color=np.sqrt(mean10_map**2+mean10_map2**2), cmap='Greys_r', density=[2,2],norm=norm)#,alpha=0.8) # color="black" #,vmin=-1600,vmax=1600) #color=mean10_map2, cmap="BR"
        #plt.show()

        #map -= np.mean(map)

        vlim = 1600 ##longitudinal slice
        #vlim = 800 ##transverse slice
        im = ax.imshow(map, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-vlim, vmax=vlim) #BR #origin='lower'

        #im = ax.imshow(map, cmap="cividis", origin="lower", alpha=1, extent=dim, vmin=0, vmax=1000)  # BR #origin='lower'

        #im = ax.imshow(np.abs(map), cmap="Greys", origin="lower", alpha=1, extent=dim)#, vmin=-vlim, vmax=vlim)

        #im = ax.imshow(map_diff, cmap="Greys", origin="lower", alpha=1, extent=dim)#, vmin=-1600, vmax=1600)  # BR #origin='lower'
        #im = ax.imshow(map_mean, cmap="rocket", origin="lower", alpha=1, extent=dim,vmin=0,vmax=600)#, vmin=-1600, vmax=1600)  # BR #origin='lower'

        #factor = 3556/15728
        #dim = [-dimx*factor, dimx*factor, -dimy*factor, dimy*factor]

        print("dim",dim)

        ###im = ax.imshow(map_mean[6085:9641,6085:9641], cmap="rocket", origin="lower", alpha=1, extent=dim, vmin=0,vmax=600)  # , vmin=-1600, vmax=1600)  # BR #origin='lower'
        #show_radii()

        ratioh = map.shape[1] / map.shape[0]

        cb = fig.colorbar(im, ax=ax,orientation='horizontal',fraction=0.047*ratioh)#pad=0.04, shrink=0.95 ) #orientation='horizontal'
        cb.set_label(r'$v_x~\mathrm{[km~s^{-1}]}$',size=14)#, size='large')
        cb.ax.tick_params(labelsize=14)

        #cb.set_label(r'$||\vec{v}||~[km/s]$', size='large')

        ratiov = map.shape[0] / map.shape[1]

        cbs = fig.colorbar(stream.lines, ax=ax, orientation='vertical',fraction=0.047*ratiov)#, pad=0.14, shrink=0.95)
        cbs.set_label(r'$||\vec{v}||~\mathrm{[km~s^{-1}]}$',size=14)#, size='large')
        cbs.ax.tick_params(labelsize=14)

        #cb.ax.tick_params(labelsize=14)
        #cb.ax.set_title('$T~\mathrm{[K]}$', fontsize=14)

        #ax.get_xaxis().set_visible(False)
        #ax.get_yaxis().set_visible(False)

        #ax.set_xlabel("x [Mpc]")
        #ax.set_ylabel("y [Mpc]")



        #cb.set_label(r'$||\vec{v}||~[km/s]$', size='large')

        #cb.set_label(r'$|v_{x0.08}-v_{x0.04}|~[km/s]$', size='large')

    elif (type == "f"):

        #vx_virgo = -509.1301
        #vy_virgo = 228.9488
        #vz_virgo = -131.9249
        fig, ax = plt.subplots()

        nbin = 10

        #x, y = np.meshgrid(np.linspace(-dimx, dimx, int(nx/nbin)), np.linspace(-dimy, dimy, int(ny/nbin)))

        #np.save("mesh_x_left_fil",x)
        #np.save("mesh_y_left_fil",y)

        x = np.load("mesh_x.npy")
        y = np.load("mesh_y.npy")

        #print(np.shape(x))
        #print(np.shape(y))

        #sys.exit()

        #x = np.load("mesh_x_left_fil.npy")
        #y = np.load("mesh_y_left_fil.npy")

        print("mean calculation")

        #map -= vz_virgo
        #map2 -= vy_virgo
        #map3 -= vz_virgo

        #map_mean = np.sqrt(map**2 + map2**2)# + map3**2)

        #np.save("mean_v_2D_map.npy",map_mean)

        #map_mean = np.load("mean_v_map.npy")

        #print("max=",np.max(map_mean))
        #print("min=",np.min(map_mean))

        #print(np.shape(map))
        #print(np.shape(map2))

        #mean10_map = np.array([[np.mean(map[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for j in range(int(ny / nbin))] for i in range(int(nx / nbin))])

        #mean10_map = np.array([[np.nanmean(map[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for j in range(int(ny / nbin))] for i in range(int(nz / nbin))])

        print("mean10_map done")

        #np.save("stream_fz_left_fil_minus_8Mpc.npy",mean10_map)

        #mean10_map2 = np.array([[np.mean(map2[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for j in range(int(ny / nbin))] for i in range(int(nx / nbin))])

        #mean10_map2 = np.array([[np.nanmean(map2[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for j in range(int(ny / nbin))] for i in range(int(nz / nbin))])

        print("mean10_map2 done")

        #np.save("stream_fy_left_fil_minus_8Mpc.npy",mean10_map2)

        #mean10_map = np.load("stream_vx.npy")
        #mean10_map2 = np.load("stream_vy.npy")

        #mean10_map = np.load("stream_vx_left_fil.npy")
        #mean10_map2 = np.load("stream_vy_left_fil.npy")

        #mean10_map = np.load("stream_vz_left_fil_minus_2Mpc.npy")
        #mean10_map2 = np.load("stream_vy_left_fil_minus_2Mpc.npy")

        mean10_map = np.load("stream_fz_left_fil_minus_8Mpc.npy")
        mean10_map2 = np.load("stream_fy_left_fil_minus_8Mpc.npy")


        #print("max mean map vx =",np.max(mean10_map))
        #print("min mean map vx =",np.min(mean10_map))
        #print("max mean map vy =",np.max(mean10_map2))
        #print("min mean map vy =",np.min(mean10_map2))

        #sys.exit()

        #mean10_map3 = np.array([[np.mean(map3[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for i in range(int(nx / nbin))] for j in range(int(ny / nbin))])

        #d = np.array([[np.log10(np.sqrt((cen_x - (i + 0.5)) ** 2 + (cen_y - (j + 0.5)) ** 2) * 5.6262) for i in range(max_x)] for j in range(max_y)])

        #print("mean calculated")

        print("mean mesh calculation")

        #print(np.shape(mean10_map))
        #print(np.shape(mean10_map2))

        print("len nan", len(mean10_map[np.isnan(mean10_map) == True]), 'ratio, ', len(mean10_map[np.isnan(mean10_map) == True]) / len(mean10_map))
        print("len nan", len(mean10_map2[np.isnan(mean10_map2) == True]), 'ratio, ', len(mean10_map2[np.isnan(mean10_map2) == True]) / len(mean10_map2))

        map_mean_mesh = np.sqrt(mean10_map**2 + mean10_map2**2)

        print("len nan", len(map_mean_mesh[np.isnan(map_mean_mesh) == True]), 'ratio, ', len(map_mean_mesh[np.isnan(map_mean_mesh) == True]) / len(map_mean_mesh))

        map_mean_mesh[np.isnan(map_mean_mesh) == True] = 0

        #map_mean_mesh = np.sqrt(mean10_map[np.isnan(mean10_map) == False] ** 2 + mean10_map2[np.isnan(mean10_map2) == False] ** 2)

        #print(np.shape(map_mean_mesh))

        #map_diff = np.abs(map-map2)

        print("max mean map =",np.max(map_mean_mesh))
        print("min mean map =",np.min(map_mean_mesh))

        print("min map",np.min(map))
        print("max map",np.max(map))

        #norm = Normalize(vmin = np.min(map_mean_mesh), vmax = np.max(map_mean_mesh))

        #map[map>0]=np.log10(map[map>0])
        #map[map<=0]=-np.log10(-map[map<=0])

        #plt.hist(np.log10(map_mean_mesh[map_mean_mesh>0].flatten()), bins=100, histtype='step', color='black', linewidth=2,log=True)
        #plt.show()

        #sys.exit()

        #lw = 2 * map_mean_mesh / np.max(map_mean_mesh)

        #print("diff",map-map2)

        #print("lw",np.shape(lw))
        #sys.exit()
        #plt.quiver(x, y, mean10_map2, mean10_map3)
        #stream = ax.streamplot(x, y, mean10_map, mean10_map2, linewidth=0.5, color=map_mean_mesh, cmap='Ice', density=[4,4],norm=colors.LogNorm(vmin=1e-24, vmax=1e-19))#,norm=norm) # color="black" #,vmin=-1600,vmax=1600) #color=mean10_map2, cmap="BR"
        #plt.show()

        #map -= np.mean(map)
        im = ax.imshow(map, cmap="DryWet", origin="lower", alpha=1, extent=dim,norm=colors.SymLogNorm(linthresh=1e-26, linscale=0.3,vmin=-1e-18, vmax=1e-18, base=10))#, vmin=-1600, vmax=1600) #BR #origin='lower'

        #im = ax.imshow(map_diff, cmap="Greys", origin="lower", alpha=1, extent=dim)#, vmin=-1600, vmax=1600)  # BR #origin='lower'
        #im = ax.imshow(map_mean, cmap="rocket", origin="lower", alpha=1, extent=dim)#, vmin=-1600, vmax=1600)  # BR #origin='lower'
        #show_radii()
        cb = fig.colorbar(im, ax=ax, pad=0.04,orientation='vertical')
        cb.set_label(r'$\rho v_z~[kg~m^{-2}~s^{-1}]$', size='large')

        #cbs = fig.colorbar(stream.lines, ax=ax, orientation='vertical', pad=0.14)
        #cbs.set_label(r'$||\rho \vec{v}||~[kg~m^{-2}s^{-1}]$', size='large')


        ax.set_xlabel("x [Mpc]")
        ax.set_ylabel("y [Mpc]")

        #cb.set_label(r'$||\vec{v}||~[km/s]$', size='large')

        #cb.set_label(r'$|v_{x0.08}-v_{x0.04}|~[km/s]$', size='large')

    elif (type == "vproj"):
        fig, ax = plt.subplots(figsize=(10, 8), facecolor='white')

        vx_virgo = -509.1301
        vy_virgo = 228.9488
        vz_virgo = -131.9249

        if proj=="x":
            #map -= vx_virgo
            mid = int(ny/2)
            #min = int(mid - ny/4)
            #max = int(mid + ny/4)

        elif proj=="y":
            #map -= vy_virgo
            mid = int(nx/2)
            #min = int(mid - nx/4)
            #max = int(mid + nx/4)
        elif proj=="z":
            #map -= vz_virgo
            mid = int(nx/2)

            #min = int(mid - nx/4)
            #max = int(mid + nx/4)

        px_size = 737.441 / 2 ** lvl

        npx500 = int(0.5 / px_size)

        npx1000 = int(1.0 / px_size)

        npx1200 = int(1.2 / px_size)

        npx2250 = int(2.25 / px_size)

        print("npx2250",npx2250)

        npx=npx500

        #sys.exit()

        print("px_size",px_size,"npx",npx)


        #print("npx500",npx500)

        #sys.exit()

        min = mid - npx
        max = mid + npx

        #print("min",min,"max",max)

        v_list = map.flatten()

        #np.save("./v_list_pdf/vy_list_from_1Mpc2_map.npy", v_list)
        #print('file saved')

        #sys.exit()

        def v_hist():

            vz_list_100kpc = np.load("vz_list_from_100kpc2_map.npy")
            vz_Trange_list_100kpc = np.load("vz_Trange_list_from_100kpc2_map.npy")

            vx_list_100kpc = np.load("vx_list_from_100kpc2_map.npy")
            vx_Trange_list_100kpc = np.load("vx_Trange_list_from_100kpc2_map.npy")

            vy_list_100kpc = np.load("vy_list_from_100kpc2_map.npy")
            vy_Trange_list_100kpc = np.load("vy_Trange_list_from_100kpc2_map.npy")





            plt.hist(map.flatten(), bins=100, alpha=0.6, color='blue', label='100kpc**2, mw velocity')

            #plt.hist(vy_list_100kpc, bins=100, alpha=0.6, color='blue', label='vy, 100kpc**2')
            plt.hist(vz_list_100kpc, bins=100, alpha=0.6, color='green', label='vz, 100kpc**2, ew velocity')
            #plt.hist(vx_list_100kpc, bins=100, alpha=0.6, color='red', label='vx, 100kpc**2')

            #plt.hist(vy_Trange_list_100kpc, bins=100, alpha=0.6, color='blue', label='100kpc**2, XRISM T range')
            plt.hist(vz_Trange_list_100kpc, bins=100, alpha=0.6, color='orange', label='100kpc**2,ew velovity, XRISM T range')
            #plt.hist(vx_Trange_list_100kpc, bins=100, alpha=0.6, color='red', label='100kpc**2, XRISM T range')



            #plt.hist(v_Trange_list_100kpc, bins=100, alpha=0.6, color='green', label='100kpc**2, XRISM T range')
            plt.grid(b=None)
            plt.legend()
            plt.ylabel("Number of cells")
            plt.xlabel("V [km/s]")
            plt.yscale('log')
            #plt.title("Vz distribution on ew sightline velocity along z axis")
            plt.title("Vz distribution along z axis")
            #plt.title("Velocity distribution along sightlines, 100kpc**2 maps, XRISM T range")



            plt.show()

            sys.exit()

        #v_hist()


        map_small = map[min:max,min:max]
        #dim_small = [-0.5, 0.5, -0.5, 0.5]
        dim_small = [-npx*px_size,npx*px_size,-npx*px_size,npx*px_size]

        size=22
        #nbin = 10

        #x, y = np.meshgrid(np.linspace(-dimx, dimx, int(nx/nbin)), np.linspace(-dimy, dimy, int(ny/nbin)))

        #print("mean calculation")

        #mean10_map2 = np.array([[np.mean(map2[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for i in range(int(nx / nbin))] for j in range(int(ny / nbin))])

        #print("mean10_map2 done")

        #mean10_map3 = np.array([[np.mean(map3[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for i in range(int(nx / nbin))] for j in range(int(ny / nbin))])

        #d = np.array([[np.log10(np.sqrt((cen_x - (i + 0.5)) ** 2 + (cen_y - (j + 0.5)) ** 2) * 5.6262) for i in range(max_x)] for j in range(max_y)])

        #print("mean calculated")

        #plt.quiver(x, y, mean10_map2, mean10_map3)
        #plt.streamplot(x, y, mean10_map2, mean10_map3, color='white', linewidth=0.1, density=3)
        #plt.show()

        #map -= np.mean(map)
        #im = ax.imshow(map_small, cmap="BR", origin='lower', alpha=1, extent=dim_small, vmin=-600, vmax=600)
        im = ax.imshow(map, cmap="BR", origin='lower', alpha=1, extent=dim, vmin=-600, vmax=600) #vmin=-2000, vmax=2000
        show_radii()
        #show_quadrants()
        cb = fig.colorbar(im, ax=ax)

        cb.ax.tick_params(labelsize=size)
        ax.set_xlabel("Mpc",size=size)
        ax.set_ylabel("Mpc",size=size)
        plt.yticks(fontsize=size)
        plt.xticks(fontsize=size)

        if proj=="x":
            cb.set_label('$v_{x}~[km~s^{-1}]$', size=size)
        elif proj=="y":
            cb.set_label('$v_{y}~[km~s^{-1}]$', size=size)
        elif proj=="z":
            cb.set_label('$v_{z}~[km~s^{-1}]$', size=size)

        #print('test')

        #cb.set_label('$v_{z}~[km/s]$', size=size)

        #show_radii()

    elif (type == "d"):
        fig, ax = plt.subplots(figsize=(10, 8), facecolor='white')

        vx_virgo = -509.1301
        vy_virgo = 228.9488
        vz_virgo = -131.9249

        if proj=="x":
            #map -= vx_virgo
            mid = int(ny/2)
            #min = int(mid - ny/4)
            #max = int(mid + ny/4)

        elif proj=="y":
            #map -= vy_virgo
            mid = int(nx/2)
            #min = int(mid - nx/4)
            #max = int(mid + nx/4)
        elif proj=="z":
            #map -= vz_virgo
            mid = int(nx/2)

            #min = int(mid - nx/4)
            #max = int(mid + nx/4)

        px_size = 737.441 / 2 ** lvl

        npx500 = int(0.5 / px_size)

        npx1000 = int(1.0 / px_size)

        npx1200 = int(1.2 / px_size)

        npx=npx1200

        print("px_size",px_size,"npx500",npx500)


        #print("npx500",npx500)

        #sys.exit()

        min = mid - npx
        max = mid + npx

        print("min",min,"max",max)





        map_small = map[min:max,min:max]
        #dim_small = [-0.5, 0.5, -0.5, 0.5]
        dim_small = [-npx*px_size,npx*px_size,-npx*px_size,npx*px_size]

        size=22
        #nbin = 10

        #x, y = np.meshgrid(np.linspace(-dimx, dimx, int(nx/nbin)), np.linspace(-dimy, dimy, int(ny/nbin)))

        #print("mean calculation")

        #mean10_map2 = np.array([[np.mean(map2[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for i in range(int(nx / nbin))] for j in range(int(ny / nbin))])

        #print("mean10_map2 done")

        #mean10_map3 = np.array([[np.mean(map3[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for i in range(int(nx / nbin))] for j in range(int(ny / nbin))])

        #d = np.array([[np.log10(np.sqrt((cen_x - (i + 0.5)) ** 2 + (cen_y - (j + 0.5)) ** 2) * 5.6262) for i in range(max_x)] for j in range(max_y)])

        #print("mean calculated")

        #plt.quiver(x, y, mean10_map2, mean10_map3)
        #plt.streamplot(x, y, mean10_map2, mean10_map3, color='white', linewidth=0.1, density=3)
        #plt.show()

        #map -= np.mean(map)
        #im = ax.imshow(map_small, cmap="BR", origin='lower', alpha=1, extent=dim_small, vmin=-400, vmax=400)
        im = ax.imshow(map, cmap="Sunset_r", origin='lower', alpha=1, extent=dim, vmin=0, vmax=800)
        show_radii()
        #show_quadrants()
        cb = fig.colorbar(im, ax=ax,)

        cb.ax.tick_params(labelsize=size)
        ax.set_xlabel("x [Mpc]",size=size)
        ax.set_ylabel("y [Mpc]",size=size)
        plt.yticks(fontsize=size)
        plt.xticks(fontsize=size)

        if proj=="x":
            cb.set_label('$w_{x}~[km~s^{-1}]$', size=size)
        elif proj=="y":
            cb.set_label('$w_{y}~[km~s^{-1}]$', size=size)
        elif proj=="z":
            cb.set_label('$w_{z}~[km~s^{-1}]$', size=size)

        #cb.set_label('$v_{z}~[km/s]$', size=size)

        #show_radii()

    elif (type == "EM"):

        map[map == 0] = 1e-15
        #map = np.log10(map)
        map[np.isnan(map) == True] = 0

        fig, ax = plt.subplots()
        #im = ax.imshow(map, cmap="afmhot", origin='lower', alpha=1, extent=dim, vmin=-14, vmax=-1)
        im = ax.imshow(map, cmap="afmhot", origin='lower', alpha=1, extent=dim)#, vmin=, vmax=-1)# 3D: 0 6.5
        #show_radii()
        ax.axis('off')

        # divider = make_axes_locatable(ax)
        # cax = divider.new_vertical(size='5%', pad=0.6, pack_start=True)
        # fig.add_axes(cax)
        # cb = fig.colorbar(im, cax=cax, orientation='horizontal')

        # cb.set_label('$log_{10}(n_e[cm^{-2}])$',size='large')

        cb = fig.colorbar(im,orientation='horizontal')
        cb.ax.tick_params(labelsize=20)
        cb.set_label('$log_{10}(EM~[\mathrm{cm^{-6}~Mpc}])$', size=30)

    elif (type== "SD"):

        map[map == 0] = 1
        map = np.log10(map)
        map[np.isnan(map) == True] = 1.5

        fig, ax = plt.subplots()
        im = ax.imshow(map, cmap="gray", origin='lower', alpha=1, extent=dim)#, vmin=3, vmax=9)  # 3D: 0 6.5
        show_radii()
        ax.axis('off')

        # divider = make_axes_locatable(ax)
        # cax = divider.new_vertical(size='5%', pad=0.6, pack_start=True)
        # fig.add_axes(cax)
        # cb = fig.colorbar(im, cax=cax, orientation='horizontal')

        # cb.set_label('$log_{10}(n_e[cm^{-2}])$',size='large')

        cb = fig.colorbar(im,orientation='horizontal')
        #cb.set_label('$log_{10}(EM[cm^{-3}])$', size='large')
        cb.set_label('$log_{10}(SD~[\mathrm{M_{\odot}~kpc^{-2}}])$', size=30)
        cb.ax.tick_params(labelsize=20)

    elif (type == "K"):
        fig, ax = plt.subplots()
        ##test in K :
        #map /= (kb / 1.602e-16)

        map -= np.mean(map)

        #map=np.log10(map)


        im=ax.imshow(map, cmap="coolwarm", origin='lower', alpha=1, extent=dim, vmin=-5, vmax=5)#, vmin=0, vmax=6)  "hot"
        #show_radii()
        #ax.axis('off')


        cb = fig.colorbar(im, ax=ax)
        #divider = make_axes_locatable(ax)
        #cax = divider.new_vertical(size='5%', pad=0.6, pack_start=True)
        #fig.add_axes(cax)
        #cb = fig.colorbar(im, cax=cax, orientation='horizontal')

        #cb.set_label('$T[keV]$',size='large')
        #cb.set_label('$log(T[K])$', size='large')
        cb.set_label('$K\,[keV\,cm^2]$', size='large')

    elif (type=='T-v'):
        vx_virgo = -509.1301
        vy_virgo = 228.9488
        vz_virgo = -131.9249
        fig, ax = plt.subplots()

        nbin = 10

        # x, y = np.meshgrid(np.linspace(-dimx, dimx, int(nx/nbin)), np.linspace(-dimy, dimy, int(ny/nbin)))

        # np.save("mesh_x_left_fil",x)
        # np.save("mesh_y_left_fil",y)

        x = np.load("mesh_x.npy")
        y = np.load("mesh_y.npy")

        # print(np.shape(x))
        # print(np.shape(y))

        # sys.exit()

        # x = np.load("mesh_x_left_fil.npy")
        # y = np.load("mesh_y_left_fil.npy")

        #print("mean calculation")

        #map -= vz_virgo
        #map2 -= vy_virgo
        # map3 -= vz_virgo

        # map_mean = np.sqrt(map**2 + map2**2)# + map3**2)

        # np.save("mean_v_2D_map.npy",map_mean)

        # map_mean = np.load("mean_v_map.npy")

        # print("max=",np.max(map_mean))
        # print("min=",np.min(map_mean))

        # print(np.shape(map))
        # print(np.shape(map2))

        # mean10_map = np.array([[np.mean(map[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for j in range(int(ny / nbin))] for i in range(int(nx / nbin))])

        # mean10_map = np.array([[np.nanmean(map[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for j in range(int(ny / nbin))] for i in range(int(nz / nbin))])

        #print("mean10_map done")

        # np.save("stream_vz_left_fil_minus_4Mpc.npy",mean10_map)

        # mean10_map2 = np.array([[np.mean(map2[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for j in range(int(ny / nbin))] for i in range(int(nx / nbin))])

        # mean10_map2 = np.array([[np.nanmean(map2[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for j in range(int(ny / nbin))] for i in range(int(nz / nbin))])

        #print("mean10_map2 done")

        # np.save("stream_vy_left_fil_minus_4Mpc.npy",mean10_map2)

        # mean10_map = np.load("stream_vx.npy")
        # mean10_map2 = np.load("stream_vy.npy")

        # mean10_map = np.load("stream_vx_left_fil.npy")
        # mean10_map2 = np.load("stream_vy_left_fil.npy")

        mean10_map = np.load("stream_vz_left_fil_minus_8Mpc.npy")
        mean10_map2 = np.load("stream_vy_left_fil_minus_8Mpc.npy")

        # print("max mean map vx =",np.max(mean10_map))
        # print("min mean map vx =",np.min(mean10_map))
        # print("max mean map vy =",np.max(mean10_map2))
        # print("min mean map vy =",np.min(mean10_map2))

        # sys.exit()

        # mean10_map3 = np.array([[np.mean(map3[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for i in range(int(nx / nbin))] for j in range(int(ny / nbin))])

        # d = np.array([[np.log10(np.sqrt((cen_x - (i + 0.5)) ** 2 + (cen_y - (j + 0.5)) ** 2) * 5.6262) for i in range(max_x)] for j in range(max_y)])

        # print("mean calculated")

        print("mean mesh calculation")

        # print(np.shape(mean10_map))
        # print(np.shape(mean10_map2))

        #print("len nan", len(mean10_map[np.isnan(mean10_map) == True]), 'ratio, ',
        #      len(mean10_map[np.isnan(mean10_map) == True]) / len(mean10_map))
        #print("len nan", len(mean10_map2[np.isnan(mean10_map2) == True]), 'ratio, ',
        #      len(mean10_map2[np.isnan(mean10_map2) == True]) / len(mean10_map2))

        map_mean_mesh = np.sqrt(mean10_map ** 2 + mean10_map2 ** 2)

        #print("len nan", len(map_mean_mesh[np.isnan(map_mean_mesh) == True]), 'ratio, ',
        #      len(map_mean_mesh[np.isnan(map_mean_mesh) == True]) / len(map_mean_mesh))

        map_mean_mesh[np.isnan(map_mean_mesh) == True] = 0

        # map_mean_mesh = np.sqrt(mean10_map[np.isnan(mean10_map) == False] ** 2 + mean10_map2[np.isnan(mean10_map2) == False] ** 2)

        # print(np.shape(map_mean_mesh))

        # map_diff = np.abs(map-map2)

        # print("max mean map =",np.max(map_mean_mesh))
        # print("min mean map =",np.min(map_mean_mesh))

        # norm = Normalize(vmin = np.min(map_mean_mesh), vmax = np.max(map_mean_mesh))

        # print("max mean map =", np.max(mean10_map2))
        # print("min mean map =", np.min(mean10_map2))

        norm = Normalize(vmin=0, vmax=700)

        map /= (kb / 1.602e-16)
        #map = np.log10(map)


        # plt.hist(map_mean_mesh.flatten(), bins=100, range=(0, 2000), histtype='step', color='black', linewidth=2)
        # plt.show()

        # sys.exit()

        # lw = 2 * map_mean_mesh / np.max(map_mean_mesh)

        # print("diff",map-map2)

        # print("lw",np.shape(lw))
        # sys.exit()
        # plt.quiver(x, y, mean10_map2, mean10_map3)
        stream = ax.streamplot(x, y, mean10_map, mean10_map2, linewidth=0.5,
                               color=np.sqrt(mean10_map ** 2 + mean10_map2 ** 2), cmap='Greys_r', density=[3, 3],
                               norm=norm)  # ,alpha=0.8) # color="black" #,vmin=-1600,vmax=1600) #color=mean10_map2, cmap="BR"
        # plt.show()

        # map -= np.mean(map)
        #im = ax.imshow(map, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600, vmax=1600)  # BR #origin='lower'

        im = ax.imshow(map, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))

        # im = ax.imshow(map_diff, cmap="Greys", origin="lower", alpha=1, extent=dim)#, vmin=-1600, vmax=1600)  # BR #origin='lower'
        # im = ax.imshow(map_mean, cmap="rocket", origin="lower", alpha=1, extent=dim)#, vmin=-1600, vmax=1600)  # BR #origin='lower'
        # show_radii()
        cb = fig.colorbar(im, ax=ax, pad=0.04, orientation='vertical')
        cb.set_label('$T[K]$', size='large')

        cbs = fig.colorbar(stream.lines, ax=ax, orientation='vertical', pad=0.14)
        cbs.set_label(r'$||\vec{v}||~[km/s]$', size='large')


        ax.set_xlabel("x [Mpc]")
        ax.set_ylabel("y [Mpc]")

        # cb.set_label(r'$||\vec{v}||~[km/s]$', size='large')

        # cb.set_label(r'$|v_{x0.08}-v_{x0.04}|~[km/s]$', size='large')

    elif (type=='ne-rhov'):

        fig, ax = plt.subplots()

        nbin = 10

        # x, y = np.meshgrid(np.linspace(-dimx, dimx, int(nx/nbin)), np.linspace(-dimy, dimy, int(ny/nbin)))

        # np.save("mesh_x_left_fil",x)
        # np.save("mesh_y_left_fil",y)

        x = np.load("mesh_x.npy")
        y = np.load("mesh_y.npy")

        # print(np.shape(x))
        # print(np.shape(y))

        # sys.exit()

        # x = np.load("mesh_x_left_fil.npy")
        # y = np.load("mesh_y_left_fil.npy")

        #print("mean calculation")

        #map -= vz_virgo
        #map2 -= vy_virgo
        # map3 -= vz_virgo

        # map_mean = np.sqrt(map**2 + map2**2)# + map3**2)

        # np.save("mean_v_2D_map.npy",map_mean)

        # map_mean = np.load("mean_v_map.npy")

        # print("max=",np.max(map_mean))
        # print("min=",np.min(map_mean))

        # print(np.shape(map))
        # print(np.shape(map2))

        # mean10_map = np.array([[np.mean(map[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for j in range(int(ny / nbin))] for i in range(int(nx / nbin))])

        # mean10_map = np.array([[np.nanmean(map[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for j in range(int(ny / nbin))] for i in range(int(nz / nbin))])

        #print("mean10_map done")

        # np.save("stream_vz_left_fil_minus_4Mpc.npy",mean10_map)

        # mean10_map2 = np.array([[np.mean(map2[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for j in range(int(ny / nbin))] for i in range(int(nx / nbin))])

        # mean10_map2 = np.array([[np.nanmean(map2[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for j in range(int(ny / nbin))] for i in range(int(nz / nbin))])

        #print("mean10_map2 done")

        # np.save("stream_vy_left_fil_minus_4Mpc.npy",mean10_map2)

        # mean10_map = np.load("stream_vx.npy")
        # mean10_map2 = np.load("stream_vy.npy")

        # mean10_map = np.load("stream_vx_left_fil.npy")
        # mean10_map2 = np.load("stream_vy_left_fil.npy")

        mean10_map = np.load("stream_fz_left_fil_minus_8Mpc.npy")
        mean10_map2 = np.load("stream_fy_left_fil_minus_8Mpc.npy")

        # print("max mean map vx =",np.max(mean10_map))
        # print("min mean map vx =",np.min(mean10_map))
        # print("max mean map vy =",np.max(mean10_map2))
        # print("min mean map vy =",np.min(mean10_map2))

        # sys.exit()

        # mean10_map3 = np.array([[np.mean(map3[i * nbin:(i + 1) * nbin, j * nbin:(j + 1) * nbin]) for i in range(int(nx / nbin))] for j in range(int(ny / nbin))])

        # d = np.array([[np.log10(np.sqrt((cen_x - (i + 0.5)) ** 2 + (cen_y - (j + 0.5)) ** 2) * 5.6262) for i in range(max_x)] for j in range(max_y)])

        # print("mean calculated")

        print("mean mesh calculation")

        # print(np.shape(mean10_map))
        # print(np.shape(mean10_map2))

        #print("len nan", len(mean10_map[np.isnan(mean10_map) == True]), 'ratio, ',
        #      len(mean10_map[np.isnan(mean10_map) == True]) / len(mean10_map))
        #print("len nan", len(mean10_map2[np.isnan(mean10_map2) == True]), 'ratio, ',
        #      len(mean10_map2[np.isnan(mean10_map2) == True]) / len(mean10_map2))

        map_mean_mesh = np.sqrt(mean10_map ** 2 + mean10_map2 ** 2)

        #print("len nan", len(map_mean_mesh[np.isnan(map_mean_mesh) == True]), 'ratio, ',
        #      len(map_mean_mesh[np.isnan(map_mean_mesh) == True]) / len(map_mean_mesh))

        map_mean_mesh[np.isnan(map_mean_mesh) == True] = 0

        # map_mean_mesh = np.sqrt(mean10_map[np.isnan(mean10_map) == False] ** 2 + mean10_map2[np.isnan(mean10_map2) == False] ** 2)

        # print(np.shape(map_mean_mesh))

        # map_diff = np.abs(map-map2)

        # print("max mean map =",np.max(map_mean_mesh))
        # print("min mean map =",np.min(map_mean_mesh))

        # norm = Normalize(vmin = np.min(map_mean_mesh), vmax = np.max(map_mean_mesh))

        # print("max mean map =", np.max(mean10_map2))
        # print("min mean map =", np.min(mean10_map2))

        #map /= np.mean(map)
        #map = np.log10(map)



        # plt.hist(map_mean_mesh.flatten(), bins=100, range=(0, 2000), histtype='step', color='black', linewidth=2)
        # plt.show()

        # sys.exit()

        # lw = 2 * map_mean_mesh / np.max(map_mean_mesh)

        # print("diff",map-map2)

        # print("lw",np.shape(lw))
        # sys.exit()
        # plt.quiver(x, y, mean10_map2, mean10_map3)
        stream = ax.streamplot(x, y, mean10_map, mean10_map2, linewidth=0.5, color=np.sqrt(mean10_map ** 2 + mean10_map2 ** 2), cmap='Greys', density=[3, 3], norm = colors.LogNorm(vmin=1e-24, vmax=1e-19))  # ,alpha=0.8) # color="black" #,vmin=-1600,vmax=1600) #color=mean10_map2, cmap="BR"

        im = ax.imshow(map, cmap="cividis", origin='lower', alpha=1, extent=dim,norm=colors.LogNorm(vmin=1e13,vmax=1e17))#,vmin=-2,vmax=3)  # ,vmin=14,vmax=19)#,vmin=16)#, vmin=16, vmax=21)  # 3D: 0 6.5 #filament : vmin=14, vmax=19

        # im = ax.imshow(map_diff, cmap="Greys", origin="lower", alpha=1, extent=dim)#, vmin=-1600, vmax=1600)  # BR #origin='lower'
        # im = ax.imshow(map_mean, cmap="rocket", origin="lower", alpha=1, extent=dim)#, vmin=-1600, vmax=1600)  # BR #origin='lower'
        # show_radii()
        cb = fig.colorbar(im, ax=ax, pad=0.04, orientation='vertical')
        #cb.set_label('$log(T[K])$', size='large')

        cb.set_label('$log_{10}(n_e[cm^{-2}])$', size='large')
        #cb.set_label('$log_{10}(n_e/\overline{n_e})$', size='large')

        cbs = fig.colorbar(stream.lines, ax=ax, orientation='vertical', pad=0.14)
        cbs.set_label(r'$\rho v_z~[kg~m^{-2}~s^{-1}]$', size='large')


        ax.set_xlabel("x [Mpc]")
        ax.set_ylabel("y [Mpc]")

        # cb.set_label(r'$||\vec{v}||~[km/s]$', size='large')

        # cb.set_label(r'$|v_{x0.08}-v_{x0.04}|~[km/s]$', size='large')

    elif (type=="mach"):

        load_T=1
        if load_T==1:

            file='./maps/high_res/velocity/15f15_analysis/map_high_15f15_x_map_T_sl_Tsup7_5Mpc2.bin'
            proj='x'

            h = FortranFile(file, 'r')
            nx, ny, nz = h.read_ints()
            cen_x, cen_y, cen_z = h.read_reals()
            print("nx,ny,nz", nx, ny, nz)
            if proj == "x":
                ncell = nz * ny
            elif proj == "y":
                ncell = nx * nz
            elif proj == "z":
                ncell = nx * ny
            print("ncell", ncell)

            map2 = np.zeros(ncell)

            map2 = ftp.f90_to_py.read_map_file(ncell, file, 0)

            if proj == "x":
                map2 = np.reshape(map2, (nz, ny))
                # map2 = np.reshape(map2, (nz, ny))
                # map3 = np.reshape(map3, (nz, ny))
                nl = nx
            elif proj == "y":
                map2 = np.reshape(map2, (nx, nz))
                # map2 = np.reshape(map2, (nx, nz))
                # map3 = np.reshape(map3, (nx, nz))
                nl = ny
            elif proj == "z":
                map2 = np.reshape(map2, (ny, nx))
                # map2 = np.reshape(map2, (ny, nx))
                # map3 = np.reshape(map3, (ny, nx))
                nl = nz

        fig, ax = plt.subplots(figsize=(8, 8), facecolor='white')

        #print("shape map2",np.shape(map2))

        map=np.abs(map)

        print("T", map2)
        print("mean T", np.mean(map2))
        map2*=(kev/kb)
        print("T", map2)
        print("mean T", np.mean(map2))
        # convert to K
        #c_s = 1170
        gamma=5/3
        c_s = np.sqrt((gamma*kb*map2)/(mu*mp))

        print("cs",c_s)
        print("mean cs",np.mean(c_s))

        #sys.exit()

        map /= (c_s/1e3)  # convert to km/s

        #print("max map",np.max(map), "min map",np.min(map))
        #sys.exit()
        ##test in K :
        #map /= (kb / 1.602e-16)
        # map=np.log10(map)

        #map = np.log10(map)

        map_above_1 = np.ma.masked_less(map, 1)  # Mask values below 1
        map_below_1 = np.ma.masked_greater(map, 1)

        ratio = map.shape[1] / map.shape[0]

        im2 = ax.imshow(map, cmap="Blues_r", origin='lower', alpha=1, extent=dim, vmin=0.01, vmax=1)  # , vmin=0, vmax=6)  "hot" "flare"

        cb2 = fig.colorbar(im2, ax=ax, orientation="horizontal",pad=0.15,fraction=0.047 * ratio)
        #cb2.set_label(r'$M=\frac{\sigma}{c_s}$', size=14)
        cb2.set_label(r'$M=\frac{v}{c_s}$', size=14)
        #im1 = ax.imshow(map_above_1, cmap="Blues", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1, vmax=100))  # , vmin=0, vmax=6)  "hot" "flare"
        # show_radii()
        # ax.axis('off')
        #im2 = ax.imshow(map_below_1, cmap="Mono_r", origin='lower', alpha=1, extent=dim, vmin=0.1, vmax=1)  # , vmin=0, vmax=6)  "hot" "flare"
        plt.sca(ax)
        show_radii()
        size=22
        cb2.ax.tick_params(labelsize=size)
        ax.set_xlabel("Mpc", size=size)
        ax.set_ylabel("Mpc", size=size)
        plt.yticks(fontsize=size)
        plt.xticks(fontsize=size)


        #ratio = map.shape[1] / map.shape[0]

        #
        #cb2 = fig.colorbar(im2, ax=ax, orientation="horizontal",fraction=0.047 * ratio,pad=0.2)  # , label='$T~\mathrm{[K]}$',size="large")#,pad=0.1)
        #cb2.ax.tick_params(labelsize=14)
        #cb2.set_label(r'$M=\frac{v}{c_s}$ (0.1 to 1, linear scale)', size=14)

        #cb1 = fig.colorbar(im1, ax=ax, orientation="horizontal",fraction=0.047 * ratio,pad=0.15)  # , label='$T~\mathrm{[K]}$',size="large")#,pad=0.1)
        #cb1.ax.tick_params(labelsize=14)
        #cb1.set_label(r'$M=\frac{v}{c_s}$ (1 to 100, log scale)', size=14)



        #cb.ax.xaxis.set_major_locator(LogLocator(base=10))  # Base 10 log scale
        #cb.ax.xaxis.set_major_formatter(LogFormatter())
        #cbar_ticks = LogLocator(base=10)  # Base 10 log scale
        #cb.set_ticks(cbar_ticks)
        #cb.set_ticklabels([0.1, 0.3, 1, 3, 10])  # Set specific tick labels for clarity

    elif (type=="lvl"):
        fig, ax = plt.subplots(figsize=(8, 8), facecolor='white')

        map = 737441 / (2 ** map)  # convert to kpc

        list=np.array([737441/2**i for i in range(12,22)])
        #list=np.round(list,decimals=1)
        #print("list",list)
        #sys.exit()
        im = ax.imshow(map, cmap="Crest", origin='lower', alpha=1, extent=dim,norm=colors.LogNorm(vmin=0.3, vmax=200))  # ,norm=colors.Normalize(vmin=1.7, vmax=12))  # , vmin=0, vmax=6)  "hot"
        #show_radii()
        # ax.axis('off')

        #ratio = map.shape[1] / map.shape[0]


        #

        size = 22
        cb = fig.colorbar(im,ticks=list)
        #cb.ax.set_yticks([])
        cb.ax.minorticks_off()
        cb.ax.set_yticks(list)
        cb.ax.set_yticklabels(["{:4.2f}".format(i) for i in list])#, ax=ax, orientation="horizontal", fraction=0.047 * ratio)  # , label='$T~\mathrm{[K]}$',size="large")#,pad=0.1)
        cb.set_label('Mean cell size [kpc]', size=size)
        cb.ax.tick_params(labelsize=20)
        ax.set_xlabel("Mpc", size=size)
        ax.set_ylabel("Mpc", size=size)
        plt.yticks(fontsize=size)
        plt.xticks(fontsize=size)
        #titre="Transverse slice: $\Delta x_{cen}$=-8 Mpc"
        titre="Longitudinal slice"



    #if (proj == "x"):
    #    plt.xlabel("y (Mpc)")
    #    plt.ylabel("z (Mpc)")
    #elif (proj == "y"):
    #    plt.xlabel("x (Mpc)")
    #    plt.ylabel("z (Mpc)")
    #elif (proj == "z"):
    #    plt.xlabel("x (Mpc)")
    #    plt.ylabel("y (Mpc)")
    #else:
    #    plt.xlabel("(Mpc)")
    #    plt.ylabel("(Mpc)")

    #print("cen_z",cen_z)
    #sys.exit()

    if proj == "z":
        ax.set_xlabel("x [Mpc]")
        ax.set_ylabel("y [Mpc]")

    elif proj == "y":
        ax.set_xlabel("x [Mpc]")
        ax.set_ylabel("z [Mpc]")

    elif proj == "x":
        ax.set_xlabel("z [Mpc]")
        ax.set_ylabel("y [Mpc]")

    #pos_z = cen_z/1e3 +2.3076327260135868
    pos_x = cen_x/1e3 +11.348717971989485

    plt.grid(b=None)

    depth = nl * (737.441 / 2 ** lvl)
    #titre = titre + ", depth = " + str(round(depth,3)) + " Mpc" #,

    #titre = titre + ", depth = 0.021 Mpc"  # ,
    #titre = "$\Delta x_{cen}$=" +str(round(pos_x,3)) + " Mpc"
    #titre = "Absolute velocity difference between 0.08 and 0.02 Mpc"

    #titre = "Compton-y"

    #titre=

    print("title : ",titre)

    size=22

    plt.title(titre, fontsize=size+6)
    #plt.xticks(fontsize=14)
    #plt.yticks(fontsize=14)
    #if proj=="z":
    #    plt.xlabel("Mpc", fontsize=14)
    #    plt.ylabel("y [Mpc]", fontsize=14)
    #elif proj=="y":
    #   plt.xlabel("x [Mpc]", fontsize=14)
    #    plt.ylabel("z [Mpc]", fontsize=14)
    #elif proj=="x":
    #    plt.xlabel("z [Mpc]", fontsize=14)
    #    plt.ylabel("y [Mpc]", fontsize=14)

    ax.set_xlabel("Mpc", size=size)
    ax.set_ylabel("Mpc", size=size)


    #titre = "./maps/high_res/filament/T_slices/x/savedfigs/T_dist_xvirgo="+str(round(pos_x,3))

    #plt.savefig(titre + ".png", dpi=300)

    #print('fig saved')

    if savecond==1:

        plt.savefig(savefile,format="png", dpi=600)
        print("fig saved")

    #sys.exit()

    plt.show()

    sys.exit()

#mapf90("T", 15, "x",'./maps/high_res/velocity/15f15_analysis/map_high_15f15_x_map_T_sl_Tsup7_5Mpc2.bin','',1,1,"./papers_plots/paper_4_vel/T_sl_x.png")
#mapf90("vproj", 15, "x",'./maps/high_res/velocity/15f15_analysis/map_high_15f15_x_map_vx_ew_Tsup7_5Mpc2.bin','',1,1,"./papers_plots/paper_4_vel/vx_ew.png")
#mapf90("ne", 15, "x",'./maps/high_res/velocity/15f15_analysis/map_high_15f15_x_map_ne_mw_Tsup5_gal8.5_10Mpc2.bin','',1,0,"./papers_plots/paper_4_vel/P_mw.png")

#mapf90("P", 15, "z",'./maps/high_res/velocity/15f15_analysis/map_high_15f15_cen_map_P_mw_Tsup_0.2kev_10Mpc2.bin','',1,0,"./papers_plots/paper_4_vel/P_mw.png")
#mapf90("P", 15, "z",'./maps/high_res/velocity/15f15_analysis/map_high_15f15_cen_map_P_mw_Tsup_0.2kev_15Mpc2.bin','',1,0,"./papers_plots/paper_4_vel/P_mw.png")
#mapf90("P", 15, "z",'./maps/high_res/velocity/15f15_analysis/map_high_15f15_cen_map_P_mw_Tsup_8e-3kev_15Mpc2_gc_1e8.5.bin','',1,0,"./papers_plots/paper_4_vel/P_mw.png")


#mapf90("d", 15, "z",'./maps/high_res/velocity/15f15_analysis/random_proj/map_high_15f15_rand_10_map_vd_ew_Tsup7_5Mpc2.bin',"",1,0,"")

#mapf90("lvl", 19, "z",'./maps/high_res/filament/map_high_19_xyz_full_map_lvl_z_d0.02.bin',"",1,0,"")

#mapf90("lvl", 19, "x",'./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_lvl_x_d0.02.bin',"",1,0,"")

#sys.exit()
#mapf90("mach", 15, "x",'./maps/high_res/velocity/15f15_analysis/map_high_15f15_x_map_vx_ew_Tsup7_5Mpc2.bin',"",1,1,"./papers_plots/paper_4_vel/mach_vx_ew.png")
#mapf90("vproj", 15, "z",'./maps/high_res/velocity/15f15_analysis/random_proj/map_high_15f15_rand_20_map_v_ew_Tsup7_5Mpc2.bin',"",1,0,"./papers_plots/paper_4_vel/vdx_ew.png")
#mapf90("d", 16, "x",'./maps/high_res/velocity/16f16_analysis/map_high_16f16_x_map_vdx_mw_Tsup7_5Mpc2.bin',"",1,1,"./papers_plots/paper_4_vel/vdx_mw.png")
#mapf90("d", 16, "y",'./maps/high_res/velocity/16f16_analysis/map_high_16f16_y_map_vdy_ew_Tsup7_5Mpc2.bin',"y",1,1,"./papers_plots/paper_4_vel/vdy_ew.png")
#mapf90("d", 16, "y",'./maps/high_res/velocity/16f16_analysis/map_high_16f16_y_map_vdy_mw_Tsup7_5Mpc2.bin',"",1,1,"./papers_plots/paper_4_vel/vdy_mw.png")
#mapf90("d", 16, "z",'./maps/high_res/velocity/16f16_analysis/map_high_16f16_z_map_vdz_ew_Tsup7_5Mpc2.bin',"z",1,1,"./papers_plots/paper_4_vel/vdz_ew.png")
#mapf90("d", 16, "z",'./maps/high_res/velocity/16f16_analysis/map_high_16f16_z_map_vdz_mw_Tsup7_5Mpc2.bin',"",1,1,"./papers_plots/paper_4_vel/vdz_mw.png")
#mapf90("d", 16, "z",'./maps/high_res/velocity/16f16_analysis/map_high_16f16_cen_map_vdz_ew_Tsup7_5Mpc2.bin',"cen",1,1,"./papers_plots/paper_4_vel/vdcen_ew.png")
#mapf90("d", 16, "z",'./maps/high_res/velocity/16f16_analysis/map_high_16f16_cen_map_vdz_mw_Tsup7_5Mpc2.bin',"",1,1,"./papers_plots/paper_4_vel/vdcen_mw.png")

#mapf90("vz", 16, "z",'./maps/high_res/velocity/map_high_16_z_map_vz.bin',"Sightline velocity, lvl 16",1)

#mapf90("vz", 16, "z",'./maps/high_res/velocity/map_high_16_z_map_vz_ew.bin',"Sightline velocity, emission-weighted, lvl 16",1)

#mapf90("vz", 17, "z",'./maps/high_res/velocity/map_high_17_z_map_vz.bin',"Sightline velocity, lvl 17",1)

#mapf90("vz", 17, "z",'./maps/high_res/velocity/map_high_17f17_z_map_vz.bin',"Sightline velocity, lvl 17 from sim 17",1)

#mapf90("vz", 17, "z",'./maps/high_res/velocity/map_high_17f17_z_map_vz_ew.bin',"Sightline velocity, emission-weighted, lvl 17 from sim 17",1)

#mapf90("vz", 19, "z",'./maps/high_res/velocity/map_high_19_z_core_map_vz.bin',"Sightline velocity, lvl 19",1)

#mapf90("vz", 21, "x",'./maps/high_res/velocity/map_high_21_x_core_map_vx.bin',"Sightline velocity, lvl 21",1)

#mapf90("d", 21, "z",'./maps/high_res/velocity/map_high_21_z_map_vdz_ew.bin',"Sightline velocity dispersion, lvl 21",1)

#mapf90("vz", 16, "x",'./maps/high_res/velocity/map_high_16f16_x_map_vx_ew.bin',"Sightline velocity, emission weighted lvl 16",1)

#mapf90("vz", 15, "z",'./maps/high_res/velocity/map_high_15f15_z_map_vz_ew_MWlos.bin',"EW velocity along MW sightline",1)

#mapf90("d", 16, "y",'./maps/high_res/velocity/map_high_16f16_analysis/map_high_16f16_y_map_vdy_mw_Tsup7_5Mpc2.bin',"MW vel disp y LOS, 5Mpc2",1)

#mapf90("d", 21, "z",'./maps/high_res/velocity/map_high_21_z_map_vdz_ew.bin',"Sightline velocity dispersion, lvl 21",1)

#mapf90("vz", 21, "z",'./maps/high_res/velocity/map_high_21_z_map_vz_ew_100kpc2_XRISM_T_range.bin',"EW velocity along z sightline, 100kpc**2 FoV, XRISM T range",1)

#mapf90("T", 21, "z",'./maps/high_res/velocity/map_high_21_z_map_T_mw_100kpc2.bin',"MW temperature along z sightline, 100kpc**2 FoV",1)

#mapf90("ne", 19, "z",'./maps/high_res/filament/map_high_19_xyz_full_map_ne_z_d0.02.bin',"Electron Density",1)

#mapf90("mach", 19, "z",'./maps/high_res/filament/map_high_19_xyz_full_map_mach_z_d0.02.bin',"Mach number",1)

#sys.exit()

#mapf90("vz", 16, "x",'./maps/high_res/velocity/map_high_16f16_x_map_vx_ew_5Mpc2.bin',"EW vx along x sightline",1)


#mapf90("T", 21, "z",'./maps/high_res/velocity/map_high_21_z_map_T_mw_100kpc2_10kpc_slice.bin',"MW temperature along z sightline, 10kpc slice",1)

#mapf90("ne", 21, "z",'./maps/high_res/velocity/map_high_21_z_map_ne_mw_100kpc2_10kpc_slice.bin',"MW Electron density along z sightline, 10kpc slice",1)


#mapf90("vz", 21, "z",'./maps/high_res/velocity/map_high_21_z_map_vz_ew_100kpc2.bin',"EW velocity along z sightline",1)


#sys.exit()

#mapf90("T", 19, "z",'./maps/high_res/filament/map_high_19_xyz_full_map_T_z_d0.02.bin',"Temperature",1)




#mapf90("ne", 19, "z",'./maps/high_res/filament/map_high_19_xyz_full_map_ne_z_d0.08.bin',"Electron Density",1)

#mapf90("P", 19, "z",'./maps/high_res/filament/map_high_19_xyz_full_map_P_z_d0.08.bin',"Pressure",1)

#'./maps/high_res/filament/map_high_19_xyz_full_map_ne_z_d0.08.bin'

#mapf90("v", 19, "z",'./maps/high_res/filament/map_high_19_xyz_full_map_vx_z_d0.02.bin',"Velocity field",1)

#mapf90("v", 15, "z",'./maps/high_res/filament/map_high_15_xyz_full_map_vx_z_d0.04.bin',"Velocity field",1)

#mapf90("v", 21, "z",'./maps/high_res/filament/map_high_21_xyz_left_map_vx_z_d0.08.bin',"V_x",1)

#mapf90("v", 19, "x",'./maps/high_res/filament/map_high_19_xyz_left_minus_0Mpc_map_vz_x_d0.02.bin',"Velocity field",2) #Velocity field

#mapf90("f", 19, "x",'./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_fz_x_d0.02.bin',"Flux field",1) #Flux field

###mapf90("ne-rhov", 19, "x",'./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_ne_x_d0.02.bin',"Electron density",1) #Electrond density

#mapf90("T", 19, "x",'./maps/high_res/filament/map_high_19_xyz_left_plus_2Mpc_map_T_x_d0.02.bin',"Velocity field",1) #Temperature

#mapf90("v", 21, "z",'maps/high_res/filament/map_high_21_xyz_mid3_bottom_right_map_vd_z_d0.08.bin',"v",1)

#mapf90("T", 19, "z",'./maps/high_res/filament/T_slices/map_high_21_xyz_left_map_T_z_d0.08.bin',"T",1)


#mapf90("P", 19, "y",'./maps/high_res/P_maps/map_high_19_y_P_los.bin',"Pressure",0)

#map_high_19_fil_large_map_ne_y.bin


#mapf90("SD", 19, "x",'./maps/high_res/map_high_19_cen_SD_bar_los.bin',"Total mass surface density",0)

#mapf90("P", 21, "z",'./maps/high_res/velocity/map_high_21_z_core_map_P.bin',"Pressure stacked along the z sightline",1)

#mapf90("P", 21, "z",'./maps/high_res/velocity/map_high_21_z_pm_500_kpc_map_P_d0.3.bin',"Pressure stacked along the z sightline",1)

#mapf90("v", 21, "z",'./maps/high_res/velocity/map_high_21_z_pm_500_kpc_map_vx_d0.02.bin',"vx",1)

#mapf90("v", 21, "z",'./maps/high_res/velocity/map_high_21_z_pm_1.3_Mpc_map_vz_d0.02.bin',"vz",1)

#mapf90("vz", 21, "y",'./maps/high_res/velocity/map_high_21_y_core_map_vy.bin',"Vy stacked along the y sightline",1)

#mapf90("vz", 21, "z",'./maps/high_res/velocity/map_high_21_z_core170kpc_4kev_cut_map_vz.bin',"Vz along the cen line of sight",1)

#sys.exit()


def map_3d_virgo_movie():
    d = FortranFile('virgo_xyz_dm_low_res.dat', 'r')

    print("ouverture du fichier dm")

    n= d.read_ints()
    x = d.read_reals()
    y = d.read_reals()
    z = d.read_reals()
    vx = d.read_reals()
    vy = d.read_reals()
    vz = d.read_reals()
    m = d.read_reals()

    unit_l = 0.227550518265197E+28

    x = x * (3.08567758128E21 / unit_l) + 0.5
    y = y * (3.08567758128E21 / unit_l) + 0.5
    z = z * (3.08567758128E21 / unit_l) + 0.5

    cond = np.logical_and(x > 0.45, np.logical_and(x < 0.53, np.logical_and(y > 0.47, np.logical_and(y < 0.55,
                                                                                                    np.logical_and(
                                                                                                        z > 0.46,
                                                                                                        z < 0.54)))))

    x = x[cond]
    y = y[cond]
    z = z[cond]

    #print("len(x)", len(x))

    step = 50
    n = int(len(x) / step)

    xplot = np.array([x[step * i] for i in range(n)])
    yplot = np.array([y[step * i] for i in range(n)])
    zplot = np.array([z[step * i] for i in range(n)])

    print("fin lecture fichier dm")

    Virgo = np.array([0.48461068, 0.50809848, 0.49687076])
    MW = np.array([0.5, 0.5, 0.5])
    fil= np.array([0.497,0.4984,0.4996])
    filx=np.array([0.4839,0.5087,0.51])
    fily = np.array([0.476,0.497,0.4953])
    opposfil=2*Virgo-fil
    opposMW=2*Virgo-MW
    opposfilx=2*Virgo-filx
    opposfily = 2 * Virgo - fily
    mainfil=np.array([opposfil,fil])
    vmw=np.array([opposMW,MW])
    vfilx=np.array([opposfilx,filx])
    vfily = np.array([opposfily, fily])
    #rint(test[:,0])
    #sys.exit()

    xl=0.492
    xu=0.496
    yl=0.501
    yu=0.504
    zl=0.496
    zu=0.502

    #section = "zoom"
    section = "rotation"


    for i in range(100):


        fig = plt.figure(figsize=(10,10))
        ax = plt.axes(projection='3d')
        #if section == "zoom":
        ax.set_axis_off()

        #cond=np.logical_and(xplot>xl,np.logical_and(xplot<xu,np.logical_and(yplot>yl,np.logical_and(yplot<yu,np.logical_and(zplot>zl,zplot<zu)))))
        #ax.scatter3D(xplot, yplot, zplot, s=2, c='black')
        #ax.scatter3D(Virgo[0], Virgo[1], Virgo[2], s=150, c='red')
        #ax.scatter3D(MW[0], MW[1], MW[2], s=80, c='blue')
        #ax.scatter3D(fil[0], fil[1], fil[2], s=100, c='red')
        ax.plot(mainfil[:, 0], mainfil[:, 1], mainfil[:, 2], color='blue', lw=8, label="Fil",ls="solid") #dotted
        ax.plot(vfilx[:, 0], vfilx[:, 1], vfilx[:, 2], color='green', lw=8, label="Filx",ls="solid") #dashed
        ax.plot(vfily[:, 0], vfily[:, 1], vfily[:, 2], color='dimgrey', lw=8, label="Fily",ls="solid") #dashdot
        #ax.plot(vmw[:, 0], vmw[:, 1], vmw[:, 2], color='red', lw=8, label="Cen",ls="solid")

        len_x = 0.4981 - 0.4711
        len_y = 0.52159 - 0.49459
        len_z = 0.5103 - 0.4833

        if section == "rotation":

            #ax.set_xlim(0.4711, 0.4981)
            #ax.set_ylim(0.49459, 0.52159)
            #ax.set_zlim(0.4833, 0.5103)

            zoom = 0.0125 * 19

            ax.set_xlim(0.4711 + len_x * zoom, 0.4981 - len_x * zoom)
            ax.set_ylim(0.49459 + len_y * zoom, 0.52159 - len_y * zoom)
            ax.set_zlim(0.4833 + len_z * zoom, 0.5103 - len_z * zoom)

        elif section == "zoom":

            zoom = 0.0125*i
            #lim=0.25
            ax.set_xlim(0.4711 + len_x * zoom, 0.4981 - len_x * zoom)
            ax.set_ylim(0.49459 + len_y * zoom, 0.52159 - len_y * zoom)
            ax.set_zlim(0.4833 + len_z * zoom, 0.5103 - len_z * zoom)

        #ax.set_xlim(0.4711 + len_x*0.1*i, 0.4981 - len_x*0.1*i)
        #ax.set_ylim(0.49459 + len_y*0.1*i, 0.52159 - len_y*0.1*i)
        #@ax.set_zlim(0.4833 + len_z*0.1*i, 0.5103 - len_z*0.1*i)


        #elev = 30 - 0.22*i #rotat 1
        #azim = -60 + 0.32*i #rotat 1

        elev = 8 + 0.9*i
        azim = -38

        if section == "rotation":
            ax.view_init(elev = elev, azim = azim)
        elif section == "zoom":
            ax.view_init(elev = 8, azim = -128)

        #RADIUS = 2.0  # Control this value.
        #ax.set_xlim3d(-RADIUS / 2, RADIUS / 2)
        #ax.set_zlim3d(-RADIUS / 2, RADIUS / 2)
        #ax.set_ylim3d(-RADIUS / 2, RADIUS / 2)

        #ax.set_xlim(xl,xu)
        #ax.set_ylim(yl, yu)
        #ax.set_zlim(zl, zu)

        # Second legend 'imaginary' lines
        line_dot = mlines.Line2D([], [], color='blue', linestyle='dotted', linewidth=4,label="Fil")
        line_dashed = mlines.Line2D([], [], color='green', linestyle='--',linewidth=4,label="Filx")
        line_dashdot = mlines.Line2D([], [], color='dimgrey', linestyle='dashdot', linewidth=4,label="Fily")
        line_solid = mlines.Line2D([], [], color='red', linestyle='-', linewidth=4,label='Cen')

        #ax.legend(handles=[line_dot,line_dashed,line_dashdot,line_solid])#,fontsize=20)
        if section == "rotation":
            ax.legend(handles=[line_dot, line_dashed, line_dashdot, line_solid])#,fontsize=20)

        j = 0
        j +=i

        if section == "rotation":
            if j<10:
                filename = 'Virgo_movie/im_rotat_4_00'+ str(j) + '.png'
            else:
                filename = 'Virgo_movie/im_rotat_4_0' + str(j) + '.png'

        elif section == "zoom" :
            if j<10:
                filename = 'Virgo_movie/im_zoom_4_00'+ str(j) + '.png'
            else:
                filename = 'Virgo_movie/im_zoom_4_0' + str(j) + '.png'

        plt.savefig(filename, format='png', bbox_inches='tight')
        if i==0:
            plt.show()
            sys.exit()
        print(i)
    #sys.exit()
        #plt.show()
    print("all figures saved")
    sys.exit()

def map_3d():
    f90_data = 1
    if f90_data==1:
        d = FortranFile('./virgo_xyz_files/virgo_xyz_dm_low_res.dat', 'r')
        #d = FortranFile('./virgo_xyz_files/virgo_xyz_dm_high_res_fil_los.dat', 'r')

        print("ouverture du fichier dm")

        n= d.read_ints()
        x = d.read_reals()
        y = d.read_reals()
        z = d.read_reals()
        vx = d.read_reals()
        vy = d.read_reals()
        vz = d.read_reals()
        m = d.read_reals()


    else:
        pos=np.load("dm_pos_low_ellip.npy")
        #print(np.shape(pos))
        #print(pos)
        x = pos[0,:]
        y = pos[1,:]
        z = pos[2,:]
        print("x",x)

        #sys.exit()

    unit_l = 0.227550518265197E+28

    x = x * (3.08567758128E21 / unit_l) + 0.5
    y = y * (3.08567758128E21 / unit_l) + 0.5
    z = z * (3.08567758128E21 / unit_l) + 0.5

    cond = np.logical_and(x > 0.45, np.logical_and(x < 0.53, np.logical_and(y > 0.47, np.logical_and(y < 0.55,np.logical_and(z > 0.46,z < 0.54)))))

    #lim = 5000
    lim = 10000
    #cond = np.logical_and(x > -lim, np.logical_and(x < lim, np.logical_and(y > -lim, np.logical_and(y < lim,np.logical_and(z > -lim,z < lim)))))




    x = x[cond]
    y = y[cond]
    z = z[cond]

    #print("len(x)", len(x))

    step = 50
    n = int(len(x) / step)

    xplot = np.array([x[step * i] for i in range(n)])
    yplot = np.array([y[step * i] for i in range(n)])
    zplot = np.array([z[step * i] for i in range(n)])

    print("fin lecture fichier dm")

    Virgo = np.array([0.48461068, 0.50809848, 0.49687076])
    MW = np.array([0.5, 0.5, 0.5])
    fil= np.array([0.497,0.4984,0.4996])
    filx=np.array([0.4839,0.5087,0.51])
    fily = np.array([0.476,0.497,0.4953])
    opposfil=2*Virgo-fil
    opposMW=2*Virgo-MW
    opposfilx=2*Virgo-filx
    opposfily = 2 * Virgo - fily
    mainfil=np.array([opposfil,fil])
    vmw=np.array([opposMW,MW])
    vfilx=np.array([opposfilx,filx])
    vfily = np.array([opposfily, fily])
    #rint(test[:,0])
    #sys.exit()

    xl = 0.492
    xu = 0.496
    yl = 0.501
    yu = 0.504
    zl = 0.496
    zu = 0.502

    x_axis_m = np.array([0.47,0.50809848, 0.49687076])
    x_axis_p = np.array([0.5,0.50809848, 0.49687076])
    vx=np.array([x_axis_p,x_axis_m])
    y_axis_m = np.array([0.48461068,0.495, 0.49687076])
    y_axis_p = np.array([0.48461068,0.525, 0.49687076])
    vy=np.array([y_axis_p,y_axis_m])
    z_axis_m = np.array([0.48461068, 0.50809848,0.48])
    z_axis_p = np.array([0.48461068, 0.50809848,0.51])
    vz=np.array([z_axis_p,z_axis_m])



    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection='3d')

    #cond=np.logical_and(xplot>xl,np.logical_and(xplot<xu,np.logical_and(yplot>yl,np.logical_and(yplot<yu,np.logical_and(zplot>zl,zplot<zu)))))
    ax.scatter3D(xplot, yplot, zplot, s=2, c='black',alpha = 0.4)
    #ax.scatter3D(Virgo[0], Virgo[1], Virgo[2], s=150, c='red')
    #ax.scatter3D(MW[0], MW[1], MW[2], s=80, c='blue')
    #ax.scatter3D(fil[0], fil[1], fil[2], s=100, c='red')
    #ax.plot(mainfil[:, 0], mainfil[:, 1], mainfil[:, 2], color='blue', lw=8, label="Fil",ls="dotted") #dotted
    #ax.plot(vfilx[:, 0], vfilx[:, 1], vfilx[:, 2], color='green', lw=8, label="Filx",ls="dashed") #dashed
    #ax.plot(vfily[:, 0], vfily[:, 1], vfily[:, 2], color='dimgrey', lw=8, label="Fily",ls="dashdot") #dashdot


    #line_dot = mlines.Line2D([], [], color='blue', linestyle='dotted', linewidth=4, label="Fil")
    #line_dashed = mlines.Line2D([], [], color='green', linestyle='--', linewidth=4, label="Filx")
    #line_dashdot = mlines.Line2D([], [], color='dimgrey', linestyle='dashdot', linewidth=4, label="Fily")
    #line_solid = mlines.Line2D([], [], color='red', linestyle='-', linewidth=4, label='Cen')

    ax.plot(vx[:, 0], vx[:, 1], vx[:, 2], color='red', lw=8, label="x", ls="solid")  # dotted
    ax.plot(vy[:, 0], vy[:, 1], vy[:, 2], color='blue', lw=8, label="y", ls="solid")  # dashed
    ax.plot(vz[:, 0], vz[:, 1], vz[:, 2], color='green', lw=8, label="z", ls="solid")  # dashdot

    ax.plot(vmw[:, 0], vmw[:, 1], vmw[:, 2], color='deeppink', lw=8, label="Cen", ls="solid")

    #line_dot = mlines.Line2D([], [], color='red', linestyle='solid', linewidth=4, label="x")
    #line_dashed = mlines.Line2D([], [], color='blue', linestyle='solid', linewidth=4, label="y")
    #line_dashdot = mlines.Line2D([], [], color='green', linestyle='solid', linewidth=4, label="z")
    #line_solid = mlines.Line2D([], [], color='mediumvioletred', linestyle='solid', linewidth=4, label='Cen')

    #len_x = 0.4981 - 0.4711
    #len_y = 0.52159 - 0.49459
    #len_z = 0.5103 - 0.4833


    ax.set_xlim(0.4711, 0.4981)
    ax.set_ylim(0.49459, 0.52159)
    ax.set_zlim(0.4833, 0.5103)

    #ax.set_xlim(-lim,lim)
    #ax.set_ylim(-lim,lim)
    #ax.set_zlim(-lim,lim)




    #ax.legend(handles=[line_dot, line_dashed, line_dashdot, line_solid])#,fontsize=20)
    ax.legend()
    ax.set_xlabel("x [Mpc]")
    ax.set_ylabel("y [Mpc]")
    ax.set_zlabel("z [Mpc]")

    ticks_box=np.linspace(0.471,0.4981,5)
    ticks=np.linspace(-10,10,5)
    labels_mpc=[str(i) for i in ticks]
    #labels[5]="0"
    #labels[0]="-10"
    #labels[-1]="10"
    ax.set_xticks(ticks_box)
    ax.set_xticklabels(labels_mpc)

    ticks_box_y = np.linspace(0.4944, 0.5215, 5)
    ax.set_yticks(ticks_box_y)
    ax.set_yticklabels(labels_mpc)

    ticks_box_z = np.linspace(0.4833, 0.5104, 5)
    ax.set_zticks(ticks_box_z)
    ax.set_zticklabels(labels_mpc)

    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

    def plot_ellipsoid():

        fact = 2787.5

        a = (782149.68011333 / 782149.68011333) * fact
        b = (432041.40085049 / 782149.68011333) * fact
        c = (518642.86421945 / 782149.68011333) * fact

        rx = 1 / a
        ry = 1 / b
        rz = 1 / c

        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        x = a * np.outer(np.cos(u), np.sin(v))
        y = b * np.outer(np.sin(u), np.sin(v))
        z = c * np.outer(np.ones_like(u), np.cos(v))

        ax.plot_surface(x, y, z, rstride=4, cstride=4, color='b', alpha=0.5)

    #plot_ellipsoid()

    plt.show()
    sys.exit()

map_3d()
sys.exit()

def map_4_p():

    def show_radii_a():
        xpx = np.linspace(-10,10,4000)
        ypx = np.linspace(-10,10,4000)

        Xpx, Ypx = np.meshgrid(xpx, ypx)

        # F = (Xpx-gal_pos[i,0])**2 + (Ypx-gal_pos[i,1])**2 - (r200dm[i]*(10/1e4))**2
        #F = (Xpx - x[i]) ** 2 + (Ypx - y[i]) ** 2 - (r200dm[i] * (10 / 1e4)) ** 2

        F = Xpx**2 + Ypx**2 - (2024*(10/1e4))**2  #rvir
        rvir=a.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9)
        a.clabel(rvir,rvir.levels,inline=True,fmt="$R_{vir}$",fontsize=20)

        F = Xpx ** 2 + Ypx ** 2 - (5031 * (10 / 1e4)) ** 2  #rzv
        rzv=a.contour(Xpx, Ypx, F, [0], colors='white',ls='dashed', linewidths=2, alpha=0.9)
        a.clabel(rzv, rzv.levels, inline=True, fmt="$R_{zv}$", fontsize=20)

    def show_radii_b():
        xpx = np.linspace(-10,10,4000)
        ypx = np.linspace(-10,10,4000)

        Xpx, Ypx = np.meshgrid(xpx, ypx)

        # F = (Xpx-gal_pos[i,0])**2 + (Ypx-gal_pos[i,1])**2 - (r200dm[i]*(10/1e4))**2
        #F = (Xpx - x[i]) ** 2 + (Ypx - y[i]) ** 2 - (r200dm[i] * (10 / 1e4)) ** 2

        F = Xpx**2 + Ypx**2 - (2024*(10/1e4))**2  #rvir
        rvir=b.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9)
        b.clabel(rvir,rvir.levels,inline=True,fmt="$R_{vir}$",fontsize=20)

        F = Xpx ** 2 + Ypx ** 2 - (5031 * (10 / 1e4)) ** 2  #rzv
        rzv=b.contour(Xpx, Ypx, F, [0], colors='white',ls='dashed', linewidths=2, alpha=0.9)
        b.clabel(rzv, rzv.levels, inline=True, fmt="$R_{zv}$", fontsize=20)

    def show_radii_c():
        xpx = np.linspace(-10,10,4000)
        ypx = np.linspace(-10,10,4000)

        Xpx, Ypx = np.meshgrid(xpx, ypx)

        # F = (Xpx-gal_pos[i,0])**2 + (Ypx-gal_pos[i,1])**2 - (r200dm[i]*(10/1e4))**2
        #F = (Xpx - x[i]) ** 2 + (Ypx - y[i]) ** 2 - (r200dm[i] * (10 / 1e4)) ** 2

        F = Xpx**2 + Ypx**2 - (2024*(10/1e4))**2  #rvir
        rvir=c.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9)
        c.clabel(rvir,rvir.levels,inline=True,fmt="$R_{vir}$",fontsize=20)

        F = Xpx ** 2 + Ypx ** 2 - (5031 * (10 / 1e4)) ** 2  #rzv
        rzv=c.contour(Xpx, Ypx, F, [0], colors='white',ls='dashed', linewidths=2, alpha=0.9)
        c.clabel(rzv, rzv.levels, inline=True, fmt="$R_{zv}$", fontsize=20)

    def show_radii_d():
        xpx = np.linspace(-10,10,4000)
        ypx = np.linspace(-10,10,4000)

        Xpx, Ypx = np.meshgrid(xpx, ypx)

        # F = (Xpx-gal_pos[i,0])**2 + (Ypx-gal_pos[i,1])**2 - (r200dm[i]*(10/1e4))**2
        #F = (Xpx - x[i]) ** 2 + (Ypx - y[i]) ** 2 - (r200dm[i] * (10 / 1e4)) ** 2

        F = Xpx**2 + Ypx**2 - (2024*(10/1e4))**2  #rvir
        rvir=d.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9)
        d.clabel(rvir,rvir.levels,inline=True,fmt="$R_{vir}$",fontsize=20)

        F = Xpx ** 2 + Ypx ** 2 - (5031 * (10 / 1e4)) ** 2  #rzv
        rzv=d.contour(Xpx, Ypx, F, [0], colors='white',ls='dashed', linewidths=2, alpha=0.9)
        d.clabel(rzv, rzv.levels, inline=True, fmt="$R_{zv}$", fontsize=20)

    hcen = FortranFile('./maps/high_res/P_maps/map_high_19_cen_P_los.bin', 'r')
    map_cen = hcen.read_reals()
    px = np.int_(np.sqrt(len(map_cen)))
    print("nbr px per line", px)
    map_cen = np.reshape(map_cen, (px, px))
    map_cen = np.log10(map_cen)
    map_cen[np.isnan(map_cen) == True] = -8

    hfil = FortranFile('./maps/high_res/P_maps/map_high_19_fil_P_los.bin', 'r')
    map_fil = hfil.read_reals()
    map_fil = np.reshape(map_fil, (px, px))
    map_fil = np.log10(map_fil)
    map_fil[np.isnan(map_fil) == True] = -8

    hfilx = FortranFile('./maps/high_res/P_maps/map_high_19_filx_P_los.bin', 'r')
    map_filx = hfilx.read_reals()
    map_filx = np.reshape(map_filx, (px, px))
    map_filx = np.log10(map_filx)
    map_filx[np.isnan(map_filx) == True] = -8

    hfily = FortranFile('./maps/high_res/P_maps/map_high_19_fily_P_los.bin', 'r')
    map_fily = hfily.read_reals()
    map_fily = np.reshape(map_fily, (px, px))
    map_fily = np.log10(map_fily)
    map_fily[np.isnan(map_fily) == True] = -8

    dim=11.0615

    dim = [-dim, dim, -dim, dim]

    fontprops = fm.FontProperties(size=8)

    #scalebar = AnchoredSizeBar(a.transData,2, '2 Mpc', 'lower center',pad=0.1,color='white',frameon=False,size_vertical=1,fontproperties=fontprops)
    #a.add_artist(scalebar)

    print("maps loaded")

    fig, ((a,c),(b,d)) = plt.subplots(2,2)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])


    im_cen = a.imshow(map_cen, cmap="inferno", origin='lower', alpha=1, extent=dim, vmin=-8, vmax=-1.5)
    a.yaxis.set_ticklabels([])
    a.xaxis.set_ticklabels([])
    a.axis("off")
    scalebar = AnchoredSizeBar(a.transData, 5, '5 Mpc', 'lower left', pad=0.1, color='white', frameon=False,size_vertical=0.2, fontproperties=fontprops,label_top=True)
    a.add_artist(scalebar)
    a.set_title("Cen")
    a.margins(x=0)
    show_radii_a()

    im_fil = b.imshow(map_fil, cmap="inferno", origin='lower', alpha=1, extent=dim, vmin=-8, vmax=-1.5)
    b.yaxis.set_ticklabels([])
    b.xaxis.set_ticklabels([])
    b.axis("off")
    scalebar = AnchoredSizeBar(b.transData, 5, '5 Mpc', 'lower left', pad=0.1, color='white', frameon=False,size_vertical=0.2, fontproperties=fontprops,label_top=True)
    b.add_artist(scalebar)
    #b.set_title("Fil")
    b.set_xlabel("Fil")
    show_radii_b()

    im_filx = c.imshow(map_filx, cmap="inferno", origin='lower', alpha=1, extent=dim, vmin=-8, vmax=-1.5)
    c.yaxis.set_ticklabels([])
    c.xaxis.set_ticklabels([])
    c.axis("off")
    scalebar = AnchoredSizeBar(c.transData, 5, '5 Mpc', 'lower left', pad=0.1, color='white', frameon=False,size_vertical=0.2, fontproperties=fontprops,label_top=True)
    c.add_artist(scalebar)
    c.set_title("Filx")
    show_radii_c()

    im_fily = d.imshow(map_fily, cmap="inferno", origin='lower', alpha=1, extent=dim, vmin=-8, vmax=-1.5)
    d.yaxis.set_ticklabels([])
    d.set_xlabel("Fily")
    d.xaxis.set_ticklabels([])
    d.axis("off")
    scalebar = AnchoredSizeBar(d.transData, 5, '5 Mpc', 'lower left', pad=0.1, color='white', frameon=False,size_vertical=0.2, fontproperties=fontprops,label_top=True)
    d.add_artist(scalebar)
    #d.set_title("Fily")
    show_radii_d()

    cb=fig.colorbar(im_fily, cax=cbar_ax)

    #cb=fig.colorbar(im_fily,ax=d)
    cb.set_label('$log_{10}(P[keV.cm^{-3}])$',size='large')

    #fig.tight_layout(pad=0)

    #plt.subplots_adjust(wspace=0.0001, hspace=0.01)

    plt.show()

def sphere_3d():
    # from mpl_toolkits.mplot3d import Axes3D  # Not needed with Matplotlib 3.6.3
    #import matplotlib.pyplot as plt
    #import numpy as np

    fig = plt.figure(figsize=plt.figaspect(1))  # Square figure
    ax = fig.add_subplot(111, projection='3d')

    #a=(-3)**2
    #b=70**2
    #c=133*2

    #a=4
    #b=2
    #c=1

    #a=7.33
    #b=3.19
    #c=5.06

    #a = 7.82
    #b = 4.32
    #c = 5.18

    fact = 2787.5

    a = (782149.68011333 / 782149.68011333) * fact
    b = (432041.40085049 / 782149.68011333) * fact
    c = (518642.86421945 / 782149.68011333) * fact


    print("coeffs",a,b,c)

    #a = a * (3.08567758128E21 / unit_l)
    #b = b * (3.08567758128E21 / unit_l)
    #c = c * (3.08567758128E21 / unit_l)

    print("coeffs",a,b,c)

    #sys.exit()


    coefs = (a,b,c)  # Coefficients in a0/c x**2 + a1/c y**2 + a2/c z**2 = 1
    # Radii corresponding to the coefficients:

    #print("rx,ry,rz",rx,ry,rz)

    #rx, ry, rz = 1 / np.sqrt(coefs)
    rx = 1 / a
    ry = 1 / b
    rz = 1 / c

    print("rx,ry,rz", rx, ry, rz)

    print("sqrt",np.sqrt(coefs))

    #Virgo = np.array([0.48461068, 0.50809848, 0.49687076]) #before rotations
    Virgo = np.array([0.48533832, 0.49623066, 0.4908876])



    #sys.exit()

    # Set of all spherical angles:
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    # Cartesian coordinates that correspond to the spherical angles:
    # (this is the equation of an ellipsoid):
    #x = Virgo[0] + a * np.outer(np.cos(u), np.sin(v))
    #y = Virgo[1] + b * np.outer(np.sin(u), np.sin(v))
    #z = Virgo[2] + c * np.outer(np.ones_like(u), np.cos(v))

    x = a * np.outer(np.cos(u), np.sin(v))
    y = b * np.outer(np.sin(u), np.sin(v))
    z = c * np.outer(np.ones_like(u), np.cos(v))

    #x = a * np.outer(np.cos(v), np.cos(u))
    #y = b * np.outer(np.cos(v), np.sin(u))
    #z = c * np.outer(np.ones_like(v), np.sin(v))

    ax.plot_surface(x, y, z, rstride=4, cstride=4, color='b', alpha=0.1)

    x=np.zeros(10000)
    y=np.zeros(10000)
    z=np.zeros(10000)

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(-0.5*np.pi, 0.5*np.pi, 100)

    for i in range(100):
        for j in range(100):
            k=100*i+j
            x[k] = a * np.cos(v[i]) * np.cos(u[j])
            y[k] = b * np.cos(v[i]) * np.sin(u[j])
            z[k] = c * np.sin(v[i])

    #ax.scatter3D(x,y,z,s=1,alpha=0.5,c='green')

    #plt.show()
    #sys.exit()



    # Plot:


    # Adjustment of the axes, so that they all have the same span:
    max_radius = max(a, b, c)
    for axis in 'xyz':
        getattr(ax, 'set_{}lim'.format(axis))((-max_radius, max_radius))

    test=np.random.rand(3,100)
    test=test*2-1
    xt=test[0,:]*a
    yt=test[1,:]*b
    zt=test[2,:]*c

    #print("xt",xt)
    #print("yt", yt)
    #print("zt", zt)

    at,bt,ct=np.sqrt(coefs)

    ell=(xt/a)**2+(yt/b)**2+(zt/c)**2

    #print("ell",ell)

    #ax.scatter3D(xt[ell>1],yt[ell>1],zt[ell>1],c='red')
    #ax.scatter3D(xt[ell < 1], yt[ell < 1], zt[ell < 1], c='black')

    plt.show()

def plot_div_v():
    #div_v = np.load("./filament/left_fil_vel_ops/div_v_meanz.npy")

    #div_v = np.load("./filament/left_fil_transverse_vel_ops/div_v_minus_2Mpc.npy")

    div_v = np.load("./filament/right_fil_transverse_vel_ops/div_v_plus_2Mpc.npy")

    #div_v = np.load("./filament/left_fil_transverse_vel_ops/div_f_minus_8Mpc.npy")

    #div_v = np.load("./filament/div_v_full_map.npy")

    #div_v = np.load("./velocity/div_v_core_map_z_proj.npy")



    print("min div_v",np.min(div_v))
    print("max div_v",np.max(div_v))

    #sys.exit()

    #plt.hist(div_v[div_v>0].flatten(), bins=100, log=True)
    #plt.xscale('log')
    #plt.show()
    #sys.exit()


    div_v = div_v.T
    div_v = np.flip(div_v, axis=0)

    nx, ny = div_v.shape
    lvl=19

    dimx = (nx / 2) * (737.441 / 2 ** lvl)
    dimy = (ny / 2) * (737.441 / 2 ** lvl)
    dim = [-dimx, dimx, -dimy, dimy]

    fig, ax = plt.subplots(facecolor='white',figsize=(12,12))

    #im = ax.imshow(div_v, cmap="PuOr", norm=colors.SymLogNorm(linthresh=1e-25, linscale=0.3,vmin=-1e-22, vmax=1e-22, base=10), extent=dim)

    #im = ax.imshow(np.abs(div_v), cmap="Reds", norm=colors.LogNorm(vmin=1e-24, vmax=1e-22), extent=dim)

    im = ax.imshow(div_v, cmap="seismic", vmin=-20, vmax=20,extent=dim)  # , origin='lower')

    #plt.title("2D velocity field divergence")
    #plt.title("2D flux field divergence")

    #cb = fig.colorbar(im, ax=ax)
    ax.set_xlabel("x [Mpc]")
    ax.set_ylabel("y [Mpc]")

    #cb.set_label(r'$\vec{\nabla}.\vec{v}~[km.s^{-1}.kpc^{-1}]$', size='large')

    #cb.set_label(r'$\vec{\nabla}.\rho\vec{v}[kg~m^{-2}~s^{-1}~kpc^{-1}]$', size='large')

    plt.grid(b=None)

    ax.get_xaxis().set_visible(False)
    #ax.get_yaxis().set_visible(False)

    plt.savefig('./papers_plots/paper_3_turb_fil/div_v_plus2.png', dpi=600, format='png')



    print("plot saved")

    plt.show()

def plot_rot_v():
    #rot_v = np.load("./filament/left_fil_vel_ops/rot_vz_meanz.npy")
    #rot_v = np.load("./filament/left_fil_transverse_vel_ops/rot_vx_minus_6Mpc.npy")

    rot_v = np.load("./filament/right_fil_transverse_vel_ops/rot_vx_plus_2Mpc.npy")

    #rot_v = np.load("./filament/rot_vz_full_map.npy")

    #rot_v = np.load("./velocity/rot_vz_core_map_z_proj.npy")

    rot_v = rot_v.T
    rot_v = np.flip(rot_v, axis=0)

    nx, ny = rot_v.shape
    lvl = 19

    dimx = (nx / 2) * (737.441 / 2 ** lvl)
    dimy = (ny / 2) * (737.441 / 2 ** lvl)
    dim = [-dimx, dimx, -dimy, dimy]

    fig, ax = plt.subplots(facecolor='white',figsize=(10,10))

    im = ax.imshow(rot_v, cmap="seismic", vmin=-20, vmax=20,extent=dim)  # , origin='lower')
    #plt.title("Divergence of velocity (mean over z)")

    cb = fig.colorbar(im, ax=ax)
    ax.set_xlabel("x [Mpc]")
    ax.set_ylabel("y [Mpc]")

    cb.set_label(r'$\omega_x~[km.s^{-1}.kpc^{-1}]$', size='large')

    plt.title("Vorticity in x direction")

    plt.grid(b=None)

    plt.savefig('./papers_plots/paper_3_turb_fil/w_x_plus2.png', dpi=300, format='png')

    print("plot saved")

    plt.show()

#plot_div_v()
#plot_rot_v()

def large_fil_plots():


    # load files

    #print("loading maps")

    #t8 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_T_x_d0.02.bin','x')
    #t6 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_6Mpc_map_T_x_d0.02.bin','x')
    #t4 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_4Mpc_map_T_x_d0.02.bin','x')
    #t2 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_2Mpc_map_T_x_d0.02.bin','x')
    #t0 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_0Mpc_map_T_x_d0.02.bin','x')

    #print("T maps loaded")

    def load_map(file, proj):
        h = FortranFile(file, 'r')

        nx, ny, nz = h.read_ints()
        cen_x, cen_y, cen_z = h.read_reals()

        if proj == "x":
            ncell = nz * ny
        elif proj == "y":
            ncell = nx * nz
        elif proj == "z":
            ncell = nx * ny

        map = np.zeros(ncell)

        map = ftp.f90_to_py.read_map_file(ncell, file, 0)

        if proj == "x":
            map = np.reshape(map, (nz, ny))
            #map2 = np.reshape(map2, (nz, ny))
            #map3 = np.reshape(map3, (nz, ny))
            nl = nx
        elif proj == "y":
            map = np.reshape(map, (nx, nz))
            #map2 = np.reshape(map2, (nx, nz))
            #map3 = np.reshape(map3, (nx, nz))
            nl = ny
        elif proj == "z":
            map = np.reshape(map, (ny, nx))
            #map2 = np.reshape(map2, (ny, nx))
            #map3 = np.reshape(map3, (ny, nx))
            nl = nz

        return map

    def left_fil_T_plot():
        print("loading maps")

        t8 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_T_x_d0.02.bin', 'x')
        t6 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_6Mpc_map_T_x_d0.02.bin', 'x')
        t4 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_4Mpc_map_T_x_d0.02.bin', 'x')
        t2 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_2Mpc_map_T_x_d0.02.bin', 'x')
        t0 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_0Mpc_map_T_x_d0.02.bin', 'x')

        print("maps loaded")

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        f, axs = plt.subplots(1, 5, figsize=(16, 4), constrained_layout=True, facecolor='white')
        # gs = gridspec.GridSpec(4, 4, wspace=0.001, hspace=0)

        fs = 14

        title = "$\Delta x_{cen}$= -8Mpc"

        # First row : temperature

        # axs[0, 0].plot(x, np.sin(x))
        # ax = fig.add_subplot(gs[0, 0])
        # plt.sca(axs[0, 0])
        plt.sca(axs[0])
        title = "$\Delta x_{cen}$= -8Mpc"
        plt.title(title, fontsize=fs)
        t8 /= (kb / 1.602e-16)
        plt.imshow(t8, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        plt.xticks([])
        plt.ylabel("y [Mpc]", fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 1])
        # plt.sca(axs[0, 1])
        plt.sca(axs[1])
        title = "$\Delta x_{cen}$= -6Mpc"
        plt.title(title, fontsize=fs)
        t6 /= (kb / 1.602e-16)
        plt.imshow(t6, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 2])
        # plt.sca(axs[0, 2])
        plt.sca(axs[2])
        title = "$\Delta x_{cen}$= -4Mpc"
        plt.title(title, fontsize=fs)
        t4 /= (kb / 1.602e-16)
        plt.imshow(t4, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 3])
        plt.sca(axs[3])
        title = "$\Delta x_{cen}$= -2Mpc"
        plt.title(title, fontsize=fs)
        t2 /= (kb / 1.602e-16)
        im = plt.imshow(t2, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.14)  # , orientation="horizontal")
        # cb.set_label('T[K]', size='large')
        plt.sca(axs[4])
        title = "$\Delta x_{cen}$= 0Mpc"
        plt.title(title, fontsize=fs)
        t0 /= (kb / 1.602e-16)
        im = plt.imshow(t0, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)

        # ratiov = t0.shape[0] / t0.shape[1]

        # print("ratiov",ratiov)

        # sys.exit()

        cb = f.colorbar(im, ax=axs.ravel().tolist(), shrink=0.6, pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)
        cb.set_label('$T~\mathrm{[K]}$', size=fs+2)

        plt.savefig('./papers_plots/paper_3_turb_fil/T_filament_left.png', dpi=600, format='png')

        print("plot saved")

        plt.show()

        sys.exit()

    def left_fil_v_plot():

        vz8 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_vz_x_d0.02.bin', 'x')
        vz6 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_6Mpc_map_vz_x_d0.02.bin', 'x')
        vz4 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_4Mpc_map_vz_x_d0.02.bin', 'x')
        vz2 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_2Mpc_map_vz_x_d0.02.bin', 'x')
        vz0 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_0Mpc_map_vz_x_d0.02.bin', 'x')

        print("vz maps loaded")

        vz_virgo = -131.9249

        vz8 -= vz_virgo
        vz6 -= vz_virgo
        vz4 -= vz_virgo
        vz2 -= vz_virgo
        vz0 -= vz_virgo

        mean10_map_8 = np.load("stream_vz_left_fil_minus_8Mpc.npy")
        mean10_map2_8 = np.load("stream_vy_left_fil_minus_8Mpc.npy")
        mean10_map_6 = np.load("stream_vz_left_fil_minus_6Mpc.npy")
        mean10_map2_6 = np.load("stream_vy_left_fil_minus_6Mpc.npy")
        mean10_map_4 = np.load("stream_vz_left_fil_minus_4Mpc.npy")
        mean10_map2_4 = np.load("stream_vy_left_fil_minus_4Mpc.npy")
        mean10_map_2 = np.load("stream_vz_left_fil_minus_2Mpc.npy")
        mean10_map2_2 = np.load("stream_vy_left_fil_minus_2Mpc.npy")
        mean10_map_0 = np.load("stream_vz_left_fil_minus_0Mpc.npy")
        mean10_map0_2 = np.load("stream_vy_left_fil_minus_0Mpc.npy")

        x = np.load("mesh_x.npy")
        y = np.load("mesh_y.npy")

        print("streamlines loaded")

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        # Create a figure and a 4x4 grid of subplots
        # fig, axs = plt.subplots(4, 4, figsize=(16, 16), facecolor='white')
        # fig = plt.figure(figsize=(16, 16), facecolor='white', constrained_layout=True)

        # Second row : vz + streams

        f, axs = plt.subplots(1, 5, figsize=(16, 4), constrained_layout=True, facecolor='white')

        fs = 14

        map_mean_mesh_8 = np.sqrt(mean10_map_8 ** 2 + mean10_map2_8 ** 2)
        map_mean_mesh_8[np.isnan(map_mean_mesh_8) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 0])
        plt.sca(axs[0])

        stream = plt.streamplot(x, y, mean10_map_8, mean10_map2_8, linewidth=0.3, color=map_mean_mesh_8, cmap='Greys_r',
                                density=[1, 1.5], norm=norm)
        im = plt.imshow(vz8, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks([])
        plt.ylabel("y [Mpc]", fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)

        map_mean_mesh_6 = np.sqrt(mean10_map_6 ** 2 + mean10_map2_6 ** 2)
        map_mean_mesh_6[np.isnan(map_mean_mesh_6) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 1])

        plt.sca(axs[1])

        stream = plt.streamplot(x, y, mean10_map_6, mean10_map2_6, linewidth=0.4, color=map_mean_mesh_6, cmap='Greys_r',
                                density=[1, 1.5], norm=norm)
        im = plt.imshow(vz6, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        map_mean_mesh_4 = np.sqrt(mean10_map_4 ** 2 + mean10_map2_4 ** 2)
        map_mean_mesh_4[np.isnan(map_mean_mesh_4) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 2])
        plt.sca(axs[2])

        stream = plt.streamplot(x, y, mean10_map_4, mean10_map2_4, linewidth=0.4, color=map_mean_mesh_4, cmap='Greys_r',
                                density=[1, 1.5], norm=norm)
        im = plt.imshow(vz4, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        map_mean_mesh_2 = np.sqrt(mean10_map_2 ** 2 + mean10_map2_2 ** 2)
        map_mean_mesh_2[np.isnan(map_mean_mesh_2) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 3])
        plt.sca(axs[3])

        stream = plt.streamplot(x, y, mean10_map_2, mean10_map2_2, linewidth=0.4, color=map_mean_mesh_2, cmap='Greys_r',
                                density=[1, 1.5], norm=norm)
        im = plt.imshow(vz2, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)

        plt.sca(axs[4])

        stream = plt.streamplot(x, y, mean10_map_0, mean10_map0_2, linewidth=0.4, color=map_mean_mesh_2, cmap='Greys_r',
                                density=[1, 1.5], norm=norm)
        im = plt.imshow(vz0, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # ax.get_xaxis().set_visible(False)
        # ax.get_yaxis().set_visible(False)

        # cb = fig.colorbar(im, ax=ax,pad = 0.14, fraction=0.03 )  # , orientation="horizontal")
        # cb.set_label('vz[km/s]', size='large')

        # cbs = fig.colorbar(stream.lines, ax=ax, orientation='vertical', pad=0.24, fraction=0.03)
        # cbs.set_label(r'$||\vec{v}||~[km/s]$', size='large')

        #cb = f.colorbar(im, ax=axs.ravel().tolist(), label='$v_z~\mathrm{[km~s^{-1}]}$', shrink=0.6, pad=0.01)  # ,size='small')#,shrink=0.6)
        #cb.ax.tick_params(labelsize=fs)

        cb = f.colorbar(stream.lines, ax=axs.ravel().tolist(), label=r'$||\vec{v}||~[km/s]$', shrink=0.6, pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/v_filament_left_2.png', dpi=600, format='png')

        print("plot saved")

        plt.show()
        sys.exit()

    def left_fil_div_plot():

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        div_v_8 = np.load("./filament/left_fil_transverse_vel_ops/div_v_minus_8Mpc.npy")
        div_v_6 = np.load("./filament/left_fil_transverse_vel_ops/div_v_minus_6Mpc.npy")
        div_v_4 = np.load("./filament/left_fil_transverse_vel_ops/div_v_minus_4Mpc.npy")
        div_v_2 = np.load("./filament/left_fil_transverse_vel_ops/div_v_minus_2Mpc.npy")
        div_v_0 = np.load("./filament/left_fil_transverse_vel_ops/div_v_minus_0Mpc.npy")

        print("div_v maps loaded")

        div_v_8 = div_v_8.T
        div_v_8 = np.flip(div_v_8, axis=0)
        div_v_6 = div_v_6.T
        div_v_6 = np.flip(div_v_6, axis=0)
        div_v_4 = div_v_4.T
        div_v_4 = np.flip(div_v_4, axis=0)
        div_v_2 = div_v_2.T
        div_v_2 = np.flip(div_v_2, axis=0)
        div_v_0 = div_v_0.T
        div_v_0 = np.flip(div_v_0, axis=0)


        #Third row : div_v

        #ax = fig.add_subplot(gs[2, 0])

        f, axs = plt.subplots(1, 5, figsize=(16, 4), constrained_layout=True, facecolor='white')

        fs = 14

        plt.sca(axs[0])

        im = plt.imshow(div_v_8, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks([])
        plt.ylabel("y [Mpc]",fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(b=None)
        #plt.get_xaxis().set_visible(False)

        #ax = fig.add_subplot(gs[2, 1])

        plt.sca(axs[1])

        im = plt.imshow(div_v_6, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        #plt.get_xaxis().set_visible(False)
        #plt.get_yaxis().set_visible(False)

        #ax = fig.add_subplot(gs[2, 2])

        plt.sca(axs[2])

        im = plt.imshow(div_v_4, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        #ax.get_xaxis().set_visible(False)
        #ax.get_yaxis().set_visible(False)

        #ax = fig.add_subplot(gs[2, 3])

        plt.sca(axs[3])

        im = plt.imshow(div_v_2, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks([])
        plt.yticks([])

        plt.grid(b=None)
        #plt.get_xaxis().set_visible(False)
        #plt.get_yaxis().set_visible(False)


        #cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.14)#,pad = 0.04 )  # , orientation="horizontal")
        #cb.set_label(r'$\vec{\nabla}.\vec{v}[km.s^{-1}.kpc^{-1}]$', size='large')

        plt.sca(axs[4])

        im = plt.imshow(div_v_0, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks([])
        plt.yticks([])

        plt.grid(b=None)

        cb = f.colorbar(im, ax=axs.ravel().tolist(), label=r'$\vec{\nabla}.\vec{v}~[km~s^{-1}~kpc^{-1}]$', shrink=0.6,
                            pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/div_v_filament_left.png', dpi=600, format='png')

        print('fig saved')

        #Fourth row : rot_v

        #ax = fig.add_subplot(gs[3, 0])

        plt.show()
        sys.exit()

    def left_fil_rot_v_plot():

        rot_v_8 = np.load("./filament/left_fil_transverse_vel_ops/rot_vx_minus_8Mpc.npy")
        rot_v_6 = np.load("./filament/left_fil_transverse_vel_ops/rot_vx_minus_6Mpc.npy")
        rot_v_4 = np.load("./filament/left_fil_transverse_vel_ops/rot_vx_minus_4Mpc.npy")
        rot_v_2 = np.load("./filament/left_fil_transverse_vel_ops/rot_vx_minus_2Mpc.npy")
        rot_v_0 = np.load("./filament/left_fil_transverse_vel_ops/rot_vx_minus_0Mpc.npy")

        print("rot_v maps loaded")

        rot_v_8 = rot_v_8.T
        rot_v_8 = np.flip(rot_v_8, axis=0)
        rot_v_6 = rot_v_6.T
        rot_v_6 = np.flip(rot_v_6, axis=0)
        rot_v_4 = rot_v_4.T
        rot_v_4 = np.flip(rot_v_4, axis=0)
        rot_v_2 = rot_v_2.T
        rot_v_2 = np.flip(rot_v_2, axis=0)
        rot_v_0 = rot_v_0.T
        rot_v_0 = np.flip(rot_v_0, axis=0)


        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        fs = 14

        f, axs = plt.subplots(1, 5, figsize=(16, 4), constrained_layout=True, facecolor='white')

        plt.sca(axs[0])


        im = plt.imshow(rot_v_8, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xlabel("z [Mpc]",fontsize=fs)
        plt.ylabel("y [Mpc]", fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.grid(b=None)

        #ax = fig.add_subplot(gs[3, 1])

        plt.sca(axs[1])

        im = plt.imshow(rot_v_6, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xlabel("z [Mpc]",fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        #plt.get_yaxis().set_visible(False)
        plt.xticks(fontsize=fs)

        #ax = fig.add_subplot(gs[3, 2])

        plt.sca(axs[2])

        im = plt.imshow(rot_v_4, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xlabel("z [Mpc]",fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        #plt.get_yaxis().set_visible(False)
        plt.xticks(fontsize=fs)

        #ax = fig.add_subplot(gs[3, 3])

        plt.sca(axs[3])

        im = plt.imshow(rot_v_2, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xlabel("z [Mpc]",fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        #plt.get_yaxis().set_visible(False)
        plt.xticks(fontsize=fs)

        plt.sca(axs[4])

        im = plt.imshow(rot_v_0, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xlabel("z [Mpc]",fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        #plt.get_yaxis().set_visible(False)
        plt.xticks(fontsize=fs)

        #cb = fig.colorbar(im, ax=ax,fraction=0.03, pad=0.14)#,pad = 0.04 )  # , orientation="horizontal")
        #cb.set_label(r'$\omega_x [km.s^{-1}.kpc^{-1}]$', size='large')

        cb = f.colorbar(im, ax=axs.ravel().tolist(),label=r'$\omega_x~[km~s^{-1}~kpc^{-1}]$', shrink=0.6,
                                pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/rot_v_filament_left.png', dpi=600, format='png')

        print('fig saved')




        # Adjust the layout to prevent overlap
        #plt.tight_layout()

        #plt.subplots_adjust(wspace=0.01, hspace=0.01)

        print("showing plot")

        # Show the plot
        plt.show()

        sys.exit()

    def right_fil_T_plot():
        print("loading maps")

        #t8 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_T_x_d0.02.bin', 'x')
        t6 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_6Mpc_map_T_x_d0.02.bin', 'x')
        t4 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_4Mpc_map_T_x_d0.02.bin', 'x')
        t2 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_2Mpc_map_T_x_d0.02.bin', 'x')
        t0 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_0Mpc_map_T_x_d0.02.bin', 'x')

        print("maps loaded")

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        f, axs = plt.subplots(1, 4, figsize=(16, 4), constrained_layout=True, facecolor='white')
        # gs = gridspec.GridSpec(4, 4, wspace=0.001, hspace=0)

        fs = 14

        title = "$\Delta x_{cen}$= -8Mpc"

        # First row : temperature

        # axs[0, 0].plot(x, np.sin(x))
        # ax = fig.add_subplot(gs[0, 0])
        # plt.sca(axs[0, 0])
        plt.sca(axs[0])
        title = "$\Delta x_{cen}$= 0Mpc"
        plt.title(title, fontsize=fs)
        t0 /= (kb / 1.602e-16)
        plt.imshow(t0, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        plt.xticks([])
        plt.ylabel("y [Mpc]", fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 1])
        # plt.sca(axs[0, 1])
        plt.sca(axs[1])
        title = "$\Delta x_{cen}$= 2Mpc"
        plt.title(title, fontsize=fs)
        t2 /= (kb / 1.602e-16)
        plt.imshow(t2, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 2])
        # plt.sca(axs[0, 2])
        plt.sca(axs[2])
        title = "$\Delta x_{cen}$= 4Mpc"
        plt.title(title, fontsize=fs)
        t4 /= (kb / 1.602e-16)
        plt.imshow(t4, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 3])
        plt.sca(axs[3])
        title = "$\Delta x_{cen}$= 6Mpc"
        plt.title(title, fontsize=fs)
        t6 /= (kb / 1.602e-16)
        im = plt.imshow(t6, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.14)  # , orientation="horizontal")
        # cb.set_label('T[K]', size='large')
        #plt.sca(axs[4])
        #title = "$\Delta x_{cen}$= 0Mpc"
        #plt.title(title, fontsize=fs)
        #t0 /= (kb / 1.602e-16)
        #im = plt.imshow(t0, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        #plt.xticks([])
        #plt.yticks([])
        #plt.grid(b=None)

        # ratiov = t0.shape[0] / t0.shape[1]

        # print("ratiov",ratiov)

        # sys.exit()

        cb = f.colorbar(im, ax=axs.ravel().tolist(), label='T[k]', shrink=0.7, pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)
        cb.set_label('T[k]', size=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/T_filament_right.png', dpi=600, format='png')

        print("plot saved")

        plt.show()

        sys.exit()

    def right_fil_v_plot():

        #vz8 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_vz_x_d0.02.bin', 'x')
        vz6 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_6Mpc_map_vz_x_d0.02.bin', 'x')
        vz4 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_4Mpc_map_vz_x_d0.02.bin', 'x')
        vz2 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_2Mpc_map_vz_x_d0.02.bin', 'x')
        vz0 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_0Mpc_map_vz_x_d0.02.bin', 'x')

        print("vz maps loaded")

        vz_virgo = -131.9249

        #vz8 -= vz_virgo
        vz6 -= vz_virgo
        vz4 -= vz_virgo
        vz2 -= vz_virgo
        vz0 -= vz_virgo

        #mean10_map_8 = np.load("stream_vz_left_fil_minus_8Mpc.npy")
        #mean10_map2_8 = np.load("stream_vy_left_fil_minus_8Mpc.npy")
        mean10_map_6 = np.load("stream_vz_right_fil_plus_6Mpc.npy")
        mean10_map2_6 = np.load("stream_vy_right_fil_plus_6Mpc.npy")
        mean10_map_4 = np.load("stream_vz_right_fil_plus_4Mpc.npy")
        mean10_map2_4 = np.load("stream_vy_right_fil_plus_4Mpc.npy")
        mean10_map_2 = np.load("stream_vz_right_fil_plus_2Mpc.npy")
        mean10_map2_2 = np.load("stream_vy_right_fil_plus_2Mpc.npy")
        mean10_map_0 = np.load("stream_vz_left_fil_minus_0Mpc.npy")
        mean10_map0_2 = np.load("stream_vy_left_fil_minus_0Mpc.npy")

        x = np.load("mesh_x.npy")
        y = np.load("mesh_y.npy")

        print("streamlines loaded")

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        # Create a figure and a 4x4 grid of subplots
        # fig, axs = plt.subplots(4, 4, figsize=(16, 16), facecolor='white')
        # fig = plt.figure(figsize=(16, 16), facecolor='white', constrained_layout=True)

        # Second row : vz + streams

        f, axs = plt.subplots(1, 4, figsize=(16, 4), constrained_layout=True, facecolor='white')

        fs = 14

        #map_mean_mesh_8 = np.sqrt(mean10_map_8 ** 2 + mean10_map2_8 ** 2)
        #map_mean_mesh_8[np.isnan(map_mean_mesh_8) == True] = 0

        #norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 0])
        #plt.sca(axs[0])

        #stream = plt.streamplot(x, y, mean10_map_8, mean10_map2_8, linewidth=0.3, color=map_mean_mesh_8, cmap='Greys_r',
        #                        density=[1, 1.5], norm=norm)
        #im = plt.imshow(vz8, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
        #                vmax=1600)  # BR #origin='lower'
        #plt.xticks([])
        #plt.ylabel("y [Mpc]", fontsize=fs)
        #plt.yticks(fontsize=fs)
        #plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)

        map_mean_mesh_6 = np.sqrt(mean10_map_6 ** 2 + mean10_map2_6 ** 2)
        map_mean_mesh_6[np.isnan(map_mean_mesh_6) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 1])

        plt.sca(axs[3])

        stream = plt.streamplot(x, y, mean10_map_6, mean10_map2_6, linewidth=0.4, color=map_mean_mesh_6, cmap='Greys_r',
                                density=[1, 1.5], norm=norm)
        im = plt.imshow(vz6, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        map_mean_mesh_4 = np.sqrt(mean10_map_4 ** 2 + mean10_map2_4 ** 2)
        map_mean_mesh_4[np.isnan(map_mean_mesh_4) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 2])
        plt.sca(axs[2])

        stream = plt.streamplot(x, y, mean10_map_4, mean10_map2_4, linewidth=0.4, color=map_mean_mesh_4, cmap='Greys_r',
                                density=[1, 1.5], norm=norm)
        im = plt.imshow(vz4, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        map_mean_mesh_2 = np.sqrt(mean10_map_2 ** 2 + mean10_map2_2 ** 2)
        map_mean_mesh_2[np.isnan(map_mean_mesh_2) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 3])
        plt.sca(axs[1])

        stream = plt.streamplot(x, y, mean10_map_2, mean10_map2_2, linewidth=0.4, color=map_mean_mesh_2, cmap='Greys_r',
                                density=[1, 1.5], norm=norm)
        im = plt.imshow(vz2, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)

        plt.sca(axs[0])

        map_mean_mesh_0 = np.sqrt(mean10_map_0 ** 2 + mean10_map0_2 ** 2)
        map_mean_mesh_0[np.isnan(map_mean_mesh_0) == True] = 0

        stream = plt.streamplot(x, y, mean10_map_0, mean10_map0_2, linewidth=0.4, color=map_mean_mesh_0, cmap='Greys_r',
                                density=[1, 1.5], norm=norm)
        im = plt.imshow(vz0, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks([])
        #plt.yticks([])
        plt.grid(b=None)
        plt.ylabel("y [Mpc]", fontsize=fs)
        plt.yticks(fontsize=fs)
        # ax.get_xaxis().set_visible(False)
        # ax.get_yaxis().set_visible(False)

        # cb = fig.colorbar(im, ax=ax,pad = 0.14, fraction=0.03 )  # , orientation="horizontal")
        # cb.set_label('vz[km/s]', size='large')

        # cbs = fig.colorbar(stream.lines, ax=ax, orientation='vertical', pad=0.24, fraction=0.03)
        # cbs.set_label(r'$||\vec{v}||~[km/s]$', size='large')

        cb = f.colorbar(im, ax=axs.ravel().tolist(), label='$v_z~\mathrm{[km~s^{-1}]}$', shrink=0.7, pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)

        #cb = f.colorbar(stream.lines, ax=axs.ravel().tolist(), label=r'$||\vec{v}||~[km/s]$', shrink=0.7, pad=0.01)  # ,size='small')#,shrink=0.6)
        #cb.ax.tick_params(labelsize=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/v_filament_right.png', dpi=600, format='png')

        print("plot saved")

        plt.show()
        sys.exit()
        
    def right_fil_div_plot():

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        #div_v_8 = np.load("./filament/left_fil_transverse_vel_ops/div_v_minus_8Mpc.npy")
        div_v_6 = np.load("./filament/right_fil_transverse_vel_ops/div_v_plus_6Mpc.npy")
        div_v_4 = np.load("./filament/right_fil_transverse_vel_ops/div_v_plus_4Mpc.npy")
        div_v_2 = np.load("./filament/right_fil_transverse_vel_ops/div_v_plus_2Mpc.npy")
        div_v_0 = np.load("./filament/left_fil_transverse_vel_ops/div_v_minus_0Mpc.npy")

        print("div_v maps loaded")

        #div_v_8 = div_v_8.T
        #div_v_8 = np.flip(div_v_8, axis=0)
        div_v_6 = div_v_6.T
        div_v_6 = np.flip(div_v_6, axis=0)
        div_v_4 = div_v_4.T
        div_v_4 = np.flip(div_v_4, axis=0)
        div_v_2 = div_v_2.T
        div_v_2 = np.flip(div_v_2, axis=0)
        div_v_0 = div_v_0.T
        div_v_0 = np.flip(div_v_0, axis=0)


        #Third row : div_v

        #ax = fig.add_subplot(gs[2, 0])

        f, axs = plt.subplots(1, 4, figsize=(16, 4), constrained_layout=True, facecolor='white')

        fs = 14

        plt.sca(axs[0])

        im = plt.imshow(div_v_0, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks([])
        plt.ylabel("y [Mpc]",fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(b=None)
        #plt.get_xaxis().set_visible(False)

        #ax = fig.add_subplot(gs[2, 1])

        plt.sca(axs[1])

        im = plt.imshow(div_v_2, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        #plt.get_xaxis().set_visible(False)
        #plt.get_yaxis().set_visible(False)

        #ax = fig.add_subplot(gs[2, 2])

        plt.sca(axs[2])

        im = plt.imshow(div_v_4, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        #ax.get_xaxis().set_visible(False)
        #ax.get_yaxis().set_visible(False)

        #ax = fig.add_subplot(gs[2, 3])

        plt.sca(axs[3])

        im = plt.imshow(div_v_6, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks([])
        plt.yticks([])

        plt.grid(b=None)
        #plt.get_xaxis().set_visible(False)
        #plt.get_yaxis().set_visible(False)


        #cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.14)#,pad = 0.04 )  # , orientation="horizontal")
        #cb.set_label(r'$\vec{\nabla}.\vec{v}[km.s^{-1}.kpc^{-1}]$', size='large')

        #plt.sca(axs[4])

        #im = plt.imshow(div_v_0, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        #plt.xticks([])
        #plt.yticks([])

        #plt.grid(b=None)

        cb = f.colorbar(im, ax=axs.ravel().tolist(), label=r'$\vec{\nabla}.\vec{v}~[km~s^{-1}~kpc^{-1}]$', shrink=0.7,
                            pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/div_v_filament_right.png', dpi=600, format='png')

        print('fig saved')

        #Fourth row : rot_v

        #ax = fig.add_subplot(gs[3, 0])

        plt.show()
        sys.exit()

    def right_fil_rot_v_plot():

        #rot_v_8 = np.load("./filament/right_fil_transverse_vel_ops/rot_vx_minus_8Mpc.npy")
        rot_v_6 = np.load("./filament/right_fil_transverse_vel_ops/rot_vx_plus_6Mpc.npy")
        rot_v_4 = np.load("./filament/right_fil_transverse_vel_ops/rot_vx_plus_4Mpc.npy")
        rot_v_2 = np.load("./filament/right_fil_transverse_vel_ops/rot_vx_plus_2Mpc.npy")
        rot_v_0 = np.load("./filament/left_fil_transverse_vel_ops/rot_vx_minus_0Mpc.npy")

        print("rot_v maps loaded")

        #rot_v_8 = rot_v_8.T
        #rot_v_8 = np.flip(rot_v_8, axis=0)
        rot_v_6 = rot_v_6.T
        rot_v_6 = np.flip(rot_v_6, axis=0)
        rot_v_4 = rot_v_4.T
        rot_v_4 = np.flip(rot_v_4, axis=0)
        rot_v_2 = rot_v_2.T
        rot_v_2 = np.flip(rot_v_2, axis=0)
        rot_v_0 = rot_v_0.T
        rot_v_0 = np.flip(rot_v_0, axis=0)

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        fs = 14

        f, axs = plt.subplots(1, 4, figsize=(16, 4), constrained_layout=True, facecolor='white')

        plt.sca(axs[0])

        im = plt.imshow(rot_v_0, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xlabel("z [Mpc]", fontsize=fs)
        plt.ylabel("y [Mpc]", fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.grid(b=None)

        # ax = fig.add_subplot(gs[3, 1])

        plt.sca(axs[1])

        im = plt.imshow(rot_v_2, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xlabel("z [Mpc]", fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_yaxis().set_visible(False)
        plt.xticks(fontsize=fs)

        # ax = fig.add_subplot(gs[3, 2])

        plt.sca(axs[2])

        im = plt.imshow(rot_v_4, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xlabel("z [Mpc]", fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_yaxis().set_visible(False)
        plt.xticks(fontsize=fs)

        # ax = fig.add_subplot(gs[3, 3])

        plt.sca(axs[3])

        im = plt.imshow(rot_v_6, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xlabel("z [Mpc]", fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_yaxis().set_visible(False)
        plt.xticks(fontsize=fs)

        #plt.sca(axs[4])

        #im = plt.imshow(rot_v_0, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        #plt.xlabel("z [Mpc]", fontsize=fs)
        #plt.yticks([])
        #plt.grid(b=None)
        # plt.get_yaxis().set_visible(False)
        #plt.xticks(fontsize=fs)

        # cb = fig.colorbar(im, ax=ax,fraction=0.03, pad=0.14)#,pad = 0.04 )  # , orientation="horizontal")
        # cb.set_label(r'$\omega_x [km.s^{-1}.kpc^{-1}]$', size='large')

        cb = f.colorbar(im, ax=axs.ravel().tolist(), label=r'$\omega_x~[km~s^{-1}~kpc^{-1}]$', shrink=0.7,
                        pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/rot_v_filament_right.png', dpi=600, format='png')

        print('fig saved')

        # Adjust the layout to prevent overlap
        # plt.tight_layout()

        # plt.subplots_adjust(wspace=0.01, hspace=0.01)

        print("showing plot")

        # Show the plot
        plt.show()

        sys.exit()

    def two_fil_v_plot():

        h = FortranFile('./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_vz_x_d0.02.bin', 'r')

        nx, ny, nz = h.read_ints()

        vzm8 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_vz_x_d0.02.bin', 'x')
        vzm6 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_6Mpc_map_vz_x_d0.02.bin', 'x')
        vzm4 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_4Mpc_map_vz_x_d0.02.bin', 'x')
        vzm2 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_2Mpc_map_vz_x_d0.02.bin', 'x')
        vz0 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_0Mpc_map_vz_x_d0.02.bin', 'x')

        vz6 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_6Mpc_map_vz_x_d0.02.bin', 'x')
        vz4 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_4Mpc_map_vz_x_d0.02.bin', 'x')
        vz2 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_2Mpc_map_vz_x_d0.02.bin', 'x')


        print("vz maps loaded")


        vz_virgo = -131.9249

        vzm8 -= vz_virgo
        vzm6 -= vz_virgo
        vzm4 -= vz_virgo
        vzm2 -= vz_virgo
        vz0 -= vz_virgo

        vz6 -= vz_virgo
        vz4 -= vz_virgo
        vz2 -= vz_virgo


        mean10_map_m8 = np.load("stream_vz_left_fil_minus_8Mpc.npy")
        mean10_map2_m8 = np.load("stream_vy_left_fil_minus_8Mpc.npy")
        mean10_map_m6 = np.load("stream_vz_left_fil_minus_6Mpc.npy")
        mean10_map2_m6 = np.load("stream_vy_left_fil_minus_6Mpc.npy")
        mean10_map_m4 = np.load("stream_vz_left_fil_minus_4Mpc.npy")
        mean10_map2_m4 = np.load("stream_vy_left_fil_minus_4Mpc.npy")
        mean10_map_m2 = np.load("stream_vz_left_fil_minus_2Mpc.npy")
        mean10_map2_m2 = np.load("stream_vy_left_fil_minus_2Mpc.npy")
        mean10_map_0 = np.load("stream_vz_left_fil_minus_0Mpc.npy")
        mean10_map2_0 = np.load("stream_vy_left_fil_minus_0Mpc.npy")

        mean10_map_6 = np.load("stream_vz_right_fil_plus_6Mpc.npy")
        mean10_map2_6 = np.load("stream_vy_right_fil_plus_6Mpc.npy")
        mean10_map_4 = np.load("stream_vz_right_fil_plus_4Mpc.npy")
        mean10_map2_4 = np.load("stream_vy_right_fil_plus_4Mpc.npy")
        mean10_map_2 = np.load("stream_vz_right_fil_plus_2Mpc.npy")
        mean10_map2_2 = np.load("stream_vy_right_fil_plus_2Mpc.npy")

        x = np.load("mesh_x.npy")
        y = np.load("mesh_y.npy")

        print("streamlines loaded")

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        # Create a figure and a 4x4 grid of subplots
        # fig, axs = plt.subplots(4, 4, figsize=(16, 16), facecolor='white')
        # fig = plt.figure(figsize=(16, 16), facecolor='white', constrained_layout=True)

        # Second row : vz + streams

        f, axs = plt.subplots(1, 8, figsize=(16, 4), constrained_layout=True, facecolor='white')

        fs = 14

        size_px = 737.441/2**19

        cf_1Mpc = 1/size_px

        print("size_px",size_px)

        print("cf_1Mpc",cf_1Mpc)

        #sys.exit()

        map_mean_mesh_m8 = np.sqrt(mean10_map_m8 ** 2 + mean10_map2_m8 ** 2)
        map_mean_mesh_m8[np.isnan(map_mean_mesh_m8) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 0])
        plt.sca(axs[0])

        stream = plt.streamplot(x, y, mean10_map_m8, mean10_map2_m8, linewidth=0.3, color=map_mean_mesh_m8, cmap='Greys_r', density=[0.5, 0.5], norm=norm)
        im = plt.imshow(vzm8, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,vmax=1600)  # BR #origin='lower'

        anchor_z_5 = -1
        anchor_y_5 = 1
        length_5 = 5
        plt.gca().add_patch(Rectangle((anchor_z_5,anchor_y_5),length_5,length_5,edgecolor='black',facecolor='none',lw=1,ls="solid"))

        anchor_z_2 = 0.5
        anchor_y_2 = 2.5
        length_2 = 2
        plt.gca().add_patch(Rectangle((anchor_z_2,anchor_y_2),length_2,length_2,edgecolor='black',facecolor='none',lw=1,ls="dashed"))

        anchor_z_1 = 1
        anchor_y_1 = 3
        length_1 = 1
        plt.gca().add_patch(Rectangle((anchor_z_1,anchor_y_1),length_1,length_1,edgecolor='black',facecolor='none',lw=1,ls="dotted"))

        plt.axvline(x=1.5, color='black', linestyle='--', lw=1)
        plt.axhline(y=3.5, color='black', linestyle='--', lw=1)

        plt.xlabel("z [Mpc]", fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.ylabel("y [Mpc]", fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        plt.title("$\Delta x_{cen}$= -8Mpc", fontsize=fs)


        map_mean_mesh_m6 = np.sqrt(mean10_map_m6 ** 2 + mean10_map2_m6 ** 2)
        map_mean_mesh_m6[np.isnan(map_mean_mesh_m6) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 1])

        plt.sca(axs[1])

        stream = plt.streamplot(x, y, mean10_map_m6, mean10_map2_m6, linewidth=0.2, color=map_mean_mesh_m6, cmap='Greys_r',density=[0.5, 0.5], norm=norm)
        im = plt.imshow(vzm6, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'

        anchor_z_5 = -1.5
        anchor_y_5 = -1
        length_5 = 5
        plt.gca().add_patch(Rectangle((anchor_z_5,anchor_y_5),length_5,length_5,edgecolor='black',facecolor='none',lw=1,ls="solid"))

        anchor_z_2 = 0
        anchor_y_2 = 0.5
        length_2 = 2
        plt.gca().add_patch(Rectangle((anchor_z_2,anchor_y_2),length_2,length_2,edgecolor='black',facecolor='none',lw=1,ls="dashed"))

        anchor_z_1 = 0.5
        anchor_y_1 = 1
        length_1 = 1
        plt.gca().add_patch(Rectangle((anchor_z_1,anchor_y_1),length_1,length_1,edgecolor='black',facecolor='none',lw=1,ls="dotted"))

        plt.axvline(x=1, color='black', linestyle='--', lw=1)
        plt.axhline(y=1.5, color='black', linestyle='--', lw=1)

        plt.xlabel("z [Mpc]", fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        plt.title("$\Delta x_{cen}$= -6Mpc", fontsize=fs)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        map_mean_mesh_m4 = np.sqrt(mean10_map_m4 ** 2 + mean10_map2_m4 ** 2)
        map_mean_mesh_m4[np.isnan(map_mean_mesh_m4) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 2])
        plt.sca(axs[2])

        stream = plt.streamplot(x, y, mean10_map_m4, mean10_map2_m4, linewidth=0.2, color=map_mean_mesh_m4, cmap='Greys_r',density=[0.5, 0.5], norm=norm)
        im = plt.imshow(vzm4, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'

        anchor_z_5 = -2
        anchor_y_5 = -1
        length_5 = 5
        plt.gca().add_patch(
            Rectangle((anchor_z_5, anchor_y_5), length_5, length_5, edgecolor='black', facecolor='none', lw=1,
                      ls="solid"))

        anchor_z_2 = -0.5
        anchor_y_2 = 0.5
        length_2 = 2
        plt.gca().add_patch(
            Rectangle((anchor_z_2, anchor_y_2), length_2, length_2, edgecolor='black', facecolor='none', lw=1,
                      ls="dashed"))

        anchor_z_1 = 0
        anchor_y_1 = 1
        length_1 = 1
        plt.gca().add_patch(
            Rectangle((anchor_z_1, anchor_y_1), length_1, length_1, edgecolor='black', facecolor='none', lw=1,
                      ls="dotted"))

        plt.axvline(x=0.5, color='black', linestyle='--', lw=1)
        plt.axhline(y=1.5, color='black', linestyle='--', lw=1)

        plt.xlabel("z [Mpc]", fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        plt.title("$\Delta x_{cen}$= -4Mpc", fontsize=fs)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        map_mean_mesh_m2 = np.sqrt(mean10_map_m2 ** 2 + mean10_map2_m2 ** 2)
        map_mean_mesh_m2[np.isnan(map_mean_mesh_m2) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 3])
        plt.sca(axs[3])

        stream = plt.streamplot(x, y, mean10_map_m2, mean10_map2_m2, linewidth=0.2, color=map_mean_mesh_m2, cmap='Greys_r',density=[0.5, 0.5], norm=norm)
        im = plt.imshow(vzm2, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'

        anchor_z_5 = -2.5
        anchor_y_5 = -2.5
        length_5 = 5
        plt.gca().add_patch(
            Rectangle((anchor_z_5, anchor_y_5), length_5, length_5, edgecolor='black', facecolor='none', lw=1,
                      ls="solid"))

        anchor_z_2 = -1
        anchor_y_2 = -1
        length_2 = 2
        plt.gca().add_patch(
            Rectangle((anchor_z_2, anchor_y_2), length_2, length_2, edgecolor='black', facecolor='none', lw=1,
                      ls="dashed"))

        anchor_z_1 = -0.5
        anchor_y_1 = -0.5
        length_1 = 1
        plt.gca().add_patch(
            Rectangle((anchor_z_1, anchor_y_1), length_1, length_1, edgecolor='black', facecolor='none', lw=1,
                      ls="dotted"))

        plt.axvline(x=0, color='black', linestyle='--', lw=1)
        plt.axhline(y=0, color='black', linestyle='--', lw=1)

        plt.xlabel("z [Mpc]", fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        plt.title("$\Delta x_{cen}$= -2Mpc", fontsize=fs)

        map_mean_mesh_0 = np.sqrt(mean10_map_0 ** 2 + mean10_map2_0 ** 2)
        map_mean_mesh_0[np.isnan(map_mean_mesh_0) == True] = 0

        plt.sca(axs[4])

        stream = plt.streamplot(x, y, mean10_map_0, mean10_map2_0, linewidth=0.2, color=map_mean_mesh_0, cmap='Greys_r',density=[0.5, 0.5], norm=norm)
        im = plt.imshow(vz0, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'

        anchor_z_5 = -2.5
        anchor_y_5 = -2.5
        length_5 = 5
        plt.gca().add_patch(
            Rectangle((anchor_z_5, anchor_y_5), length_5, length_5, edgecolor='black', facecolor='none', lw=1,
                      ls="solid"))

        anchor_z_2 = -1
        anchor_y_2 = -1
        length_2 = 2
        plt.gca().add_patch(
            Rectangle((anchor_z_2, anchor_y_2), length_2, length_2, edgecolor='black', facecolor='none', lw=1,
                      ls="dashed"))

        anchor_z_1 = -0.5
        anchor_y_1 = -0.5
        length_1 = 1
        plt.gca().add_patch(
            Rectangle((anchor_z_1, anchor_y_1), length_1, length_1, edgecolor='black', facecolor='none', lw=1,
                      ls="dotted"))

        plt.axvline(x=0, color='black', linestyle='--', lw=1)
        plt.axhline(y=0, color='black', linestyle='--', lw=1)

        plt.xlabel("z [Mpc]", fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        plt.title("$\Delta x_{cen}$= 0Mpc", fontsize=fs)
        # ax.get_xaxis().set_visible(False)
        # ax.get_yaxis().set_visible(False)

        # cb = fig.colorbar(im, ax=ax,pad = 0.14, fraction=0.03 )  # , orientation="horizontal")
        # cb.set_label('vz[km/s]', size='large')

        # cbs = fig.colorbar(stream.lines, ax=ax, orientation='vertical', pad=0.24, fraction=0.03)
        # cbs.set_label(r'$||\vec{v}||~[km/s]$', size='large')

        #cb = f.colorbar(im, ax=axs.ravel().tolist(), label='$v_z~\mathrm{[km~s^{-1}]}$', shrink=0.6, pad=0.01)  # ,size='small')#,shrink=0.6)
        #cb.ax.tick_params(labelsize=fs)

        map_mean_mesh_6 = np.sqrt(mean10_map_6 ** 2 + mean10_map2_6 ** 2)
        map_mean_mesh_6[np.isnan(map_mean_mesh_6) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 1])

        plt.sca(axs[7])

        stream = plt.streamplot(x, y, mean10_map_6, mean10_map2_6, linewidth=0.2, color=map_mean_mesh_6, cmap='Greys_r', density=[0.5, 0.5], norm=norm)
        im = plt.imshow(vz6, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'

        anchor_z_5 = -3
        anchor_y_5 = -5
        length_5 = 5
        plt.gca().add_patch(
            Rectangle((anchor_z_5, anchor_y_5), length_5, length_5, edgecolor='black', facecolor='none', lw=1,
                      ls="solid"))

        anchor_z_2 = -1.5
        anchor_y_2 = -3.5
        length_2 = 2
        plt.gca().add_patch(
            Rectangle((anchor_z_2, anchor_y_2), length_2, length_2, edgecolor='black', facecolor='none', lw=1,
                      ls="dashed"))

        anchor_z_1 = -1
        anchor_y_1 = -3
        length_1 = 1
        plt.gca().add_patch(
            Rectangle((anchor_z_1, anchor_y_1), length_1, length_1, edgecolor='black', facecolor='none', lw=1,
                      ls="dotted"))

        plt.axvline(x=-0.5, color='black', linestyle='--', lw=1)
        plt.axhline(y=-2.5, color='black', linestyle='--', lw=1)

        plt.xlabel("z [Mpc]", fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        plt.title("$\Delta x_{cen}$= 6Mpc", fontsize=fs)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        map_mean_mesh_4 = np.sqrt(mean10_map_4 ** 2 + mean10_map2_4 ** 2)
        map_mean_mesh_4[np.isnan(map_mean_mesh_4) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 2])
        plt.sca(axs[6])

        stream = plt.streamplot(x, y, mean10_map_4, mean10_map2_4, linewidth=0.2, color=map_mean_mesh_4, cmap='Greys_r', density=[0.5, 0.5], norm=norm)
        im = plt.imshow(vz4, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'

        anchor_z_5 = -2.5
        anchor_y_5 = -3.5
        length_5 = 5
        plt.gca().add_patch(
            Rectangle((anchor_z_5, anchor_y_5), length_5, length_5, edgecolor='black', facecolor='none', lw=1,
                      ls="solid"))

        anchor_z_2 = -1
        anchor_y_2 = -2
        length_2 = 2
        plt.gca().add_patch(
            Rectangle((anchor_z_2, anchor_y_2), length_2, length_2, edgecolor='black', facecolor='none', lw=1,
                      ls="dashed"))

        anchor_z_1 = -0.5
        anchor_y_1 = -1.5
        length_1 = 1
        plt.gca().add_patch(
            Rectangle((anchor_z_1, anchor_y_1), length_1, length_1, edgecolor='black', facecolor='none', lw=1,
                      ls="dotted"))

        plt.axvline(x=0, color='black', linestyle='--', lw=1)
        plt.axhline(y=-1, color='black', linestyle='--', lw=1)

        plt.xlabel("z [Mpc]", fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        plt.title("$\Delta x_{cen}$= 4Mpc", fontsize=fs)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        map_mean_mesh_2 = np.sqrt(mean10_map_2 ** 2 + mean10_map2_2 ** 2)
        map_mean_mesh_2[np.isnan(map_mean_mesh_2) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 3])
        plt.sca(axs[5])

        stream = plt.streamplot(x, y, mean10_map_2, mean10_map2_2, linewidth=0.2, color=map_mean_mesh_2, cmap='Greys_r',density=[0.5, 0.5], norm=norm)
        im = plt.imshow(vz2, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'

        anchor_z_5 = -2.5
        anchor_y_5 = -2.5
        length_5 = 5
        plt.gca().add_patch(
            Rectangle((anchor_z_5, anchor_y_5), length_5, length_5, edgecolor='black', facecolor='none', lw=1,
                      ls="solid"))

        anchor_z_2 = -1
        anchor_y_2 = -1
        length_2 = 2
        plt.gca().add_patch(
            Rectangle((anchor_z_2, anchor_y_2), length_2, length_2, edgecolor='black', facecolor='none', lw=1,
                      ls="dashed"))

        anchor_z_1 = -0.5
        anchor_y_1 = -0.5
        length_1 = 1
        plt.gca().add_patch(
            Rectangle((anchor_z_1, anchor_y_1), length_1, length_1, edgecolor='black', facecolor='none', lw=1,
                      ls="dotted"))

        plt.axvline(x=0, color='black', linestyle='--', lw=1)
        plt.axhline(y=0, color='black', linestyle='--', lw=1)

        plt.xlabel("z [Mpc]", fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        plt.title("$\Delta x_{cen}$= 2Mpc", fontsize=fs)

        ratiov = (vz2.shape[0] / vz2.shape[1])/8

        print("ratiov",ratiov,"shape",vz2.shape[0],vz2.shape[1])

        #sys.exit()

        cb = f.colorbar(im, ax=axs.ravel().tolist(), label=r'$v_z~[km/s]$', fraction=0.047 * ratiov, pad=0.01)
        cb.ax.tick_params(labelsize=fs)

        #plt.savefig('./papers_plots/paper_3_turb_fil/v_filament_left_2.png', dpi=600, format='png')

        #print("plot saved")

        plt.show()
        sys.exit()

    def left_fil_T_plot_2():
        print("loading maps")

        t8 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_T_x_d0.02.bin', 'x')
        t6 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_6Mpc_map_T_x_d0.02.bin', 'x')
        t4 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_4Mpc_map_T_x_d0.02.bin', 'x')
        t2 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_2Mpc_map_T_x_d0.02.bin', 'x')
        t0 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_0Mpc_map_T_x_d0.02.bin', 'x')

        print("maps loaded")

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        f, axs = plt.subplots(1, 4, figsize=(16, 8), constrained_layout=True, facecolor='white')
        # gs = gridspec.GridSpec(4, 4, wspace=0.001, hspace=0)

        fs = 14

        #title = "$\Delta x_{cen}$= -8Mpc"

        # First row : temperature

        # axs[0, 0].plot(x, np.sin(x))
        # ax = fig.add_subplot(gs[0, 0])
        # plt.sca(axs[0, 0])
        plt.sca(axs[0])
        title = "$\Delta x_{cen}$= -8Mpc"
        plt.title(title, fontsize=fs)
        t8 /= (kb / 1.602e-16)
        plt.imshow(t8, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        plt.xticks([])
        plt.ylabel("y [Mpc]", fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 1])
        # plt.sca(axs[0, 1])
        plt.sca(axs[1])
        title = "$\Delta x_{cen}$= -6Mpc"
        plt.title(title, fontsize=fs)
        t6 /= (kb / 1.602e-16)
        plt.imshow(t6, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        plt.xticks([])
        plt.yticks([])
        #plt.ylabel("y [Mpc]", fontsize=fs)
        #plt.yticks(fontsize=fs)
        #plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 2])
        # plt.sca(axs[0, 2])
        plt.sca(axs[2])
        title = "$\Delta x_{cen}$= -4Mpc"
        plt.title(title, fontsize=fs)
        t4 /= (kb / 1.602e-16)
        plt.imshow(t4, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 3])
        plt.sca(axs[3])
        title = "$\Delta x_{cen}$= -2Mpc"
        plt.title(title, fontsize=fs)
        t2 /= (kb / 1.602e-16)
        im = plt.imshow(t2, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.14)  # , orientation="horizontal")
        # cb.set_label('T[K]', size='large')
        #plt.sca(axs[4])
        #title = "$\Delta x_{cen}$= 0Mpc"
        #plt.title(title, fontsize=fs)
        #t0 /= (kb / 1.602e-16)
        #im = plt.imshow(t0, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        #plt.xticks([])
        #plt.yticks([])
        #plt.grid(b=None)

        # ratiov = t0.shape[0] / t0.shape[1]

        # print("ratiov",ratiov)

        # sys.exit()

        cb = f.colorbar(im, ax=axs.ravel().tolist(), shrink=0.5, pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)
        cb.set_label('$T~\mathrm{[K]}$', size=fs+2)

        plt.savefig('./papers_plots/paper_3_turb_fil/T_filament_left_2.png', dpi=600, format='png')

        print("plot saved")

        plt.show()

        sys.exit()

    def left_fil_ne_plot_2():
        print("loading maps")

        ne8 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_ne_x_d0.02.bin', 'x')
        ne6 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_6Mpc_map_ne_x_d0.02.bin', 'x')
        ne4 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_4Mpc_map_ne_x_d0.02.bin', 'x')
        ne2 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_2Mpc_map_ne_x_d0.02.bin', 'x')
        #t0 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_0Mpc_map_ne_x_d0.02.bin', 'x')

        print("maps loaded")

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        f, axs = plt.subplots(1, 4, figsize=(16, 8), constrained_layout=True, facecolor='white')
        # gs = gridspec.GridSpec(4, 4, wspace=0.001, hspace=0)

        fs = 14

        #title = "$\Delta x_{cen}$= -8Mpc"

        # First row : temperature

        # axs[0, 0].plot(x, np.sin(x))
        # ax = fig.add_subplot(gs[0, 0])
        # plt.sca(axs[0, 0])

        plt.sca(axs[0])
        title = "$\Delta x_{cen}$= -8Mpc"
        plt.title(title, fontsize=fs)
        #t8 /= (kb / 1.602e-16)
        plt.imshow(ne8, cmap="Solar", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=10**(-7.5), vmax=10**(-3)))
        #plt.imshow(ne8, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        #plt.xticks([])
        plt.xlabel("x [Mpc]", fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.ylabel("y [Mpc]", fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 1])
        # plt.sca(axs[0, 1])


        plt.sca(axs[1])
        title = "$\Delta x_{cen}$= -6Mpc"
        plt.title(title, fontsize=fs)
        #t6 /= (kb / 1.602e-16)
        plt.imshow(ne6, cmap="Solar", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=10**(-7.5), vmax=10**(-3)))# norm=colors.LogNorm(vmin=10 ** (-7.5), vmax=10 ** (-4)))
        plt.xlabel("x [Mpc]", fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks([])
        #plt.ylabel("y [Mpc]", fontsize=fs)
        #plt.yticks(fontsize=fs)
        #plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 2])
        # plt.sca(axs[0, 2])
        plt.sca(axs[2])
        title = "$\Delta x_{cen}$= -4Mpc"
        plt.title(title, fontsize=fs)
        #t4 /= (kb / 1.602e-16)
        plt.imshow(ne4, cmap="Solar", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=10**(-7.5), vmax=10**(-3)))#, norm=colors.LogNorm(vmin=10 ** (-7.5), vmax=10 ** (-4)))
        #plt.xticks([])
        plt.xlabel("x [Mpc]", fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 3])
        plt.sca(axs[3])
        title = "$\Delta x_{cen}$= -2Mpc"
        plt.title(title, fontsize=fs)
        #t2 /= (kb / 1.602e-16)
        im = plt.imshow(ne2, cmap="Solar", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=10**(-7.5), vmax=10**(-3)))#, norm=colors.LogNorm(vmin=10 ** (-7.5), vmax=10 ** (-4)))
        #plt.xticks([])
        plt.xlabel("x [Mpc]", fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.14)  # , orientation="horizontal")
        # cb.set_label('T[K]', size='large')
        #plt.sca(axs[4])
        #title = "$\Delta x_{cen}$= 0Mpc"
        #plt.title(title, fontsize=fs)
        #t0 /= (kb / 1.602e-16)
        #im = plt.imshow(t0, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        #plt.xticks([])
        #plt.yticks([])
        #plt.grid(b=None)

        # ratiov = t0.shape[0] / t0.shape[1]

        # print("ratiov",ratiov)

        # sys.exit()

        cb = f.colorbar(im, ax=axs.ravel().tolist(), shrink=0.5, pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)
        cb.set_label('$n_e~\mathrm{[cm^{-3}]}$', size=fs+2)

        plt.savefig('./papers_plots/paper_3_turb_fil/ne_filament_left_2.png', dpi=600, format='png')

        print("plot saved")

        plt.show()

        sys.exit()

    def left_fil_v_plot_2():

        vz8 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_vz_x_d0.02.bin', 'x')
        vz6 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_6Mpc_map_vz_x_d0.02.bin', 'x')
        vz4 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_4Mpc_map_vz_x_d0.02.bin', 'x')
        vz2 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_2Mpc_map_vz_x_d0.02.bin', 'x')
        #vz0 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_0Mpc_map_vz_x_d0.02.bin', 'x')

        print("vz maps loaded")

        vz_virgo = -131.9249

        vz8 -= vz_virgo
        vz6 -= vz_virgo
        vz4 -= vz_virgo
        vz2 -= vz_virgo
        #vz0 -= vz_virgo

        mean10_map_8 = np.load("stream_vz_left_fil_minus_8Mpc.npy")
        mean10_map2_8 = np.load("stream_vy_left_fil_minus_8Mpc.npy")
        mean10_map_6 = np.load("stream_vz_left_fil_minus_6Mpc.npy")
        mean10_map2_6 = np.load("stream_vy_left_fil_minus_6Mpc.npy")
        mean10_map_4 = np.load("stream_vz_left_fil_minus_4Mpc.npy")
        mean10_map2_4 = np.load("stream_vy_left_fil_minus_4Mpc.npy")
        mean10_map_2 = np.load("stream_vz_left_fil_minus_2Mpc.npy")
        mean10_map2_2 = np.load("stream_vy_left_fil_minus_2Mpc.npy")
        #mean10_map_0 = np.load("stream_vz_left_fil_minus_0Mpc.npy")
        #mean10_map0_2 = np.load("stream_vy_left_fil_minus_0Mpc.npy")

        x = np.load("mesh_x.npy")
        y = np.load("mesh_y.npy")

        print("streamlines loaded")

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        # Create a figure and a 4x4 grid of subplots
        # fig, axs = plt.subplots(4, 4, figsize=(16, 16), facecolor='white')
        # fig = plt.figure(figsize=(16, 16), facecolor='white', constrained_layout=True)

        # Second row : vz + streams

        f, axs = plt.subplots(1, 4, figsize=(16, 8), constrained_layout=True, facecolor='white')

        fs = 14

        map_mean_mesh_8 = np.sqrt(mean10_map_8 ** 2 + mean10_map2_8 ** 2)
        map_mean_mesh_8[np.isnan(map_mean_mesh_8) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 0])
        plt.sca(axs[0])

        stream = plt.streamplot(x, y, mean10_map_8, mean10_map2_8, linewidth=0.3, color=map_mean_mesh_8, cmap='Greys_r',
                                density=[2, 3], norm=norm)
        im = plt.imshow(vz8, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks([])
        plt.ylabel("y [Mpc]", fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)

        map_mean_mesh_6 = np.sqrt(mean10_map_6 ** 2 + mean10_map2_6 ** 2)
        map_mean_mesh_6[np.isnan(map_mean_mesh_6) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 1])

        plt.sca(axs[1])

        stream = plt.streamplot(x, y, mean10_map_6, mean10_map2_6, linewidth=0.4, color=map_mean_mesh_6, cmap='Greys_r',
                                density=[2, 3], norm=norm)
        im = plt.imshow(vz6, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks([])
        plt.yticks([])
        #plt.ylabel("y [Mpc]", fontsize=fs)
        #plt.yticks(fontsize=fs)
        #plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        map_mean_mesh_4 = np.sqrt(mean10_map_4 ** 2 + mean10_map2_4 ** 2)
        map_mean_mesh_4[np.isnan(map_mean_mesh_4) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 2])
        plt.sca(axs[2])

        stream = plt.streamplot(x, y, mean10_map_4, mean10_map2_4, linewidth=0.4, color=map_mean_mesh_4, cmap='Greys_r',
                                density=[2, 3], norm=norm)
        im = plt.imshow(vz4, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        map_mean_mesh_2 = np.sqrt(mean10_map_2 ** 2 + mean10_map2_2 ** 2)
        map_mean_mesh_2[np.isnan(map_mean_mesh_2) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 3])
        plt.sca(axs[3])

        stream = plt.streamplot(x, y, mean10_map_2, mean10_map2_2, linewidth=0.4, color=map_mean_mesh_2, cmap='Greys_r',
                                density=[2, 3], norm=norm)
        im = plt.imshow(vz2, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)

        #plt.sca(axs[4])

        #stream = plt.streamplot(x, y, mean10_map_0, mean10_map0_2, linewidth=0.4, color=map_mean_mesh_2, cmap='Greys_r',
        #                        density=[1, 1.5], norm=norm)
        #im = plt.imshow(vz0, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
        #                vmax=1600)  # BR #origin='lower'
        #plt.xticks([])
        #plt.yticks([])
        #plt.grid(b=None)
        # ax.get_xaxis().set_visible(False)
        # ax.get_yaxis().set_visible(False)

        # cb = fig.colorbar(im, ax=ax,pad = 0.14, fraction=0.03 )  # , orientation="horizontal")
        # cb.set_label('vz[km/s]', size='large')

        # cbs = fig.colorbar(stream.lines, ax=ax, orientation='vertical', pad=0.24, fraction=0.03)
        # cbs.set_label(r'$||\vec{v}||~[km/s]$', size='large')

        #cb = f.colorbar(im, ax=axs.ravel().tolist(), label='$v_z~\mathrm{[km~s^{-1}]}$', shrink=0.6, pad=0.01)  # ,size='small')#,shrink=0.6)
        #cb.ax.tick_params(labelsize=fs)

        cb = f.colorbar(stream.lines, ax=axs.ravel().tolist(), label=r'$||\vec{v}||~[km/s]$', shrink=0.5, pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/v_filament_left_2.png', dpi=600, format='png')

        print("plot saved")

        plt.show()
        sys.exit()

    def left_fil_div_plot_2():

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        div_v_8 = np.load("./filament/left_fil_transverse_vel_ops/div_v_minus_8Mpc.npy")
        div_v_6 = np.load("./filament/left_fil_transverse_vel_ops/div_v_minus_6Mpc.npy")
        div_v_4 = np.load("./filament/left_fil_transverse_vel_ops/div_v_minus_4Mpc.npy")
        div_v_2 = np.load("./filament/left_fil_transverse_vel_ops/div_v_minus_2Mpc.npy")
        #div_v_0 = np.load("./filament/left_fil_transverse_vel_ops/div_v_minus_0Mpc.npy")

        print("div_v maps loaded")

        div_v_8 = div_v_8.T
        div_v_8 = np.flip(div_v_8, axis=0)
        div_v_6 = div_v_6.T
        div_v_6 = np.flip(div_v_6, axis=0)
        div_v_4 = div_v_4.T
        div_v_4 = np.flip(div_v_4, axis=0)
        div_v_2 = div_v_2.T
        div_v_2 = np.flip(div_v_2, axis=0)
        #div_v_0 = div_v_0.T
        #div_v_0 = np.flip(div_v_0, axis=0)


        #Third row : div_v

        #ax = fig.add_subplot(gs[2, 0])

        f, axs = plt.subplots(1, 4, figsize=(16, 8), constrained_layout=True, facecolor='white')

        fs = 14

        plt.sca(axs[0])

        im = plt.imshow(div_v_8, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks([])
        plt.ylabel("y [Mpc]",fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(b=None)
        #plt.get_xaxis().set_visible(False)

        #ax = fig.add_subplot(gs[2, 1])

        plt.sca(axs[1])

        im = plt.imshow(div_v_6, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks([])
        plt.yticks([])
        #plt.ylabel("y [Mpc]",fontsize=fs)
        #plt.yticks(fontsize=fs)
        plt.grid(b=None)
        #plt.get_xaxis().set_visible(False)
        #plt.get_yaxis().set_visible(False)

        #ax = fig.add_subplot(gs[2, 2])

        plt.sca(axs[2])

        im = plt.imshow(div_v_4, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        #ax.get_xaxis().set_visible(False)
        #ax.get_yaxis().set_visible(False)

        #ax = fig.add_subplot(gs[2, 3])

        plt.sca(axs[3])

        im = plt.imshow(div_v_2, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks([])
        plt.yticks([])

        plt.grid(b=None)
        #plt.get_xaxis().set_visible(False)
        #plt.get_yaxis().set_visible(False)


        #cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.14)#,pad = 0.04 )  # , orientation="horizontal")
        #cb.set_label(r'$\vec{\nabla}.\vec{v}[km.s^{-1}.kpc^{-1}]$', size='large')

        #plt.sca(axs[4])

        #im = plt.imshow(div_v_0, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        #plt.xticks([])
        #plt.yticks([])

        #plt.grid(b=None)

        cb = f.colorbar(im, ax=axs.ravel().tolist(), label=r'$\vec{\nabla}.\vec{v}~[km~s^{-1}~kpc^{-1}]$', shrink=0.5,
                            pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/div_v_filament_left_2.png', dpi=600, format='png')

        print('fig saved')

        #Fourth row : rot_v

        #ax = fig.add_subplot(gs[3, 0])

        plt.show()
        sys.exit()

    def left_fil_rot_v_plot_2():

        rot_v_8 = np.load("./filament/left_fil_transverse_vel_ops/rot_vx_minus_8Mpc.npy")
        rot_v_6 = np.load("./filament/left_fil_transverse_vel_ops/rot_vx_minus_6Mpc.npy")
        rot_v_4 = np.load("./filament/left_fil_transverse_vel_ops/rot_vx_minus_4Mpc.npy")
        rot_v_2 = np.load("./filament/left_fil_transverse_vel_ops/rot_vx_minus_2Mpc.npy")
        #rot_v_0 = np.load("./filament/left_fil_transverse_vel_ops/rot_vx_minus_0Mpc.npy")

        print("rot_v maps loaded")

        rot_v_8 = rot_v_8.T
        rot_v_8 = np.flip(rot_v_8, axis=0)
        rot_v_6 = rot_v_6.T
        rot_v_6 = np.flip(rot_v_6, axis=0)
        rot_v_4 = rot_v_4.T
        rot_v_4 = np.flip(rot_v_4, axis=0)
        rot_v_2 = rot_v_2.T
        rot_v_2 = np.flip(rot_v_2, axis=0)
        #rot_v_0 = rot_v_0.T
        #rot_v_0 = np.flip(rot_v_0, axis=0)


        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        fs = 14

        f, axs = plt.subplots(1, 4, figsize=(16, 8), constrained_layout=True, facecolor='white')

        plt.sca(axs[0])


        im = plt.imshow(rot_v_8, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xlabel("z [Mpc]",fontsize=fs)
        plt.ylabel("y [Mpc]", fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.grid(b=None)

        #ax = fig.add_subplot(gs[3, 1])

        plt.sca(axs[1])

        im = plt.imshow(rot_v_6, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xlabel("z [Mpc]",fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        #plt.get_yaxis().set_visible(False)
        plt.xticks(fontsize=fs)

        #ax = fig.add_subplot(gs[3, 2])

        plt.sca(axs[2])

        im = plt.imshow(rot_v_4, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xlabel("z [Mpc]",fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        #plt.get_yaxis().set_visible(False)
        plt.xticks(fontsize=fs)

        #ax = fig.add_subplot(gs[3, 3])

        plt.sca(axs[3])

        im = plt.imshow(rot_v_2, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xlabel("z [Mpc]",fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        #plt.get_yaxis().set_visible(False)
        plt.xticks(fontsize=fs)

        #plt.sca(axs[4])

        #im = plt.imshow(rot_v_0, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        #plt.xlabel("z [Mpc]",fontsize=fs)
        #plt.yticks([])
        #plt.grid(b=None)
        #plt.get_yaxis().set_visible(False)
        #plt.xticks(fontsize=fs)

        #cb = fig.colorbar(im, ax=ax,fraction=0.03, pad=0.14)#,pad = 0.04 )  # , orientation="horizontal")
        #cb.set_label(r'$\omega_x [km.s^{-1}.kpc^{-1}]$', size='large')

        cb = f.colorbar(im, ax=axs.ravel().tolist(),label=r'$\omega_x~[km~s^{-1}~kpc^{-1}]$', shrink=0.5,
                                pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/rot_v_filament_left_2.png', dpi=600, format='png')

        print('fig saved')




        # Adjust the layout to prevent overlap
        #plt.tight_layout()

        #plt.subplots_adjust(wspace=0.01, hspace=0.01)

        print("showing plot")

        # Show the plot
        plt.show()

        sys.exit()

    def right_fil_T_plot_2():
        print("loading maps")

        #t8 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_T_x_d0.02.bin', 'x')
        t6 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_6Mpc_map_T_x_d0.02.bin', 'x')
        t4 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_4Mpc_map_T_x_d0.02.bin', 'x')
        t2 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_2Mpc_map_T_x_d0.02.bin', 'x')
        #t0 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_0Mpc_map_T_x_d0.02.bin', 'x')

        print("maps loaded")

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        f, axs = plt.subplots(1, 3, figsize=(16, 8), constrained_layout=True, facecolor='white')
        # gs = gridspec.GridSpec(4, 4, wspace=0.001, hspace=0)

        fs = 14

        #title = "$\Delta x_{cen}$= -8Mpc"

        # First row : temperature

        # axs[0, 0].plot(x, np.sin(x))
        # ax = fig.add_subplot(gs[0, 0])
        # plt.sca(axs[0, 0])
        #plt.sca(axs[0])
        #title = "$\Delta x_{cen}$= 0Mpc"
        #plt.title(title, fontsize=fs)
        #t0 /= (kb / 1.602e-16)
        #plt.imshow(t0, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        #plt.xticks([])
        #plt.ylabel("y [Mpc]", fontsize=fs)
        #plt.yticks(fontsize=fs)
        #plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 1])
        # plt.sca(axs[0, 1])
        plt.sca(axs[0])
        title = "$\Delta x_{cen}$= 2Mpc"
        plt.title(title, fontsize=fs)
        t2 /= (kb / 1.602e-16)
        plt.imshow(t2, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        plt.xticks([])
        plt.ylabel("y [Mpc]", fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 2])
        # plt.sca(axs[0, 2])
        plt.sca(axs[1])
        title = "$\Delta x_{cen}$= 4Mpc"
        plt.title(title, fontsize=fs)
        t4 /= (kb / 1.602e-16)
        plt.imshow(t4, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 3])
        plt.sca(axs[2])
        title = "$\Delta x_{cen}$= 6Mpc"
        plt.title(title, fontsize=fs)
        t6 /= (kb / 1.602e-16)
        im = plt.imshow(t6, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.14)  # , orientation="horizontal")
        # cb.set_label('T[K]', size='large')
        #plt.sca(axs[4])
        #title = "$\Delta x_{cen}$= 0Mpc"
        #plt.title(title, fontsize=fs)
        #t0 /= (kb / 1.602e-16)
        #im = plt.imshow(t0, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        #plt.xticks([])
        #plt.yticks([])
        #plt.grid(b=None)

        # ratiov = t0.shape[0] / t0.shape[1]

        # print("ratiov",ratiov)

        # sys.exit()

        cb = f.colorbar(im, ax=axs.ravel().tolist(), label='T[k]', shrink=0.5, pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)
        cb.set_label('T[k]', size=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/T_filament_right_2.png', dpi=600, format='png')

        print("plot saved")

        plt.show()

        sys.exit()

    def right_fil_ne_plot_2():
        print("loading maps")

        #t8 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_T_x_d0.02.bin', 'x')
        ne6 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_6Mpc_map_ne_x_d0.02.bin', 'x')
        ne4 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_4Mpc_map_ne_x_d0.02.bin', 'x')
        ne2 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_2Mpc_map_ne_x_d0.02.bin', 'x')
        #t0 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_0Mpc_map_T_x_d0.02.bin', 'x')

        print("maps loaded")

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        f, axs = plt.subplots(1, 3, figsize=(16, 8), constrained_layout=True, facecolor='white')
        # gs = gridspec.GridSpec(4, 4, wspace=0.001, hspace=0)

        fs = 14

        #title = "$\Delta x_{cen}$= -8Mpc"

        #plt.imshow(ne6, cmap="cividis", origin='lower', alpha=1, extent=dim,
        #           norm=colors.LogNorm(vmin=10 ** (-7.5), vmax=10 ** (-4)))
        # plt.xticks([])
        #plt.xlabel("x [Mpc]", fontsize=fs)
        #plt.xticks(fontsize=fs)
        #plt.ylabel("y [Mpc]", fontsize=fs)
        #plt.yticks(fontsize=fs)
        # plt.yticks([])

        # First row : temperature

        # axs[0, 0].plot(x, np.sin(x))
        # ax = fig.add_subplot(gs[0, 0])
        # plt.sca(axs[0, 0])
        #plt.sca(axs[0])
        #title = "$\Delta x_{cen}$= 0Mpc"
        #plt.title(title, fontsize=fs)
        #t0 /= (kb / 1.602e-16)
        #plt.imshow(t0, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        #plt.xticks([])
        #plt.ylabel("y [Mpc]", fontsize=fs)
        #plt.yticks(fontsize=fs)
        #plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 1])
        # plt.sca(axs[0, 1])
        plt.sca(axs[0])
        title = "$\Delta x_{cen}$= 2Mpc"
        plt.title(title, fontsize=fs)
        #title = "$\Delta x_{cen}$= 2Mpc"
        #plt.title(title, fontsize=fs)
        #t2 /= (kb / 1.602e-16)
        plt.imshow(ne2, cmap="Solar", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=10 ** (-7.5), vmax=10 ** (-3)))
        plt.xlabel("x [Mpc]", fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.ylabel("y [Mpc]", fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 2])
        # plt.sca(axs[0, 2])
        plt.sca(axs[1])
        title = "$\Delta x_{cen}$= 4Mpc"
        plt.title(title, fontsize=fs)
        #title = "$\Delta x_{cen}$= 4Mpc"
        #plt.title(title, fontsize=fs)
        #t4 /= (kb / 1.602e-16)
        plt.imshow(ne4, cmap="Solar", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=10 ** (-7.5), vmax=10 ** (-3)))
        plt.xlabel("x [Mpc]", fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # ax = fig.add_subplot(gs[0, 3])
        plt.sca(axs[2])
        title = "$\Delta x_{cen}$= 6Mpc"
        plt.title(title, fontsize=fs)
        #title = "$\Delta x_{cen}$= 6Mpc"
        #plt.title(title, fontsize=fs)
        #t6 /= (kb / 1.602e-16)
        im = plt.imshow(ne6, cmap="Solar", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=10 ** (-7.5), vmax=10 ** (-3)))
        plt.xlabel("x [Mpc]", fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.14)  # , orientation="horizontal")
        # cb.set_label('T[K]', size='large')
        #plt.sca(axs[4])
        #title = "$\Delta x_{cen}$= 0Mpc"
        #plt.title(title, fontsize=fs)
        #t0 /= (kb / 1.602e-16)
        #im = plt.imshow(t0, cmap="magma", origin='lower', alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e5, vmax=1e8))
        #plt.xticks([])
        #plt.yticks([])
        #plt.grid(b=None)

        # ratiov = t0.shape[0] / t0.shape[1]

        # print("ratiov",ratiov)

        # sys.exit()

        cb = f.colorbar(im, ax=axs.ravel().tolist(), label='T[k]', shrink=0.5, pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)
        cb.set_label('$n_e~\mathrm{[cm^{-3}]}$', size=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/ne_filament_right_2.png', dpi=600, format='png')

        print("plot saved")

        plt.show()

        sys.exit()

    def right_fil_v_plot_2():

        #vz8 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_8Mpc_map_vz_x_d0.02.bin', 'x')
        vz6 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_6Mpc_map_vz_x_d0.02.bin', 'x')
        vz4 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_4Mpc_map_vz_x_d0.02.bin', 'x')
        vz2 = load_map('./maps/high_res/filament/map_high_19_xyz_left_plus_2Mpc_map_vz_x_d0.02.bin', 'x')
        #vz0 = load_map('./maps/high_res/filament/map_high_19_xyz_left_minus_0Mpc_map_vz_x_d0.02.bin', 'x')

        print("vz maps loaded")

        vz_virgo = -131.9249

        #vz8 -= vz_virgo
        vz6 -= vz_virgo
        vz4 -= vz_virgo
        vz2 -= vz_virgo
        #vz0 -= vz_virgo

        #mean10_map_8 = np.load("stream_vz_left_fil_minus_8Mpc.npy")
        #mean10_map2_8 = np.load("stream_vy_left_fil_minus_8Mpc.npy")
        mean10_map_6 = np.load("stream_vz_right_fil_plus_6Mpc.npy")
        mean10_map2_6 = np.load("stream_vy_right_fil_plus_6Mpc.npy")
        mean10_map_4 = np.load("stream_vz_right_fil_plus_4Mpc.npy")
        mean10_map2_4 = np.load("stream_vy_right_fil_plus_4Mpc.npy")
        mean10_map_2 = np.load("stream_vz_right_fil_plus_2Mpc.npy")
        mean10_map2_2 = np.load("stream_vy_right_fil_plus_2Mpc.npy")
        #mean10_map_0 = np.load("stream_vz_left_fil_minus_0Mpc.npy")
        #mean10_map0_2 = np.load("stream_vy_left_fil_minus_0Mpc.npy")

        x = np.load("mesh_x.npy")
        y = np.load("mesh_y.npy")

        print("streamlines loaded")

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        # Create a figure and a 4x4 grid of subplots
        # fig, axs = plt.subplots(4, 4, figsize=(16, 16), facecolor='white')
        # fig = plt.figure(figsize=(16, 16), facecolor='white', constrained_layout=True)

        # Second row : vz + streams

        f, axs = plt.subplots(1, 3, figsize=(16, 8), constrained_layout=True, facecolor='white')

        fs = 14

        #map_mean_mesh_8 = np.sqrt(mean10_map_8 ** 2 + mean10_map2_8 ** 2)
        #map_mean_mesh_8[np.isnan(map_mean_mesh_8) == True] = 0

        #norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 0])
        #plt.sca(axs[0])

        #stream = plt.streamplot(x, y, mean10_map_8, mean10_map2_8, linewidth=0.3, color=map_mean_mesh_8, cmap='Greys_r',
        #                        density=[1, 1.5], norm=norm)
        #im = plt.imshow(vz8, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
        #                vmax=1600)  # BR #origin='lower'
        #plt.xticks([])
        #plt.ylabel("y [Mpc]", fontsize=fs)
        #plt.yticks(fontsize=fs)
        #plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)

        map_mean_mesh_6 = np.sqrt(mean10_map_6 ** 2 + mean10_map2_6 ** 2)
        map_mean_mesh_6[np.isnan(map_mean_mesh_6) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 1])

        plt.sca(axs[2])

        stream = plt.streamplot(x, y, mean10_map_6, mean10_map2_6, linewidth=0.4, color=map_mean_mesh_6, cmap='Greys_r',
                                density=[2, 3], norm=norm)
        im = plt.imshow(vz6, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        map_mean_mesh_4 = np.sqrt(mean10_map_4 ** 2 + mean10_map2_4 ** 2)
        map_mean_mesh_4[np.isnan(map_mean_mesh_4) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 2])
        plt.sca(axs[1])

        stream = plt.streamplot(x, y, mean10_map_4, mean10_map2_4, linewidth=0.4, color=map_mean_mesh_4, cmap='Greys_r',
                                density=[2, 3], norm=norm)
        im = plt.imshow(vz4, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        map_mean_mesh_2 = np.sqrt(mean10_map_2 ** 2 + mean10_map2_2 ** 2)
        map_mean_mesh_2[np.isnan(map_mean_mesh_2) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 3])
        plt.sca(axs[0])

        stream = plt.streamplot(x, y, mean10_map_2, mean10_map2_2, linewidth=0.4, color=map_mean_mesh_2, cmap='Greys_r',
                                density=[2, 3], norm=norm)
        im = plt.imshow(vz2, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks([])
        plt.ylabel("y [Mpc]", fontsize=fs)
        plt.yticks(fontsize=fs)

        plt.grid(b=None)

        #plt.sca(axs[0])

        #map_mean_mesh_0 = np.sqrt(mean10_map_0 ** 2 + mean10_map0_2 ** 2)
        #map_mean_mesh_0[np.isnan(map_mean_mesh_0) == True] = 0

        #stream = plt.streamplot(x, y, mean10_map_0, mean10_map0_2, linewidth=0.4, color=map_mean_mesh_0, cmap='Greys_r',
        #                        density=[1, 1.5], norm=norm)
        #im = plt.imshow(vz0, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
        #                vmax=1600)  # BR #origin='lower'
        #plt.xticks([])
        #plt.yticks([])
        #plt.grid(b=None)
        #plt.ylabel("y [Mpc]", fontsize=fs)
        #plt.yticks(fontsize=fs)
        # ax.get_xaxis().set_visible(False)
        # ax.get_yaxis().set_visible(False)

        # cb = fig.colorbar(im, ax=ax,pad = 0.14, fraction=0.03 )  # , orientation="horizontal")
        # cb.set_label('vz[km/s]', size='large')

        # cbs = fig.colorbar(stream.lines, ax=ax, orientation='vertical', pad=0.24, fraction=0.03)
        # cbs.set_label(r'$||\vec{v}||~[km/s]$', size='large')

        cb = f.colorbar(im, ax=axs.ravel().tolist(), label='$v_z~\mathrm{[km~s^{-1}]}$', shrink=0.5, pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)

        #cb = f.colorbar(stream.lines, ax=axs.ravel().tolist(), label=r'$||\vec{v}||~[km/s]$', shrink=0.7, pad=0.01)  # ,size='small')#,shrink=0.6)
        #cb.ax.tick_params(labelsize=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/v_filament_right_2.png', dpi=600, format='png')

        print("plot saved")

        plt.show()
        sys.exit()

    def right_fil_div_plot_2():
        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        # div_v_8 = np.load("./filament/left_fil_transverse_vel_ops/div_v_minus_8Mpc.npy")
        div_v_6 = np.load("./filament/right_fil_transverse_vel_ops/div_v_plus_6Mpc.npy")
        div_v_4 = np.load("./filament/right_fil_transverse_vel_ops/div_v_plus_4Mpc.npy")
        div_v_2 = np.load("./filament/right_fil_transverse_vel_ops/div_v_plus_2Mpc.npy")
        # div_v_0 = np.load("./filament/left_fil_transverse_vel_ops/div_v_minus_0Mpc.npy")

        print("div_v maps loaded")

        # div_v_8 = div_v_8.T
        # div_v_8 = np.flip(div_v_8, axis=0)
        div_v_6 = div_v_6.T
        div_v_6 = np.flip(div_v_6, axis=0)
        div_v_4 = div_v_4.T
        div_v_4 = np.flip(div_v_4, axis=0)
        div_v_2 = div_v_2.T
        div_v_2 = np.flip(div_v_2, axis=0)
        # div_v_0 = div_v_0.T
        # div_v_0 = np.flip(div_v_0, axis=0)

        # Third row : div_v

        # ax = fig.add_subplot(gs[2, 0])

        f, axs = plt.subplots(1, 3, figsize=(16, 8), constrained_layout=True, facecolor='white')

        fs = 14

        # plt.sca(axs[0])

        # im = plt.imshow(div_v_0, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        # plt.xticks([])
        # plt.ylabel("y [Mpc]", fontsize=fs)
        # plt.yticks(fontsize=fs)
        # plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)

        # ax = fig.add_subplot(gs[2, 1])

        plt.sca(axs[0])

        im = plt.imshow(div_v_2, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks([])
        plt.ylabel("y [Mpc]", fontsize=fs)
        plt.yticks(fontsize=fs)
        # plt.yticks([])
        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # ax = fig.add_subplot(gs[2, 2])

        plt.sca(axs[1])

        im = plt.imshow(div_v_4, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks([])
        plt.yticks([])
        plt.grid(b=None)
        # ax.get_xaxis().set_visible(False)
        # ax.get_yaxis().set_visible(False)

        # ax = fig.add_subplot(gs[2, 3])

        plt.sca(axs[2])

        im = plt.imshow(div_v_6, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks([])
        plt.yticks([])

        plt.grid(b=None)
        # plt.get_xaxis().set_visible(False)
        # plt.get_yaxis().set_visible(False)

        # cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.14)#,pad = 0.04 )  # , orientation="horizontal")
        # cb.set_label(r'$\vec{\nabla}.\vec{v}[km.s^{-1}.kpc^{-1}]$', size='large')

        # plt.sca(axs[4])

        # im = plt.imshow(div_v_0, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        # plt.xticks([])
        # plt.yticks([])

        # plt.grid(b=None)

        cb = f.colorbar(im, ax=axs.ravel().tolist(), label=r'$\vec{\nabla}.\vec{v}~[km~s^{-1}~kpc^{-1}]$',
                        shrink=0.5,
                        pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/div_v_filament_right_2.png', dpi=600, format='png')

        print('fig saved')

        # Fourth row : rot_v

        # ax = fig.add_subplot(gs[3, 0])

        plt.show()
        sys.exit()

    def right_fil_rot_v_plot_2():
        # rot_v_8 = np.load("./filament/right_fil_transverse_vel_ops/rot_vx_minus_8Mpc.npy")
        rot_v_6 = np.load("./filament/right_fil_transverse_vel_ops/rot_vx_plus_6Mpc.npy")
        rot_v_4 = np.load("./filament/right_fil_transverse_vel_ops/rot_vx_plus_4Mpc.npy")
        rot_v_2 = np.load("./filament/right_fil_transverse_vel_ops/rot_vx_plus_2Mpc.npy")
        #rot_v_0 = np.load("./filament/left_fil_transverse_vel_ops/rot_vx_minus_0Mpc.npy")

        print("rot_v maps loaded")

        # rot_v_8 = rot_v_8.T
        # rot_v_8 = np.flip(rot_v_8, axis=0)
        rot_v_6 = rot_v_6.T
        rot_v_6 = np.flip(rot_v_6, axis=0)
        rot_v_4 = rot_v_4.T
        rot_v_4 = np.flip(rot_v_4, axis=0)
        rot_v_2 = rot_v_2.T
        rot_v_2 = np.flip(rot_v_2, axis=0)
        #rot_v_0 = rot_v_0.T
        #rot_v_0 = np.flip(rot_v_0, axis=0)

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dim = [-dimx, dimx, -dimy, dimy]

        fs = 14

        f, axs = plt.subplots(1, 3, figsize=(16, 8), constrained_layout=True, facecolor='white')

        #plt.sca(axs[0])

        #im = plt.imshow(rot_v_0, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        #plt.xlabel("z [Mpc]", fontsize=fs)
        #plt.ylabel("y [Mpc]", fontsize=fs)
        #plt.yticks(fontsize=fs)
        #plt.xticks(fontsize=fs)
        #plt.grid(b=None)

        # ax = fig.add_subplot(gs[3, 1])

        plt.sca(axs[0])

        im = plt.imshow(rot_v_2, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xlabel("z [Mpc]", fontsize=fs)
        plt.ylabel("y [Mpc]", fontsize=fs)
        plt.yticks(fontsize=fs)
        #plt.yticks([])
        plt.grid(b=None)
        # plt.get_yaxis().set_visible(False)
        plt.xticks(fontsize=fs)

        # ax = fig.add_subplot(gs[3, 2])

        plt.sca(axs[1])

        im = plt.imshow(rot_v_4, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xlabel("z [Mpc]", fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_yaxis().set_visible(False)
        plt.xticks(fontsize=fs)

        # ax = fig.add_subplot(gs[3, 3])

        plt.sca(axs[2])

        im = plt.imshow(rot_v_6, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xlabel("z [Mpc]", fontsize=fs)
        plt.yticks([])
        plt.grid(b=None)
        # plt.get_yaxis().set_visible(False)
        plt.xticks(fontsize=fs)

        # plt.sca(axs[4])

        # im = plt.imshow(rot_v_0, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        # plt.xlabel("z [Mpc]", fontsize=fs)
        # plt.yticks([])
        # plt.grid(b=None)
        # plt.get_yaxis().set_visible(False)
        # plt.xticks(fontsize=fs)

        # cb = fig.colorbar(im, ax=ax,fraction=0.03, pad=0.14)#,pad = 0.04 )  # , orientation="horizontal")
        # cb.set_label(r'$\omega_x [km.s^{-1}.kpc^{-1}]$', size='large')

        cb = f.colorbar(im, ax=axs.ravel().tolist(), label=r'$\omega_x~[km~s^{-1}~kpc^{-1}]$', shrink=0.5,
                        pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/rot_v_filament_right_2.png', dpi=600, format='png')

        print('fig saved')

        # Adjust the layout to prevent overlap
        # plt.tight_layout()

        # plt.subplots_adjust(wspace=0.01, hspace=0.01)

        print("showing plot")

        # Show the plot
        plt.show()

        sys.exit()





    #right_fil_T_plot()

    #right_fil_v_plot()

    #two_fil_v_plot()

    #right_fil_rot_v_plot()

    #right_fil_div_plot()

    #two_fil_v_plot()

    #left_fil_rot_v_plot_2()

    right_fil_ne_plot_2()

    sys.exit()

#large_fil_plots()

def zoom_fil_plot():

    def load_map(file, proj):
        h = FortranFile(file, 'r')

        nx, ny, nz = h.read_ints()
        cen_x, cen_y, cen_z = h.read_reals()

        if proj == "x":
            ncell = nz * ny
        elif proj == "y":
            ncell = nx * nz
        elif proj == "z":
            ncell = nx * ny

        map = np.zeros(ncell)

        map = ftp.f90_to_py.read_map_file(ncell, file, 0)

        if proj == "x":
            map = np.reshape(map, (nz, ny))
            # map2 = np.reshape(map2, (nz, ny))
            # map3 = np.reshape(map3, (nz, ny))
            nl = nx
        elif proj == "y":
            map = np.reshape(map, (nx, nz))
            # map2 = np.reshape(map2, (nx, nz))
            # map3 = np.reshape(map3, (nx, nz))
            nl = ny
        elif proj == "z":
            map = np.reshape(map, (ny, nx))
            # map2 = np.reshape(map2, (ny, nx))
            # map3 = np.reshape(map3, (ny, nx))
            nl = nz

        return map

    def left_fil_zoom():

        vx = load_map('./maps/high_res/filament/map_high_19_xyz_full_map_vx_z_d0.02.bin', 'z')

        nx,ny = np.shape(vx)

        print("nx,ny",nx,ny)


        lvl = 19

        cf = (737441*0.03)/nx

        ylim = 6000/cf
        #shifty = 2000/cf

        print("ylim",ylim)

        vx = vx[int(ny/2-ylim):int(ny/2+ylim),0:int(ny/2)]

        vx_virgo = -509.1301

        vx -= vx_virgo

        mean10_map = np.load("stream_vx_full_map.npy")
        mean10_map2 = np.load("stream_vy_full_map.npy")

        print("shape mean10_map",np.shape(mean10_map))

        ylim_mean10 = ylim/10
        #shifty10 = shifty/10


        mean10_map = mean10_map[int(ny/20-ylim_mean10):int(ny/20+ylim_mean10),0:int(nx/20)]
        mean10_map2 = mean10_map2[int(ny/20-ylim_mean10):int(ny/20+ylim_mean10),0:int(nx/20)]

        #print("shape vx",np.shape(vx))

        print("shape mean10_map",np.shape(mean10_map))

        #sys.exit()

        div_v = np.load("./filament/div_v_full_map.npy")

        div_v = div_v.T
        div_v = np.flip(div_v, axis=0)
        div_v = div_v[int(ny/2-ylim):int(ny/2+ylim),0:int(ny/2)]

        #print("shape div",np.shape(div_v))

        #sys.exit()

        rot_v = np.load("./filament/rot_vz_full_map.npy")

        rot_v = rot_v.T
        rot_v = np.flip(rot_v, axis=0)
        rot_v = rot_v[int(ny/2-ylim):int(ny/2+ylim),0:int(ny/2)]


        print("maps loaded")

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dimy = 6

        dim = [-dimx, 0, -dimy, dimy]

        x = np.load("mesh_x.npy")
        y = np.load("mesh_y.npy")

        print("min max x",np.min(x),np.max(x))
        print("min max y",np.min(y),np.max(y))

        #sys.exit()

        print("shape x",np.shape(x))

        x = x[int(ny/20-ylim_mean10):int(ny/20+ylim_mean10),0:int(nx/20)]
        y = y[int(ny/20-ylim_mean10):int(ny/20+ylim_mean10),0:int(nx/20)]

        print("min max x",np.min(x),np.max(x))
        print("min max y",np.min(y),np.max(y))

        print("shape x",np.shape(x))

        #sys.exit()

        def show_radii():
            xpx = np.linspace(-10, 0, 4000)
            ypx = np.linspace(-6, 6, 4000)

            Xpx, Ypx = np.meshgrid(xpx, ypx)

            # F = (Xpx-gal_pos[i,0])**2 + (Ypx-gal_pos[i,1])**2 - (r200dm[i]*(10/1e4))**2
            # F = (Xpx - x[i]) ** 2 + (Ypx - y[i]) ** 2 - (r200dm[i] * (10 / 1e4)) ** 2

            F = Xpx ** 2 + Ypx ** 2 - (1085 * (10 / 1e4)) ** 2  # r500
            r500 = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="solid", linewidths=1, alpha=1)
            plt.clabel(r500, r500.levels, inline=True, fmt="$R_{500}$", fontsize=14, manual=[(-0.4, -1)])

            F = Xpx ** 2 + Ypx ** 2 - (2024 * (10 / 1e4)) ** 2  # rvir
            rvir = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=1, alpha=1)
            plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{vir}$", fontsize=14, manual=[(-0.4, -2)])

            # F = Xpx ** 2 + Ypx ** 2 - (2895 * (10 / 1e4)) ** 2  # r200m
            # rvir = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=1, alpha=0.9)
            # plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{200m}$", fontsize=10)

            # if (type == "ne" or type == "T" or type == "P" or type == "dm" or type == "v"):
            #    F = Xpx ** 2 + Ypx ** 2 - (2024 * (10 / 1e4)) ** 2  # rvir
            #    rvir = plt.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9)
            #    plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{vir}$", fontsize=20)

            #    F = Xpx ** 2 + Ypx ** 2 - (5031 * (10 / 1e4)) ** 2  # rzv
            #    rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', ls='dashed', linewidths=2, alpha=0.9)
            #    plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{zv}$", fontsize=20)

            # if (type == "y" or type == "SD" or type == "EM"):
            #    F = Xpx ** 2 + Ypx ** 2 - (2895 * (10 / 1e4)) ** 2  # r200m
            #    rvir = plt.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9)
            #    plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{200m}$", fontsize=20)

            ##Rsp from baryons 3D rad prof

            # F = Xpx ** 2 + Ypx ** 2 - (4890 * (10 / 1e4)) ** 2  # rsb_relaxed bar
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dotted', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,rel}$", fontsize=20)

            # F = Xpx ** 2 + Ypx ** 2 - (3890 * (10 / 1e4)) ** 2  # rsb_spherical_collapse bar
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,sph}$", fontsize=20)

            ##Rsp from DM 3D rad prof

            # F = Xpx ** 2 + Ypx ** 2 - (4300 * (10 / 1e4)) ** 2  # rsb_relaxed DM
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', ls='dashed', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,rel}$", fontsize=20)

            # F = Xpx ** 2 + Ypx ** 2 - (3400 * (10 / 1e4)) ** 2  # rsb_spherical_collapse DM
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', ls='dotted', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,sph}$", fontsize=20)

            # if (type == "EM"):  ##Rsp from EM

            #    F = Xpx ** 2 + Ypx ** 2 - (3900 * (10 / 1e4)) ** 2
            #    rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            #    plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,EM}$", fontsize=20)

            # F = Xpx ** 2 + Ypx ** 2 - (3600 * (10 / 1e4)) ** 2
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,EM}$", fontsize=20)

            # F = Xpx ** 2 + Ypx ** 2 - (1098 * (10 / 1e4)) ** 2
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$test$", fontsize=20)

            # if (type == "y"):  ##Rsp from y_c

            #    F = Xpx ** 2 + Ypx ** 2 - (4900 * (10 / 1e4)) ** 2
            #    rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            #    plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,y_c}$", fontsize=20)

            # if (type == "SD"):  ##Rsp from Surface density

            #    F = Xpx ** 2 + Ypx ** 2 - (3300 * (10 / 1e4)) ** 2
            #    rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            #    plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,SD}$", fontsize=20)  # ,manual=True)

            # F = Xpx ** 2 + Ypx ** 2 - (8000 * (10 / 1e4)) ** 2  # background lower limit
            # rbgl = plt.contour(Xpx, Ypx, F, [0], colors='grey', linewidths=0.6, alpha=0.9)
            # plt.clabel(rbgl, rbgl.levels, inline=True, fmt="$Background lower limit}$", fontsize=10)

            # F = Xpx ** 2 + Ypx ** 2 - (10000 * (10 / 1e4)) ** 2  # background upper limit
            # rbgh = plt.contour(Xpx, Ypx, F, [0], colors='grey', linewidths=0.6, alpha=0.9)
            # plt.clabel(rbgh, rbgh.levels, inline=True, fmt="$Background upper limit}$", fontsize=10)

            # plt.fill(rbgl,rbgh)

            # F = Xpx ** 2 + Ypx ** 2 - (1648 * (10 / 1e4)) ** 2 #r200
            # r200=plt.contour(Xpx, Ypx, F, [0], colors='blue', linewidths=0.6, alpha=0.9)
            # plt.clabel(r200, r200.levels, inline=True, fmt="$R_{200}$", fontsize=10)
            # F = Xpx ** 2 + Ypx ** 2 - (1085 * (10 / 1e4)) ** 2 #r500
            # r500=plt.contour(Xpx, Ypx, F, [0], colors='blue',linestyles="dashed", linewidths=0.6, alpha=0.9)
            # plt.clabel(r500, r500.levels, inline=True, fmt="$R_{500}$", fontsize=10)
            # F = Xpx ** 2 + Ypx ** 2 - (900 * (10 / 1e4)) ** 2  # bump radius
            # rbump = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
            # plt.clabel(rbump, rbump.levels, inline=True, fmt="$R_{bump}$", fontsize=10)
            # F = Xpx ** 2 + Ypx ** 2 - (850 * (10 / 1e4)) ** 2  # bump radius
            # rbump = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
            # plt.clabel(rbump, rbump.levels, inline=True, fmt="$R_{bump}$", fontsize=10)

            # F = Xpx ** 2 + Ypx ** 2 - (400 * (10 / 1e4)) ** 2  # low mw deproj radius
            # rdec = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
            # plt.clabel(rdec, rdec.levels, inline=True, fmt="$R_{decrease in deproj prof}$", fontsize=10)
            # plt.scatter(x[i], y[i], alpha=0.9, s=1, c='white')
            # plt.show()
            # sys.exit()



        f, axs = plt.subplots(1, 3, figsize=(12, 6), constrained_layout=True, facecolor='white')

        fs = 14

        map_mean_mesh = np.sqrt(mean10_map ** 2 + mean10_map2 ** 2)
        map_mean_mesh[np.isnan(map_mean_mesh) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 1])

        plt.sca(axs[0])

        stream = plt.streamplot(x, y, mean10_map, mean10_map2, linewidth=0.3, color=map_mean_mesh, cmap='Greys_r',density=[4, 4], norm=norm)
        im = plt.imshow(vx, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks(fontsize=fs)
        #plt.yticks(fontsize=fs)
        plt.grid(b=None)
        plt.yticks([])
        plt.xlabel("x [Mpc]", fontsize=fs)
        #axs[0].xaxis.set_label_position('top')
        #axs[0].xaxis.tick_top()
        plt.xticks(fontsize=fs)
        plt.title("Velocity field", fontsize=fs)

        show_radii()
        #plt.ylabel("y [Mpc]", fontsize=fs)
        #plt.tick_params(which="minor", right=False, top=True,bottom=False,left=False)#, labelbottom=False) #direction="in"
        #plt.tick_params( which='major', right=False, top=True, bottom=False,left=False, labelsize=fs, labeltop=True,labelleft=False,labelbottom=False)#,labelright=True)

        cb = plt.colorbar(im, ax=axs[0], orientation='horizontal',shrink=0.8,pad=0.01)#,fraction=0.047 * ratioh)  # pad=0.04, shrink=0.95 ) #orientation='horizontal'
        cb.set_label(r'$v_x~\mathrm{[km~s^{-1}]}$', size=fs)  # , size='large')
        cb.ax.tick_params(labelsize=fs)

        #cbs = plt.colorbar(stream.lines, ax=axs.ravel().tolist(), location='left', shrink=0.8,pad=0.1)  # ,fraction=0.047 * ratioh)  # pad=0.04, shrink=0.95 ) #orientation='horizontal'
        #cbs.set_label(r'$||\vec{v}||~\mathrm{[km~s^{-1}]}$', size=fs)  # , size='large')
        #cbs.ax.tick_params(labelsize=fs)

        #divider = make_axes_locatable(axs[0])
        #cax = divider.new_vertical(size='5%', pad=0.5)
        #f.add_axes(cax)
        #f.colorbar(im, cax=cax, orientation='horizontal')

        plt.sca(axs[1])

        im = plt.imshow(div_v, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks(fontsize=fs)
        plt.yticks([])
        plt.xlabel("x [Mpc]", fontsize=fs)
        #axs[1].xaxis.set_label_position('top')
        #axs[1].xaxis.tick_top()
        plt.xticks(fontsize=fs)
        plt.grid(b=None)
        plt.title("Divergence", fontsize=fs)

        show_radii()
        #plt.tick_params(which="minor", right=False, top=True,bottom=False,left=False, length=3)#, labelbottom=False) #direction="in"
        #plt.tick_params(which="major", right=False, top=True, bottom=False,left=False, length=5, labelsize=18, labeltop=True)#,labelright=True)

        cb = plt.colorbar(im, ax=axs[1], orientation='horizontal', shrink=0.8,pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)
        cb.set_label(r'$\vec{\nabla}.\vec{v}~[km~s^{-1}~kpc^{-1}]$',size=fs)

        #cb = f.colorbar(im, ax=axs[0],orientation='horizontal')  # ,fraction=0.047 * ratioh)  # pad=0.04, shrink=0.95 ) #orientation='horizontal'
        #cb.set_label(r'$v_x~\mathrm{[km~s^{-1}]}$', size=fs)  # , size='large')
        #cb.ax.tick_params(labelsize=fs)

        plt.sca(axs[2])

        im = plt.imshow(rot_v, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        #plt.xticks(fontsize=fs)
        #plt.yticks([])
        plt.xlabel("x [Mpc]", fontsize=fs)
        #axs[2].xaxis.set_label_position('top')
        #axs[2].xaxis.tick_top()
        plt.ylabel ("y [Mpc]", fontsize=fs)
        axs[2].yaxis.set_label_position('right')
        axs[2].yaxis.tick_right()
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(b=None)
        plt.title("Vorticity", fontsize=fs)

        show_radii()
        #plt.tick_params(which="minor", right=True, top=True,bottom=False,left=False, length=3)#, labelbottom=False) #direction="in"
        #plt.tick_params(which="major", right=True, top=True, bottom=False,left=False, length=5, labelsize=18, labeltop=True,labelright=True)

        cb = plt.colorbar(im, ax=axs[2], orientation='horizontal', shrink=0.8,pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)
        cb.set_label(r'$\omega_z~[km~s^{-1}~kpc^{-1}]$',size=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/zoom_filament_left.png', dpi=600, format='png')

        print('fig saved')

        plt.show()

        sys.exit()

    def left_fil_zoom_2():

        mach = load_map('./maps/high_res/filament/map_high_19_xyz_full_map_mach_z_d0.02.bin', 'z')

        nx,ny = np.shape(mach)

        print("nx,ny",nx,ny)


        lvl = 19

        cf = (737441*0.03)/nx

        ylim = 6000/cf
        #shifty = 2000/cf

        print("ylim",ylim)

        mach = mach[int(ny/2-ylim):int(ny/2+ylim),0:int(ny/2)]

        #vx_virgo = -509.1301

        #vx -= vx_virgo

        #mean10_map = np.load("stream_vx_full_map.npy")
        #mean10_map2 = np.load("stream_vy_full_map.npy")

        #print("shape mean10_map",np.shape(mean10_map))

        ylim_mean10 = ylim/10
        #shifty10 = shifty/10


        #mean10_map = mean10_map[int(ny/20-ylim_mean10):int(ny/20+ylim_mean10),0:int(nx/20)]
        #mean10_map2 = mean10_map2[int(ny/20-ylim_mean10):int(ny/20+ylim_mean10),0:int(nx/20)]

        #print("shape vx",np.shape(vx))

        #print("shape mean10_map",np.shape(mean10_map))

        #sys.exit()

        ne = load_map('./maps/high_res/filament/map_high_19_xyz_full_map_ne_z_d0.02.bin', 'z')

        div_v = np.load("./filament/div_v_full_map.npy")

        div_v = div_v.T
        div_v = np.flip(div_v, axis=0)
        div_v *= ne
        div_v = div_v[int(ny/2-ylim):int(ny/2+ylim),0:int(ny/2)]

        #print("shape div",np.shape(div_v))

        #sys.exit()

        rot_v = np.load("./filament/rot_vz_full_map.npy")

        rot_v = rot_v.T
        rot_v = np.flip(rot_v, axis=0)
        rot_v *= ne
        rot_v = rot_v[int(ny/2-ylim):int(ny/2+ylim),0:int(ny/2)]


        print("maps loaded")

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dimy = 6

        dim = [-dimx, 0, -dimy, dimy]

        #x = np.load("mesh_x.npy")
        #y = np.load("mesh_y.npy")

        #print("min max x",np.min(x),np.max(x))
        #print("min max y",np.min(y),np.max(y))

        #sys.exit()

        #print("shape x",np.shape(x))

        #x = x[int(ny/20-ylim_mean10):int(ny/20+ylim_mean10),0:int(nx/20)]
        #y = y[int(ny/20-ylim_mean10):int(ny/20+ylim_mean10),0:int(nx/20)]

        #print("min max x",np.min(x),np.max(x))
        #print("min max y",np.min(y),np.max(y))

        #print("shape x",np.shape(x))

        #sys.exit()

        def show_radii():
            xpx = np.linspace(-10, 0, 4000)
            ypx = np.linspace(-6, 6, 4000)

            Xpx, Ypx = np.meshgrid(xpx, ypx)

            # F = (Xpx-gal_pos[i,0])**2 + (Ypx-gal_pos[i,1])**2 - (r200dm[i]*(10/1e4))**2
            # F = (Xpx - x[i]) ** 2 + (Ypx - y[i]) ** 2 - (r200dm[i] * (10 / 1e4)) ** 2

            F = Xpx ** 2 + Ypx ** 2 - (1085 * (10 / 1e4)) ** 2  # r500
            r500 = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="solid", linewidths=1, alpha=1)
            plt.clabel(r500, r500.levels, inline=True, fmt="$R_{500}$", fontsize=14, manual=[(-0.4, -1)])

            F = Xpx ** 2 + Ypx ** 2 - (2024 * (10 / 1e4)) ** 2  # rvir
            rvir = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=1, alpha=1)
            plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{vir}$", fontsize=14, manual=[(-0.4, -2)])

            # F = Xpx ** 2 + Ypx ** 2 - (2895 * (10 / 1e4)) ** 2  # r200m
            # rvir = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=1, alpha=0.9)
            # plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{200m}$", fontsize=10)

            # if (type == "ne" or type == "T" or type == "P" or type == "dm" or type == "v"):
            #    F = Xpx ** 2 + Ypx ** 2 - (2024 * (10 / 1e4)) ** 2  # rvir
            #    rvir = plt.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9)
            #    plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{vir}$", fontsize=20)

            #    F = Xpx ** 2 + Ypx ** 2 - (5031 * (10 / 1e4)) ** 2  # rzv
            #    rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', ls='dashed', linewidths=2, alpha=0.9)
            #    plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{zv}$", fontsize=20)

            # if (type == "y" or type == "SD" or type == "EM"):
            #    F = Xpx ** 2 + Ypx ** 2 - (2895 * (10 / 1e4)) ** 2  # r200m
            #    rvir = plt.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9)
            #    plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{200m}$", fontsize=20)

            ##Rsp from baryons 3D rad prof

            # F = Xpx ** 2 + Ypx ** 2 - (4890 * (10 / 1e4)) ** 2  # rsb_relaxed bar
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dotted', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,rel}$", fontsize=20)

            # F = Xpx ** 2 + Ypx ** 2 - (3890 * (10 / 1e4)) ** 2  # rsb_spherical_collapse bar
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,sph}$", fontsize=20)

            ##Rsp from DM 3D rad prof

            # F = Xpx ** 2 + Ypx ** 2 - (4300 * (10 / 1e4)) ** 2  # rsb_relaxed DM
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', ls='dashed', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,rel}$", fontsize=20)

            # F = Xpx ** 2 + Ypx ** 2 - (3400 * (10 / 1e4)) ** 2  # rsb_spherical_collapse DM
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', ls='dotted', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,sph}$", fontsize=20)

            # if (type == "EM"):  ##Rsp from EM

            #    F = Xpx ** 2 + Ypx ** 2 - (3900 * (10 / 1e4)) ** 2
            #    rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            #    plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,EM}$", fontsize=20)

            # F = Xpx ** 2 + Ypx ** 2 - (3600 * (10 / 1e4)) ** 2
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,EM}$", fontsize=20)

            # F = Xpx ** 2 + Ypx ** 2 - (1098 * (10 / 1e4)) ** 2
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$test$", fontsize=20)

            # if (type == "y"):  ##Rsp from y_c

            #    F = Xpx ** 2 + Ypx ** 2 - (4900 * (10 / 1e4)) ** 2
            #    rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            #    plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,y_c}$", fontsize=20)

            # if (type == "SD"):  ##Rsp from Surface density

            #    F = Xpx ** 2 + Ypx ** 2 - (3300 * (10 / 1e4)) ** 2
            #    rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            #    plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,SD}$", fontsize=20)  # ,manual=True)

            # F = Xpx ** 2 + Ypx ** 2 - (8000 * (10 / 1e4)) ** 2  # background lower limit
            # rbgl = plt.contour(Xpx, Ypx, F, [0], colors='grey', linewidths=0.6, alpha=0.9)
            # plt.clabel(rbgl, rbgl.levels, inline=True, fmt="$Background lower limit}$", fontsize=10)

            # F = Xpx ** 2 + Ypx ** 2 - (10000 * (10 / 1e4)) ** 2  # background upper limit
            # rbgh = plt.contour(Xpx, Ypx, F, [0], colors='grey', linewidths=0.6, alpha=0.9)
            # plt.clabel(rbgh, rbgh.levels, inline=True, fmt="$Background upper limit}$", fontsize=10)

            # plt.fill(rbgl,rbgh)

            # F = Xpx ** 2 + Ypx ** 2 - (1648 * (10 / 1e4)) ** 2 #r200
            # r200=plt.contour(Xpx, Ypx, F, [0], colors='blue', linewidths=0.6, alpha=0.9)
            # plt.clabel(r200, r200.levels, inline=True, fmt="$R_{200}$", fontsize=10)
            # F = Xpx ** 2 + Ypx ** 2 - (1085 * (10 / 1e4)) ** 2 #r500
            # r500=plt.contour(Xpx, Ypx, F, [0], colors='blue',linestyles="dashed", linewidths=0.6, alpha=0.9)
            # plt.clabel(r500, r500.levels, inline=True, fmt="$R_{500}$", fontsize=10)
            # F = Xpx ** 2 + Ypx ** 2 - (900 * (10 / 1e4)) ** 2  # bump radius
            # rbump = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
            # plt.clabel(rbump, rbump.levels, inline=True, fmt="$R_{bump}$", fontsize=10)
            # F = Xpx ** 2 + Ypx ** 2 - (850 * (10 / 1e4)) ** 2  # bump radius
            # rbump = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
            # plt.clabel(rbump, rbump.levels, inline=True, fmt="$R_{bump}$", fontsize=10)

            # F = Xpx ** 2 + Ypx ** 2 - (400 * (10 / 1e4)) ** 2  # low mw deproj radius
            # rdec = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
            # plt.clabel(rdec, rdec.levels, inline=True, fmt="$R_{decrease in deproj prof}$", fontsize=10)
            # plt.scatter(x[i], y[i], alpha=0.9, s=1, c='white')
            # plt.show()
            # sys.exit()



        f, axs = plt.subplots(1, 3, figsize=(12, 6), constrained_layout=True, facecolor='white')

        fs = 14

        #map_mean_mesh = np.sqrt(mean10_map ** 2 + mean10_map2 ** 2)
        #map_mean_mesh[np.isnan(map_mean_mesh) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 1])

        plt.sca(axs[0])

        map_above_1 = np.ma.masked_less(mach, 1)  # Mask values below 1
        map_below_1 = np.ma.masked_greater(mach, 1)

        im2 = plt.imshow(map_below_1, cmap="Mono_r", origin='lower', alpha=1, extent=dim, vmin=0.1,
                        vmax=1)  # , vmin=0, vmax=6)  "hot" "flare"

        im1 = plt.imshow(map_above_1, cmap="Blues", origin='lower', alpha=1, extent=dim,
                        norm=colors.LogNorm(vmin=1, vmax=100))  # , vmin=0, vmax=6)  "hot" "flare"
        # show_radii()
        # ax.axis('off')


        #ratio = map.shape[1] / map.shape[0]

        cb1 = f.colorbar(im1, ax=axs[0], orientation="horizontal", pad=0.01,
                         shrink=0.95)  # fraction=0.047 * ratio, , label='$T~\mathrm{[K]}$',size="large")#,pad=0.1)
        cb1.ax.tick_params(labelsize=14)
        cb1.set_label(r'$M=\frac{v}{c_s}$', size=14)

        #
        cb2 = f.colorbar(im2, ax=axs[0], orientation="horizontal", pad=0.01, shrink=0.95)#fraction=0.047 * ratio,
                           #pad=0.2)  # , label='$T~\mathrm{[K]}$',size="large")#,pad=0.1)
        cb2.ax.tick_params(labelsize=14)
        #cb2.set_label(r'$M=\frac{v}{c_s}$ (0.1 to 1, linear scale)', size=14)



        #stream = plt.streamplot(x, y, mean10_map, mean10_map2, linewidth=0.3, color=map_mean_mesh, cmap='Greys_r',density=[4, 4], norm=norm)
        #im = plt.imshow(mach, cmap="lajolla_r", origin="lower", alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e-1, vmax=20)) #vmin=-1600,vmax=1600)  # BR #origin='lower'
        plt.xticks(fontsize=fs)
        #plt.yticks(fontsize=fs)
        plt.grid(b=None)
        plt.yticks([])
        plt.xlabel("x [Mpc]", fontsize=fs)
        #axs[0].xaxis.set_label_position('top')
        #axs[0].xaxis.tick_top()
        plt.xticks(fontsize=fs)
        plt.title("Mach number", fontsize=fs)

        show_radii()
        #plt.ylabel("y [Mpc]", fontsize=fs)
        #plt.tick_params(which="minor", right=False, top=True,bottom=False,left=False)#, labelbottom=False) #direction="in"
        #plt.tick_params( which='major', right=False, top=True, bottom=False,left=False, labelsize=fs, labeltop=True,labelleft=False,labelbottom=False)#,labelright=True)

        #cb = plt.colorbar(im, ax=axs[0], orientation='horizontal',shrink=0.8,pad=0.01)#,fraction=0.047 * ratioh)  # pad=0.04, shrink=0.95 ) #orientation='horizontal'
        #cb.set_label(r'$M=\frac{v}{c_s}$', size=fs)  # , size='large')
        #cb.ax.tick_params(labelsize=fs)

        #cbs = plt.colorbar(stream.lines, ax=axs.ravel().tolist(), location='left', shrink=0.8,pad=0.1)  # ,fraction=0.047 * ratioh)  # pad=0.04, shrink=0.95 ) #orientation='horizontal'
        #cbs.set_label(r'$||\vec{v}||~\mathrm{[km~s^{-1}]}$', size=fs)  # , size='large')
        #cbs.ax.tick_params(labelsize=fs)

        #divider = make_axes_locatable(axs[0])
        #cax = divider.new_vertical(size='5%', pad=0.5)
        #f.add_axes(cax)
        #f.colorbar(im, cax=cax, orientation='horizontal')

        plt.sca(axs[1])

        im = plt.imshow(div_v, cmap="curl", extent=dim, norm=colors.SymLogNorm(linthresh=1e-8, linscale=0.2,vmin=-2, vmax=2, base=10))#vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks(fontsize=fs)
        plt.yticks([])
        plt.xlabel("x [Mpc]", fontsize=fs)
        #axs[1].xaxis.set_label_position('top')
        #axs[1].xaxis.tick_top()
        plt.xticks(fontsize=fs)
        plt.grid(b=None)
        plt.title("Electron density-weighted divergence", fontsize=fs)

        show_radii()
        #plt.tick_params(which="minor", right=False, top=True,bottom=False,left=False, length=3)#, labelbottom=False) #direction="in"
        #plt.tick_params(which="major", right=False, top=True, bottom=False,left=False, length=5, labelsize=18, labeltop=True)#,labelright=True)

        cb = plt.colorbar(im, ax=axs[1], orientation='horizontal', shrink=0.8,pad=0.01,ticks=[-1e0,-1e-3,-1e-6,1e-6,1e-3,1e0])  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)
        cb.set_label(r'$n_e \vec{\nabla}.\vec{v}~[cm^{-3}~km~s^{-1}~kpc^{-1}]$',size=fs)

        #cb = f.colorbar(im, ax=axs[0],orientation='horizontal')  # ,fraction=0.047 * ratioh)  # pad=0.04, shrink=0.95 ) #orientation='horizontal'
        #cb.set_label(r'$v_x~\mathrm{[km~s^{-1}]}$', size=fs)  # , size='large')
        #cb.ax.tick_params(labelsize=fs)

        plt.sca(axs[2])

        im = plt.imshow(rot_v, cmap="curl", extent=dim,norm=colors.SymLogNorm(linthresh=1e-8, linscale=0.2,vmin=-2, vmax=2, base=10)) #vmin=-20, vmax=20 # , origin='lower')
        #plt.xticks(fontsize=fs)
        #plt.yticks([])
        plt.xlabel("x [Mpc]", fontsize=fs)
        #axs[2].xaxis.set_label_position('top')
        #axs[2].xaxis.tick_top()
        plt.ylabel ("y [Mpc]", fontsize=fs)
        axs[2].yaxis.set_label_position('right')
        axs[2].yaxis.tick_right()
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(b=None)
        plt.title("Electron density-weighted vorticity", fontsize=fs)

        show_radii()
        #plt.tick_params(which="minor", right=True, top=True,bottom=False,left=False, length=3)#, labelbottom=False) #direction="in"
        #plt.tick_params(which="major", right=True, top=True, bottom=False,left=False, length=5, labelsize=18, labeltop=True,labelright=True)

        cb = plt.colorbar(im, ax=axs[2], orientation='horizontal', shrink=0.8,pad=0.01,ticks=[-1e0,-1e-3,-1e-6,1e-6,1e-3,1e0])  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)
        cb.set_label(r'$n_e \omega_z~[cm^{-3}~km~s^{-1}~kpc^{-1}]$',size=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/zoom_filament_left_mach-nedivv-nerotv.png', dpi=600, format='png')

        print('fig saved')

        plt.show()

        sys.exit()


    def right_fil_zoom():

        vx = load_map('./maps/high_res/filament/map_high_19_xyz_full_map_vx_z_d0.02.bin', 'z')

        nx, ny = np.shape(vx)

        print("nx,ny", nx, ny)

        lvl = 19

        cf = (737441 * 0.03) / nx

        ylim = 6000 / cf
        shifty = 2000 / cf

        print("ylim", ylim)

        vx = vx[int(ny / 2 - ylim - shifty):int(ny / 2 + ylim - shifty), int(ny / 2)::]

        vx_virgo = -509.1301

        vx -= vx_virgo

        mean10_map = np.load("stream_vx_full_map.npy")
        mean10_map2 = np.load("stream_vy_full_map.npy")

        print("shape mean10_map", np.shape(mean10_map))

        ylim_mean10 = ylim / 10
        shifty10 = shifty / 10

        mean10_map = mean10_map[int(ny / 20 - ylim_mean10 - shifty10):int(ny / 20 + ylim_mean10 - shifty10), int(ny / 20)::]
        mean10_map2 = mean10_map2[int(ny / 20 - ylim_mean10 - shifty10):int(ny / 20 + ylim_mean10 - shifty10), int(ny / 20)::]

        # print("shape vx",np.shape(vx))

        print("shape mean10_map", np.shape(mean10_map))

        # sys.exit()

        cf = (737441 * 0.03) / nx

        ylim = 6000 / cf
        shifty = 2000 / cf

        div_v = np.load("./filament/div_v_full_map.npy")

        div_v = div_v.T
        div_v = np.flip(div_v, axis=0)
        div_v = div_v[int(ny / 2 - ylim + shifty):int(ny / 2 + ylim + shifty), int(ny / 2)::]

        # print("shape div",np.shape(div_v))

        # sys.exit()

        rot_v = np.load("./filament/rot_vz_full_map.npy")

        rot_v = rot_v.T
        rot_v = np.flip(rot_v, axis=0)
        rot_v = rot_v[int(ny / 2 - ylim + shifty):int(ny / 2 + ylim + shifty), int(ny / 2)::]

        print("maps loaded")

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dimy = 6

        dim = [0,dimx, -dimy-2, dimy-2]

        x = np.load("mesh_x.npy")
        y = np.load("mesh_y.npy")

        print("min max x", np.min(x), np.max(x))
        print("min max y", np.min(y), np.max(y))

        # sys.exit()

        print("shape x", np.shape(x))

        x = x[int(ny / 20 - ylim_mean10 - shifty10):int(ny / 20 + ylim_mean10 - shifty10), int(nx / 20)::]
        y = y[int(ny / 20 - ylim_mean10 - shifty10):int(ny / 20 + ylim_mean10 - shifty10), int(nx / 20)::]

        print("min max x", np.min(x), np.max(x))
        print("min max y", np.min(y), np.max(y))

        print("shape x", np.shape(x))

        # sys.exit()

        def show_radii():
            xpx = np.linspace(0, 10, 4000)
            ypx = np.linspace(-8, 4, 4000)

            Xpx, Ypx = np.meshgrid(xpx, ypx)

            # F = (Xpx-gal_pos[i,0])**2 + (Ypx-gal_pos[i,1])**2 - (r200dm[i]*(10/1e4))**2
            # F = (Xpx - x[i]) ** 2 + (Ypx - y[i]) ** 2 - (r200dm[i] * (10 / 1e4)) ** 2

            F = Xpx ** 2 + Ypx ** 2 - (1085 * (10 / 1e4)) ** 2  # r500
            r500 = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="solid", linewidths=1, alpha=1)
            plt.clabel(r500, r500.levels, inline=True, fmt="$R_{500}$", fontsize=10, manual=[(0.4, -1)])

            F = Xpx ** 2 + Ypx ** 2 - (2024 * (10 / 1e4)) ** 2  # rvir
            rvir = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=1, alpha=1)
            plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{vir}$", fontsize=10, manual=[(0.4, -2)])

            # F = Xpx ** 2 + Ypx ** 2 - (2895 * (10 / 1e4)) ** 2  # r200m
            # rvir = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=1, alpha=0.9)
            # plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{200m}$", fontsize=10)

            # if (type == "ne" or type == "T" or type == "P" or type == "dm" or type == "v"):
            #    F = Xpx ** 2 + Ypx ** 2 - (2024 * (10 / 1e4)) ** 2  # rvir
            #    rvir = plt.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9)
            #    plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{vir}$", fontsize=20)

            #    F = Xpx ** 2 + Ypx ** 2 - (5031 * (10 / 1e4)) ** 2  # rzv
            #    rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', ls='dashed', linewidths=2, alpha=0.9)
            #    plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{zv}$", fontsize=20)

            # if (type == "y" or type == "SD" or type == "EM"):
            #    F = Xpx ** 2 + Ypx ** 2 - (2895 * (10 / 1e4)) ** 2  # r200m
            #    rvir = plt.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9)
            #    plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{200m}$", fontsize=20)

            ##Rsp from baryons 3D rad prof

            # F = Xpx ** 2 + Ypx ** 2 - (4890 * (10 / 1e4)) ** 2  # rsb_relaxed bar
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dotted', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,rel}$", fontsize=20)

            # F = Xpx ** 2 + Ypx ** 2 - (3890 * (10 / 1e4)) ** 2  # rsb_spherical_collapse bar
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,sph}$", fontsize=20)

            ##Rsp from DM 3D rad prof

            # F = Xpx ** 2 + Ypx ** 2 - (4300 * (10 / 1e4)) ** 2  # rsb_relaxed DM
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', ls='dashed', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,rel}$", fontsize=20)

            # F = Xpx ** 2 + Ypx ** 2 - (3400 * (10 / 1e4)) ** 2  # rsb_spherical_collapse DM
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', ls='dotted', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,sph}$", fontsize=20)

            # if (type == "EM"):  ##Rsp from EM

            #    F = Xpx ** 2 + Ypx ** 2 - (3900 * (10 / 1e4)) ** 2
            #    rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            #    plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,EM}$", fontsize=20)

            # F = Xpx ** 2 + Ypx ** 2 - (3600 * (10 / 1e4)) ** 2
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,EM}$", fontsize=20)

            # F = Xpx ** 2 + Ypx ** 2 - (1098 * (10 / 1e4)) ** 2
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$test$", fontsize=20)

            # if (type == "y"):  ##Rsp from y_c

            #    F = Xpx ** 2 + Ypx ** 2 - (4900 * (10 / 1e4)) ** 2
            #    rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            #    plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,y_c}$", fontsize=20)

            # if (type == "SD"):  ##Rsp from Surface density

            #    F = Xpx ** 2 + Ypx ** 2 - (3300 * (10 / 1e4)) ** 2
            #    rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            #    plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,SD}$", fontsize=20)  # ,manual=True)

            # F = Xpx ** 2 + Ypx ** 2 - (8000 * (10 / 1e4)) ** 2  # background lower limit
            # rbgl = plt.contour(Xpx, Ypx, F, [0], colors='grey', linewidths=0.6, alpha=0.9)
            # plt.clabel(rbgl, rbgl.levels, inline=True, fmt="$Background lower limit}$", fontsize=10)

            # F = Xpx ** 2 + Ypx ** 2 - (10000 * (10 / 1e4)) ** 2  # background upper limit
            # rbgh = plt.contour(Xpx, Ypx, F, [0], colors='grey', linewidths=0.6, alpha=0.9)
            # plt.clabel(rbgh, rbgh.levels, inline=True, fmt="$Background upper limit}$", fontsize=10)

            # plt.fill(rbgl,rbgh)

            # F = Xpx ** 2 + Ypx ** 2 - (1648 * (10 / 1e4)) ** 2 #r200
            # r200=plt.contour(Xpx, Ypx, F, [0], colors='blue', linewidths=0.6, alpha=0.9)
            # plt.clabel(r200, r200.levels, inline=True, fmt="$R_{200}$", fontsize=10)
            # F = Xpx ** 2 + Ypx ** 2 - (1085 * (10 / 1e4)) ** 2 #r500
            # r500=plt.contour(Xpx, Ypx, F, [0], colors='blue',linestyles="dashed", linewidths=0.6, alpha=0.9)
            # plt.clabel(r500, r500.levels, inline=True, fmt="$R_{500}$", fontsize=10)
            # F = Xpx ** 2 + Ypx ** 2 - (900 * (10 / 1e4)) ** 2  # bump radius
            # rbump = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
            # plt.clabel(rbump, rbump.levels, inline=True, fmt="$R_{bump}$", fontsize=10)
            # F = Xpx ** 2 + Ypx ** 2 - (850 * (10 / 1e4)) ** 2  # bump radius
            # rbump = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
            # plt.clabel(rbump, rbump.levels, inline=True, fmt="$R_{bump}$", fontsize=10)

            # F = Xpx ** 2 + Ypx ** 2 - (400 * (10 / 1e4)) ** 2  # low mw deproj radius
            # rdec = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
            # plt.clabel(rdec, rdec.levels, inline=True, fmt="$R_{decrease in deproj prof}$", fontsize=10)
            # plt.scatter(x[i], y[i], alpha=0.9, s=1, c='white')
            # plt.show()
            # sys.exit()

        f, axs = plt.subplots(1, 3, figsize=(12, 6), constrained_layout=True, facecolor='white')

        fs = 14

        map_mean_mesh = np.sqrt(mean10_map ** 2 + mean10_map2 ** 2)
        map_mean_mesh[np.isnan(map_mean_mesh) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 1])

        plt.sca(axs[0])

        stream = plt.streamplot(x, y, mean10_map, mean10_map2, linewidth=0.3, color=map_mean_mesh, cmap='Greys_r',
                                density=[4, 4], norm=norm)
        im = plt.imshow(vx, cmap="BR", origin="lower", alpha=1, extent=dim, vmin=-1600,
                        vmax=1600)  # BR #origin='lower'
        plt.xticks(fontsize=fs)
        # plt.yticks(fontsize=fs)
        plt.grid(b=None)
        plt.yticks([])
        plt.xlabel("x [Mpc]", fontsize=fs)
        #axs[0].xaxis.set_label_position('top')
        #axs[0].xaxis.tick_top()
        plt.title("Velocity field", fontsize=fs)
        plt.xticks(fontsize=fs)

        show_radii()
        # plt.ylabel("y [Mpc]", fontsize=fs)
        # plt.tick_params(which="minor", right=False, top=True,bottom=False,left=False)#, labelbottom=False) #direction="in"
        # plt.tick_params( which='major', right=False, top=True, bottom=False,left=False, labelsize=fs, labeltop=True,labelleft=False,labelbottom=False)#,labelright=True)

        cb = plt.colorbar(im, ax=axs[0], orientation='horizontal', shrink=0.8,
                          pad=0.01)  # ,fraction=0.047 * ratioh)  # pad=0.04, shrink=0.95 ) #orientation='horizontal'
        cb.set_label(r'$v_x~\mathrm{[km~s^{-1}]}$', size=fs)  # , size='large')
        cb.ax.tick_params(labelsize=fs)

        # cbs = plt.colorbar(stream.lines, ax=axs.ravel().tolist(), location='left', shrink=0.8,pad=0.1)  # ,fraction=0.047 * ratioh)  # pad=0.04, shrink=0.95 ) #orientation='horizontal'
        # cbs.set_label(r'$||\vec{v}||~\mathrm{[km~s^{-1}]}$', size=fs)  # , size='large')
        # cbs.ax.tick_params(labelsize=fs)

        # divider = make_axes_locatable(axs[0])
        # cax = divider.new_vertical(size='5%', pad=0.5)
        # f.add_axes(cax)
        # f.colorbar(im, cax=cax, orientation='horizontal')

        plt.sca(axs[1])

        im = plt.imshow(div_v, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        plt.xticks(fontsize=fs)
        plt.yticks([])
        plt.xlabel("x [Mpc]", fontsize=fs)
        #axs[1].xaxis.set_label_position('top')
        #axs[1].xaxis.tick_top()
        plt.xticks(fontsize=fs)
        plt.grid(b=None)
        plt.title("Divergence", fontsize=fs)

        show_radii()
        # plt.tick_params(which="minor", right=False, top=True,bottom=False,left=False, length=3)#, labelbottom=False) #direction="in"
        # plt.tick_params(which="major", right=False, top=True, bottom=False,left=False, length=5, labelsize=18, labeltop=True)#,labelright=True)

        cb = plt.colorbar(im, ax=axs[1], orientation='horizontal', shrink=0.8, pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)
        cb.set_label(r'$\vec{\nabla}.\vec{v}~[km~s^{-1}~kpc^{-1}]$', size=fs)

        # cb = f.colorbar(im, ax=axs[0],orientation='horizontal')  # ,fraction=0.047 * ratioh)  # pad=0.04, shrink=0.95 ) #orientation='horizontal'
        # cb.set_label(r'$v_x~\mathrm{[km~s^{-1}]}$', size=fs)  # , size='large')
        # cb.ax.tick_params(labelsize=fs)

        plt.sca(axs[2])

        im = plt.imshow(rot_v, cmap="seismic", vmin=-20, vmax=20, extent=dim)  # , origin='lower')
        # plt.xticks(fontsize=fs)
        # plt.yticks([])
        plt.xlabel("x [Mpc]", fontsize=fs)
        #axs[2].xaxis.set_label_position('top')
        #axs[2].xaxis.tick_top()
        plt.ylabel("y [Mpc]", fontsize=fs)
        axs[2].yaxis.set_label_position('right')
        axs[2].yaxis.tick_right()
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(b=None)
        plt.title("Vorticity", fontsize=fs)

        show_radii()
        # plt.tick_params(which="minor", right=True, top=True,bottom=False,left=False, length=3)#, labelbottom=False) #direction="in"
        # plt.tick_params(which="major", right=True, top=True, bottom=False,left=False, length=5, labelsize=18, labeltop=True,labelright=True)

        cb = plt.colorbar(im, ax=axs[2], orientation='horizontal', shrink=0.8, pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)
        cb.set_label(r'$\omega_z~[km~s^{-1}~kpc^{-1}]$', size=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/zoom_filament_right.png', dpi=600, format='png')

        print('fig saved')

        plt.show()

        sys.exit()

    def right_fil_zoom_2():

        mach = load_map('./maps/high_res/filament/map_high_19_xyz_full_map_mach_z_d0.02.bin', 'z')

        nx, ny = np.shape(mach)

        print("nx,ny", nx, ny)

        lvl = 19

        cf = (737441 * 0.03) / nx

        ylim = 6000 / cf
        shifty = 2000 / cf

        print("ylim", ylim)

        mach = mach[int(ny / 2 - ylim - shifty):int(ny / 2 + ylim - shifty), int(ny / 2)::]

        #vx_virgo = -509.1301

        #vx -= vx_virgo

        #mean10_map = np.load("stream_vx_full_map.npy")
        #mean10_map2 = np.load("stream_vy_full_map.npy")

        #print("shape mean10_map", np.shape(mean10_map))

        ylim_mean10 = ylim / 10
        shifty10 = shifty / 10

        #mean10_map = mean10_map[int(ny / 20 - ylim_mean10 - shifty10):int(ny / 20 + ylim_mean10 - shifty10), int(ny / 20)::]
        #mean10_map2 = mean10_map2[int(ny / 20 - ylim_mean10 - shifty10):int(ny / 20 + ylim_mean10 - shifty10), int(ny / 20)::]

        # print("shape vx",np.shape(vx))

        #print("shape mean10_map", np.shape(mean10_map))

        # sys.exit()

        cf = (737441 * 0.03) / nx

        ylim = 6000 / cf
        shifty = 2000 / cf

        file="./maps/high_res/filament/map_high_19_xyz_full_map_ne_z_d0.02.bin"
        h = FortranFile(file, 'r')
        nx, ny, nz = h.read_ints()
        ncell = nx * ny
        print("ncell", ncell)
        ne = np.zeros(ncell)
        ne = ftp.f90_to_py.read_map_file(ncell, file, 0)
        ne = np.reshape(ne, (nx, ny))

        #ne = np.load('./maps/high_res/filament/map_high_19_xyz_full_map_ne_z_d0.02.bin')

        div_v = np.load("./filament/div_v_full_map.npy")

        div_v = div_v.T
        div_v = np.flip(div_v, axis=0)
        div_v *= ne
        div_v = div_v[int(ny / 2 - ylim + shifty):int(ny / 2 + ylim + shifty), int(ny / 2)::]

        #print("shape div_v",np.shape(div_v))
        #print('min & max div',np.min(div_v),np.max(div_v))

        #plt.hist(div_v.flatten(), bins=100)#, range=(-20,20))
        #plt.xscale('log')
        #plt.show()

        #sys.exit()

        # print("shape div",np.shape(div_v))

        # sys.exit()

        rot_v = np.load("./filament/rot_vz_full_map.npy")

        rot_v = rot_v.T
        rot_v = np.flip(rot_v, axis=0)
        rot_v *= ne
        rot_v = rot_v[int(ny / 2 - ylim + shifty):int(ny / 2 + ylim + shifty), int(ny / 2)::]

        print("maps loaded")

        nx = 15728
        lvl = 19

        dimx = (nx / 2) * (737.441 / 2 ** lvl)
        dimy = dimx

        dimy = 6

        dim = [0,dimx, -dimy-2, dimy-2]

        #x = np.load("mesh_x.npy")
        #y = np.load("mesh_y.npy")

        #print("min max x", np.min(x), np.max(x))
        #print("min max y", np.min(y), np.max(y))

        # sys.exit()

        #print("shape x", np.shape(x))

        #x = x[int(ny / 20 - ylim_mean10 - shifty10):int(ny / 20 + ylim_mean10 - shifty10), int(nx / 20)::]
        #y = y[int(ny / 20 - ylim_mean10 - shifty10):int(ny / 20 + ylim_mean10 - shifty10), int(nx / 20)::]

        #print("min max x", np.min(x), np.max(x))
        #print("min max y", np.min(y), np.max(y))

        #print("shape x", np.shape(x))

        # sys.exit()

        def show_radii():
            xpx = np.linspace(0, 10, 4000)
            ypx = np.linspace(-8, 4, 4000)

            Xpx, Ypx = np.meshgrid(xpx, ypx)

            # F = (Xpx-gal_pos[i,0])**2 + (Ypx-gal_pos[i,1])**2 - (r200dm[i]*(10/1e4))**2
            # F = (Xpx - x[i]) ** 2 + (Ypx - y[i]) ** 2 - (r200dm[i] * (10 / 1e4)) ** 2

            F = Xpx ** 2 + Ypx ** 2 - (1085 * (10 / 1e4)) ** 2  # r500
            r500 = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="solid", linewidths=1, alpha=1)
            plt.clabel(r500, r500.levels, inline=True, fmt="$R_{500}$", fontsize=10, manual=[(0.4, -1)])

            F = Xpx ** 2 + Ypx ** 2 - (2024 * (10 / 1e4)) ** 2  # rvir
            rvir = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=1, alpha=1)
            plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{vir}$", fontsize=10, manual=[(0.4, -2)])

            # F = Xpx ** 2 + Ypx ** 2 - (2895 * (10 / 1e4)) ** 2  # r200m
            # rvir = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=1, alpha=0.9)
            # plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{200m}$", fontsize=10)

            # if (type == "ne" or type == "T" or type == "P" or type == "dm" or type == "v"):
            #    F = Xpx ** 2 + Ypx ** 2 - (2024 * (10 / 1e4)) ** 2  # rvir
            #    rvir = plt.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9)
            #    plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{vir}$", fontsize=20)

            #    F = Xpx ** 2 + Ypx ** 2 - (5031 * (10 / 1e4)) ** 2  # rzv
            #    rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', ls='dashed', linewidths=2, alpha=0.9)
            #    plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{zv}$", fontsize=20)

            # if (type == "y" or type == "SD" or type == "EM"):
            #    F = Xpx ** 2 + Ypx ** 2 - (2895 * (10 / 1e4)) ** 2  # r200m
            #    rvir = plt.contour(Xpx, Ypx, F, [0], colors='white', linewidths=2, alpha=0.9)
            #    plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{200m}$", fontsize=20)

            ##Rsp from baryons 3D rad prof

            # F = Xpx ** 2 + Ypx ** 2 - (4890 * (10 / 1e4)) ** 2  # rsb_relaxed bar
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dotted', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,rel}$", fontsize=20)

            # F = Xpx ** 2 + Ypx ** 2 - (3890 * (10 / 1e4)) ** 2  # rsb_spherical_collapse bar
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,sph}$", fontsize=20)

            ##Rsp from DM 3D rad prof

            # F = Xpx ** 2 + Ypx ** 2 - (4300 * (10 / 1e4)) ** 2  # rsb_relaxed DM
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', ls='dashed', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,rel}$", fontsize=20)

            # F = Xpx ** 2 + Ypx ** 2 - (3400 * (10 / 1e4)) ** 2  # rsb_spherical_collapse DM
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', ls='dotted', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,sph}$", fontsize=20)

            # if (type == "EM"):  ##Rsp from EM

            #    F = Xpx ** 2 + Ypx ** 2 - (3900 * (10 / 1e4)) ** 2
            #    rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            #    plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,EM}$", fontsize=20)

            # F = Xpx ** 2 + Ypx ** 2 - (3600 * (10 / 1e4)) ** 2
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,EM}$", fontsize=20)

            # F = Xpx ** 2 + Ypx ** 2 - (1098 * (10 / 1e4)) ** 2
            # rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            # plt.clabel(rzv, rzv.levels, inline=True, fmt="$test$", fontsize=20)

            # if (type == "y"):  ##Rsp from y_c

            #    F = Xpx ** 2 + Ypx ** 2 - (4900 * (10 / 1e4)) ** 2
            #    rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            #    plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,y_c}$", fontsize=20)

            # if (type == "SD"):  ##Rsp from Surface density

            #    F = Xpx ** 2 + Ypx ** 2 - (3300 * (10 / 1e4)) ** 2
            #    rzv = plt.contour(Xpx, Ypx, F, [0], colors='white', linestyles='dashed', linewidths=2, alpha=0.9)
            #    plt.clabel(rzv, rzv.levels, inline=True, fmt="$R_{sp,SD}$", fontsize=20)  # ,manual=True)

            # F = Xpx ** 2 + Ypx ** 2 - (8000 * (10 / 1e4)) ** 2  # background lower limit
            # rbgl = plt.contour(Xpx, Ypx, F, [0], colors='grey', linewidths=0.6, alpha=0.9)
            # plt.clabel(rbgl, rbgl.levels, inline=True, fmt="$Background lower limit}$", fontsize=10)

            # F = Xpx ** 2 + Ypx ** 2 - (10000 * (10 / 1e4)) ** 2  # background upper limit
            # rbgh = plt.contour(Xpx, Ypx, F, [0], colors='grey', linewidths=0.6, alpha=0.9)
            # plt.clabel(rbgh, rbgh.levels, inline=True, fmt="$Background upper limit}$", fontsize=10)

            # plt.fill(rbgl,rbgh)

            # F = Xpx ** 2 + Ypx ** 2 - (1648 * (10 / 1e4)) ** 2 #r200
            # r200=plt.contour(Xpx, Ypx, F, [0], colors='blue', linewidths=0.6, alpha=0.9)
            # plt.clabel(r200, r200.levels, inline=True, fmt="$R_{200}$", fontsize=10)
            # F = Xpx ** 2 + Ypx ** 2 - (1085 * (10 / 1e4)) ** 2 #r500
            # r500=plt.contour(Xpx, Ypx, F, [0], colors='blue',linestyles="dashed", linewidths=0.6, alpha=0.9)
            # plt.clabel(r500, r500.levels, inline=True, fmt="$R_{500}$", fontsize=10)
            # F = Xpx ** 2 + Ypx ** 2 - (900 * (10 / 1e4)) ** 2  # bump radius
            # rbump = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
            # plt.clabel(rbump, rbump.levels, inline=True, fmt="$R_{bump}$", fontsize=10)
            # F = Xpx ** 2 + Ypx ** 2 - (850 * (10 / 1e4)) ** 2  # bump radius
            # rbump = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
            # plt.clabel(rbump, rbump.levels, inline=True, fmt="$R_{bump}$", fontsize=10)

            # F = Xpx ** 2 + Ypx ** 2 - (400 * (10 / 1e4)) ** 2  # low mw deproj radius
            # rdec = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
            # plt.clabel(rdec, rdec.levels, inline=True, fmt="$R_{decrease in deproj prof}$", fontsize=10)
            # plt.scatter(x[i], y[i], alpha=0.9, s=1, c='white')
            # plt.show()
            # sys.exit()

        f, axs = plt.subplots(1, 3, figsize=(12, 6), constrained_layout=True, facecolor='white')

        fs = 14

        #map_mean_mesh = np.sqrt(mean10_map ** 2 + mean10_map2 ** 2)
        #map_mean_mesh[np.isnan(map_mean_mesh) == True] = 0

        norm = Normalize(vmin=0, vmax=1200)

        # ax = fig.add_subplot(gs[1, 1])

        plt.sca(axs[0])

        #stream = plt.streamplot(x, y, mean10_map, mean10_map2, linewidth=0.3, color=map_mean_mesh, cmap='Greys_r',
        #                        density=[4, 4], norm=norm)
        #plt.sca(axs[0])

        map_above_1 = np.ma.masked_less(mach, 1)  # Mask values below 1
        map_below_1 = np.ma.masked_greater(mach, 1)

        im2 = plt.imshow(map_below_1, cmap="Mono_r", origin='lower', alpha=1, extent=dim, vmin=0.1,
                         vmax=1)  # , vmin=0, vmax=6)  "hot" "flare"

        im1 = plt.imshow(map_above_1, cmap="Blues", origin='lower', alpha=1, extent=dim,
                         norm=colors.LogNorm(vmin=1, vmax=100))  # , vmin=0, vmax=6)  "hot" "flare"
        # show_radii()
        # ax.axis('off')

        # ratio = map.shape[1] / map.shape[0]

        cb1 = f.colorbar(im1, ax=axs[0], orientation="horizontal", pad=0.01,
                         shrink=0.95)  # fraction=0.047 * ratio, , label='$T~\mathrm{[K]}$',size="large")#,pad=0.1)
        cb1.ax.tick_params(labelsize=14)
        cb1.set_label(r'$M=\frac{v}{c_s}$', size=14)

        #
        cb2 = f.colorbar(im2, ax=axs[0], orientation="horizontal", pad=0.01, shrink=0.95)  # fraction=0.047 * ratio,
        # pad=0.2)  # , label='$T~\mathrm{[K]}$',size="large")#,pad=0.1)
        cb2.ax.tick_params(labelsize=14)
        # cb2.set_label(r'$M=\frac{v}{c_s}$ (0.1 to 1, linear scale)', size=14)

        # stream = plt.streamplot(x, y, mean10_map, mean10_map2, linewidth=0.3, color=map_mean_mesh, cmap='Greys_r',density=[4, 4], norm=norm)
        # im = plt.imshow(mach, cmap="lajolla_r", origin="lower", alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e-1, vmax=20)) #vmin=-1600,vmax=1600)  # BR #origin='lower'
        plt.xticks(fontsize=fs)
        # plt.yticks(fontsize=fs)
        plt.grid(b=None)
        plt.yticks([])
        plt.xlabel("x [Mpc]", fontsize=fs)
        # axs[0].xaxis.set_label_position('top')
        # axs[0].xaxis.tick_top()
        plt.xticks(fontsize=fs)
        plt.title("Mach number", fontsize=fs)

        show_radii()
        #im = plt.imshow(mach, cmap="lajolla_r", origin="lower", alpha=1, extent=dim, norm=colors.LogNorm(vmin=1e-1, vmax=20))#vmin=0,vmax=20)  # BR #origin='lower'
        #plt.xticks(fontsize=fs)
        # plt.yticks(fontsize=fs)
        #plt.grid(b=None)
        #plt.yticks([])
        #plt.xlabel("x [Mpc]", fontsize=fs)
        #axs[0].xaxis.set_label_position('top')
        #axs[0].xaxis.tick_top()
        #plt.title("Mach number", fontsize=fs)
        #plt.xticks(fontsize=fs)

        #show_radii()
        # plt.ylabel("y [Mpc]", fontsize=fs)
        # plt.tick_params(which="minor", right=False, top=True,bottom=False,left=False)#, labelbottom=False) #direction="in"
        # plt.tick_params( which='major', right=False, top=True, bottom=False,left=False, labelsize=fs, labeltop=True,labelleft=False,labelbottom=False)#,labelright=True)

        #cb = plt.colorbar(im, ax=axs[0], orientation='horizontal', shrink=0.8,
        #                  pad=0.01)  # ,fraction=0.047 * ratioh)  # pad=0.04, shrink=0.95 ) #orientation='horizontal'
        #cb.set_label(r'$M=\frac{v}{c_s}$', size=fs)  # , size='large')
        #cb.ax.tick_params(labelsize=fs)

        # cbs = plt.colorbar(stream.lines, ax=axs.ravel().tolist(), location='left', shrink=0.8,pad=0.1)  # ,fraction=0.047 * ratioh)  # pad=0.04, shrink=0.95 ) #orientation='horizontal'
        # cbs.set_label(r'$||\vec{v}||~\mathrm{[km~s^{-1}]}$', size=fs)  # , size='large')
        # cbs.ax.tick_params(labelsize=fs)

        # divider = make_axes_locatable(axs[0])
        # cax = divider.new_vertical(size='5%', pad=0.5)
        # f.add_axes(cax)
        # f.colorbar(im, cax=cax, orientation='horizontal')

        plt.sca(axs[1])

        im = plt.imshow(div_v, cmap="curl", extent=dim, norm=colors.SymLogNorm(linthresh=1e-8, linscale=0.2,vmin=-2, vmax=2, base=10))#norm=colors.LogNorm(vmin=1e-1, vmax=20))vmin=-2, vmax=2, extent=dim)  # , origin='lower') cmap="seismic"
        plt.xticks(fontsize=fs)
        plt.yticks([])
        plt.xlabel("x [Mpc]", fontsize=fs)
        #axs[1].xaxis.set_label_position('top')
        #axs[1].xaxis.tick_top()
        plt.xticks(fontsize=fs)
        plt.grid(b=None)
        plt.title("Electron density-weighted divergence", fontsize=fs)

        show_radii()
        # plt.tick_params(which="minor", right=False, top=True,bottom=False,left=False, length=3)#, labelbottom=False) #direction="in"
        # plt.tick_params(which="major", right=False, top=True, bottom=False,left=False, length=5, labelsize=18, labeltop=True)#,labelright=True)

        cb = plt.colorbar(im, ax=axs[1], orientation='horizontal', shrink=0.8, pad=0.01,ticks=[-1e0,-1e-3,-1e-6,1e-6,1e-3,1e0])  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)
        cb.set_label(r'$n_e\vec{\nabla}.\vec{v}~[cm^{-3}km~s^{-1}~kpc^{-1}]$', size=fs)

        # cb = f.colorbar(im, ax=axs[0],orientation='horizontal')  # ,fraction=0.047 * ratioh)  # pad=0.04, shrink=0.95 ) #orientation='horizontal'
        # cb.set_label(r'$v_x~\mathrm{[km~s^{-1}]}$', size=fs)  # , size='large')
        # cb.ax.tick_params(labelsize=fs)

        plt.sca(axs[2])

        im = plt.imshow(rot_v, cmap="curl", extent=dim, norm=colors.SymLogNorm(linthresh=1e-8, linscale=0.2,vmin=-2, vmax=2, base=10))  # , origin='lower')
        # plt.xticks(fontsize=fs)
        # plt.yticks([])
        plt.xlabel("x [Mpc]", fontsize=fs)
        #axs[2].xaxis.set_label_position('top')
        #axs[2].xaxis.tick_top()
        plt.ylabel("y [Mpc]", fontsize=fs)
        axs[2].yaxis.set_label_position('right')
        axs[2].yaxis.tick_right()
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(b=None)
        plt.title("Electron density-weighted vorticity", fontsize=fs)

        show_radii()
        # plt.tick_params(which="minor", right=True, top=True,bottom=False,left=False, length=3)#, labelbottom=False) #direction="in"
        # plt.tick_params(which="major", right=True, top=True, bottom=False,left=False, length=5, labelsize=18, labeltop=True,labelright=True)

        cb = plt.colorbar(im, ax=axs[2], orientation='horizontal', shrink=0.8, pad=0.01,ticks=[-1e0,-1e-3,-1e-6,1e-6,1e-3,1e0])  # ,size='small')#,shrink=0.6)
        cb.ax.tick_params(labelsize=fs)
        cb.set_label(r'$n_e\omega_z~[cm^{-3}~km~s^{-1}~kpc^{-1}]$', size=fs)

        plt.savefig('./papers_plots/paper_3_turb_fil/zoom_filament_right_mach-nedivv-nerotv.png', dpi=600, format='png')

        print('fig saved')

        print('showing fig')

        plt.show()

        sys.exit()

    #left_fil_zoom_2()

    left_fil_zoom_2()

zoom_fil_plot()


def show_decomposed_vel_fields():
    def read_vel_matrix(file):

        h = FortranFile(file, 'r')
        nx, ny, nz = h.read_ints()
        print("nx,ny,nz", nx, ny, nz)
        ncell = nx * ny * nz
        print("ncell", ncell)
        # sys.exit()
        map = np.zeros(ncell)

        map = ftp.f90_to_py.read_map_file(ncell, file, 1)

        print("reshape map")

        map = np.reshape(map, (nx, ny, nz), order='F')

        print("reshape done")

        map = np.array(map, dtype=np.float32)

        return map

    # vx = read_vel_matrix('./maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_vx_2Mpc_cube.bin')
    # vy = read_vel_matrix('./maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_vy_2Mpc_cube.bin')
    # vz = read_vel_matrix('./maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_vz_2Mpc_cube.bin')

    lvl = 15

    if lvl == 16:
        # vx = read_vel_matrix('./maps/high_res/filament/3D_fields/map_3D_high_16_xyz_right_minus_0Mpc_vx_10Mpc_cube.bin')
        # vy = read_vel_matrix('./maps/high_res/filament/3D_fields/map_3D_high_16_xyz_right_minus_0Mpc_vy_10Mpc_cube.bin')
        # vz = read_vel_matrix('./maps/high_res/filament/3D_fields/map_3D_high_16_xyz_right_minus_0Mpc_vz_10Mpc_cube.bin')

        # v = np.load('./maps/high_res/filament/3D_fields/map_3D_high_16_xyz_right_minus_0Mpc_v_10Mpc_cube.bin.npy')
        # v_comp = np.load('./maps/high_res/filament/3D_fields/map_3D_high_16_xyz_right_minus_0Mpc_v_comp_10Mpc_cube.bin.npy')
        # v_sol = np.load('./maps/high_res/filament/3D_fields/map_3D_high_16_xyz_right_minus_0Mpc_v_sol_10Mpc_cube.bin.npy')

        v = np.load('./maps/high_res/filament/3D_fields/map_3D_high_16_xyz_right_minus_0Mpc_v_10x10x20Mpc_cube.bin.npy')
        v_comp = np.load(
            './maps/high_res/filament/3D_fields/map_3D_high_16_xyz_right_minus_0Mpc_v_comp_10x10x20Mpc_cube.bin.npy')
        v_sol = np.load(
            './maps/high_res/filament/3D_fields/map_3D_high_16_xyz_right_minus_0Mpc_v_sol_10x10x20Mpc_cube.bin.npy')

    if lvl == 15:
        # vx = read_vel_matrix('./maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_vx_10Mpc_cube.bin')
        # vy = read_vel_matrix('./maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_vy_10Mpc_cube.bin')
        # vz = read_vel_matrix('./maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_vz_10Mpc_cube.bin')

        # vx = read_vel_matrix('./maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_vx_10x10x20Mpc_cube_rdr15.bin')
        # vy = read_vel_matrix('./maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_vy_10x10x20Mpc_cube_rdr15.bin')
        # vz = read_vel_matrix('./maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_vz_10x10x20Mpc_cube_rdr15.bin')

        # v = np.load('./maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_v_10Mpc_cube.bin.npy')
        # v_comp = np.load('./maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_v_comp_10Mpc_cube.bin.npy')
        # v_sol = np.load('./maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_v_sol_10Mpc_cube.bin.npy')

        v = np.load('./maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_v_12x12x20Mpc_cube.bin.npy')
        v_comp = np.load(
            './maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_v_comp_12x12x20Mpc_cube.bin.npy')
        v_sol = np.load(
            './maps/high_res/filament/3D_fields/map_3D_high_15_xyz_right_minus_0Mpc_v_sol_12x12x20Mpc_cube.bin.npy')

        # v = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)

    nx, ny, nz = np.shape(v)

    print("nx,ny,nz",nx,ny,nz)

    #sys.exit()

    # lvl = 15

    dimx = (ny / 2) * (737.441 / 2 ** lvl)
    dimy = dimx

    dim = [-dimx, dimx, -dimy, dimy]

    def show_map_3(vela, titrea, velb, titreb, velc, titrec, comp):
        slice = int(nx / 2)

        mapa = vela[:, :, slice]
        mapb = velb[:, :, slice]
        mapc = velc[:, :, slice]

        mapa = mapa.T
        # mapa = np.flip(mapa, axis=0)

        mapb = mapb.T
        # mapb = np.flip(mapb, axis=0)

        mapc = mapc.T
        # mapc = np.flip(mapc, axis=0)

        if comp == "comp":

            vmin = -1000
            vmax = 1000

            colormap = "bwr"

        elif comp == "tot":

            vmin = 0
            vmax = 1000

            colormap = "cividis"

        f, axs = plt.subplots(1, 3, figsize=(20, 70), constrained_layout=True)
        plt.sca(axs[0])
        im = plt.imshow(mapa, cmap=colormap, origin="lower", alpha=1, extent=dim, vmin=vmin, vmax=vmax)
        plt.title(titrea)
        plt.xlabel("x [Mpc]")
        plt.ylabel("y [Mpc]")
        # plt.colorbar()

        # vmin = 0
        # vmax = 1
        # colormap = "binary"

        plt.sca(axs[1])
        im = plt.imshow(mapb, cmap=colormap, origin="lower", alpha=1, extent=dim, vmin=vmin, vmax=vmax)
        plt.title(titreb)
        plt.yticks([])
        plt.xlabel("x [Mpc]")

        # vmin = 0
        # vmax = 1
        # colormap = "binary"

        plt.sca(axs[2])
        im = plt.imshow(mapc, cmap=colormap, origin="lower", alpha=1, extent=dim, vmin=vmin, vmax=vmax)
        plt.title(titrec)
        plt.yticks([])
        plt.xlabel("x [Mpc]")

        f.colorbar(im, ax=axs.ravel().tolist(), label='velocity $[\mathrm{km\,s^{-1}}]$')  # ,shrink=0.6)

        # plt.subplots_adjust(wspace=0.05)#, hspace=0.3)

        plt.show()

    def show_map_nx3(vela, titrea, velb, titreb, velc, titrec, comp):
        slice = int(nx / 2)

        print("slice", slice)

        print('nx', nx)

        # slice = int(nx/2-nx/5)

        print("slice", slice)

        mapa = vela[slice, :, :]
        mapb = mapa
        mapc = mapa
        # mapb = velb[slice, :, :]
        # mapc = velc[slice, :, :]

        # mapa = vela[:, slice, :]
        # mapb = velb[:, slice, :]
        # mapc = velc[:, slice, :]

        # mapa = np.mean(vela[slice-1:slice, :, :],axis=0)

        # mapa = mapa.T
        # mapa = np.flip(mapa, axis=0)

        # mapb = mapb.T
        # mapb = np.flip(mapb, axis=0)

        # mapc = mapc.T
        # mapc = np.flip(mapc, axis=0)

        if comp == "comp":

            vmin = -1000
            vmax = 1000

            colormap = "bwr"

        elif comp == "tot":

            vmin = 0
            vmax = 1200

            colormap = "cividis"

        init = 0
        n = 8

        f, axs = plt.subplots(3, n, figsize=(12, 6), constrained_layout=True)
        for i in range(init, init + n):
            print("i", i)

            fs = 13

            title = "$\Delta x_{cen}$=" + str(round(2 * i - 8, 3)) + " Mpc"

            slice = int(nx / 2 - (4 - i) * (nx / 10))

            # print("slice", slice)

            # print('nx', nx)

            # slice = int(nx/2-nx/5)

            print("slice", slice)

            mapa = vela[slice, :, :]
            # mapb = mapa
            # mapc = mapa

            mapb = velb[slice, :, :]
            mapc = velc[slice, :, :]

            plt.sca(axs[0, i - init])
            im = plt.imshow(mapa, cmap=colormap, origin="lower", alpha=1, extent=dim, vmin=vmin, vmax=vmax)
            plt.title(title, fontsize=fs)
            # plt.xlabel("x [Mpc]")
            if i - init == 0:
                plt.ylabel("y [Mpc]", fontsize=fs)
                plt.yticks(fontsize=fs)
            else:
                plt.yticks([])
            # plt.xticks([])
            a = plt.gca()
            xax = a.axes.get_xaxis()
            xax = xax.set_visible(False)

            plt.grid(b=None)

            if i == (init + n - 1):
                plt.text(1.05, 0.5, 'Total', transform=axs[0, n - 1].transAxes, verticalalignment='center',
                         rotation=270, fontsize=fs)
            # plt.colorbar()

            # vmin = 0
            # vmax = 1
            # colormap = "binary"

            plt.sca(axs[1, i - init])
            im = plt.imshow(mapb, cmap=colormap, origin="lower", alpha=1, extent=dim, vmin=vmin, vmax=vmax)
            # plt.title(titreb)
            # plt.xticks([])
            a = plt.gca()
            xax = a.axes.get_xaxis()
            xax = xax.set_visible(False)
            plt.grid(b=None)
            if i - init == 0:
                plt.ylabel("y [Mpc]", fontsize=fs)
                plt.yticks(fontsize=fs)
            else:
                plt.yticks([])
            if i == (init + n - 1):
                plt.text(1.05, 0.5, 'Compressive', transform=axs[1, n - 1].transAxes, verticalalignment='center',
                         rotation=270, fontsize=fs)

            # vmin = 0
            # vmax = 1
            # colormap = "binary"

            plt.sca(axs[2, i - init])
            im = plt.imshow(mapc, cmap=colormap, origin="lower", alpha=1, extent=dim, vmin=vmin, vmax=vmax)
            # plt.title(titrec)
            # plt.yticks([])
            plt.xlabel("z [Mpc]", fontsize=fs)
            plt.xticks(fontsize=fs)
            plt.grid(b=None)
            if i - init == 0:
                plt.ylabel("y [Mpc]", fontsize=fs)
                plt.yticks(fontsize=fs)
            else:
                plt.yticks([])
            if i == (init + n - 1):
                plt.text(1.05, 0.5, 'Solenoidal', transform=axs[2, n - 1].transAxes, verticalalignment='center',
                         rotation=270, fontsize=fs)

        cb = f.colorbar(im, ax=axs.ravel().tolist()) # ,size='small')#,shrink=0.6)
        cb.set_label(r'$||\vec{v}||~[\mathrm{km\,s^{-1}}]$', size=fs)  # , size='large'
        cb.ax.tick_params(labelsize=fs)

        # plt.subplots_adjust(wspace=0, hspace=0)

        # f.tight_layout()

        plt.show()

    def show_map_1x3_rectangle(vela, titrea, velb, titreb, velc, titrec, comp):
        slice = int(nx / 2)

        print("slice", slice)

        print('nx', nx)

        # slice = int(nx/2-nx/5)

        print("slice", slice)

        mapa = vela[slice, :, :]
        mapb = mapa
        mapc = mapa
        # mapb = velb[slice, :, :]
        # mapc = velc[slice, :, :]

        # mapa = vela[:, slice, :]
        # mapb = velb[:, slice, :]
        # mapc = velc[:, slice, :]

        # mapa = np.mean(vela[slice-1:slice, :, :],axis=0)

        # mapa = mapa.T
        # mapa = np.flip(mapa, axis=0)

        # mapb = mapb.T
        # mapb = np.flip(mapb, axis=0)

        # mapc = mapc.T
        # mapc = np.flip(mapc, axis=0)

        if comp == "comp":

            vmin = -1000
            vmax = 1000

            colormap = "bwr"

        elif comp == "tot":

            vmin = 0
            vmax = 1200

            colormap = "cividis"

        init = 0
        n = 8

        f, axs = plt.subplots(1, n, figsize=(12, 6), constrained_layout=True,facecolor='white')
        for i in range(init, init + n):
            print("i", i)

            fs = 14

            title = "$\Delta x_{cen}$=" + str(round(2 * i - 8, 3)) + " Mpc"

            slice = int(nx / 2 - (4 - i) * (nx / 10))

            # print("slice", slice)

            # print('nx', nx)

            # slice = int(nx/2-nx/5)

            print("slice", slice)

            mapa = vela[slice, :, :]
            # mapb = mapa
            # mapc = mapa

            mapb = velb[slice, :, :]
            mapc = velc[slice, :, :]

            plt.sca(axs[i - init])
            im = plt.imshow(mapa, cmap=colormap, origin="lower", alpha=1, extent=dim, vmin=vmin, vmax=vmax)
            plt.title(title, fontsize=fs)
            # plt.xlabel("x [Mpc]")
            if i - init == 0:
                plt.ylabel("y [Mpc]", fontsize=fs)
                plt.yticks(fontsize=fs)
            else:
                plt.yticks([])
            # plt.xticks([])
            #a = plt.gca()
            #xax = a.axes.get_xaxis()
            #xax = xax.set_visible(False)

            plt.grid(b=None)

            if i == (init + n - 1):
                plt.text(1.05, 0.5, 'Total', transform=axs[n - 1].transAxes, verticalalignment='center',
                         rotation=270, fontsize=fs)
            # plt.colorbar()

            # vmin = 0
            # vmax = 1
            # colormap = "binary"

            if i==0:
                anchor_z_5 = -1
                anchor_y_5 = 1
            elif i==1:
                anchor_z_5 = -1.5
                anchor_y_5 = -1
            elif i==2:
                anchor_z_5 = -2
                anchor_y_5 = -1
            elif i==3:
                anchor_z_5 = -2.5
                anchor_y_5 = -2.5
            elif i==4:
                anchor_z_5 = -2.5
                anchor_y_5 = -2.5
            elif i==5:
                anchor_z_5 = -2.5
                anchor_y_5 = -2.5
            elif i==6:
                anchor_z_5 = -2.5
                anchor_y_5 = -3.5
            elif i==7:
                anchor_z_5 = -3
                anchor_y_5 = -5


            length_5 = 5
            plt.gca().add_patch(Rectangle((anchor_z_5, anchor_y_5), length_5, length_5, edgecolor='white', facecolor='none', lw=1, ls="solid"))

            if i==0:
                anchor_z_2 = 0.5
                anchor_y_2 = 2.5
            elif i==1:
                anchor_z_2 = 0
                anchor_y_2 = 0.5
            elif i==2:
                anchor_z_2 = -0.5
                anchor_y_2 = 0.5
            elif i==3:
                anchor_z_2 = -1
                anchor_y_2 = -1
            elif i==4:
                anchor_z_2 = -1
                anchor_y_2 = -1
            elif i==5:
                anchor_z_2 = -1
                anchor_y_2 = -1
            elif i==6:
                anchor_z_2 = -1
                anchor_y_2 = -2
            elif i==7:
                anchor_z_2 = -1.5
                anchor_y_2 = -3.5



            length_2 = 2
            plt.gca().add_patch(Rectangle((anchor_z_2, anchor_y_2), length_2, length_2, edgecolor='white', facecolor='none', lw=1, ls="dashed"))

            if i==0:
                anchor_z_1 = 1
                anchor_y_1 = 3
                plt.axvline(x=1.5, color='white', linestyle='--', lw=1)
                plt.axhline(y=3.5, color='white', linestyle='--', lw=1)

            elif i==1:
                anchor_z_1 = 0.5
                anchor_y_1 = 1
                plt.axvline(x=1, color='white', linestyle='--', lw=1)
                plt.axhline(y=1.5, color='white', linestyle='--', lw=1)


            elif i==2:
                anchor_z_1 = 0
                anchor_y_1 = 1
                plt.axvline(x=0.5, color='white', linestyle='--', lw=1)
                plt.axhline(y=1.5, color='white', linestyle='--', lw=1)
            elif i==3:
                anchor_z_1 = -0.5
                anchor_y_1 = -0.5
                plt.axvline(x=0, color='white', linestyle='--', lw=1)
                plt.axhline(y=0, color='white', linestyle='--', lw=1)
            elif i==4:
                anchor_z_1 = -0.5
                anchor_y_1 = -0.5
                plt.axvline(x=0, color='white', linestyle='--', lw=1)
                plt.axhline(y=0, color='white', linestyle='--', lw=1)
            elif i==5:
                anchor_z_1 = -0.5
                anchor_y_1 = -0.5
                plt.axvline(x=0, color='white', linestyle='--', lw=1)
                plt.axhline(y=0, color='white', linestyle='--', lw=1)
            elif i==6:
                anchor_z_1 = -0.5
                anchor_y_1 = -1.5
                plt.axvline(x=0, color='white', linestyle='--', lw=1)
                plt.axhline(y=-1, color='white', linestyle='--', lw=1)
            elif i==7:
                anchor_z_1 = -1
                anchor_y_1 = -3
                plt.axvline(x=-0.5, color='white', linestyle='--', lw=1)
                plt.axhline(y=-2.5, color='white', linestyle='--', lw=1)


            length_1 = 1
            plt.gca().add_patch(Rectangle((anchor_z_1, anchor_y_1), length_1, length_1, edgecolor='white', facecolor='none', lw=1, ls="dotted"))

            plt.xlabel("z [Mpc]", fontsize=fs)
            plt.xticks(fontsize=fs)

            #plt.grid(b=None)

        ratiov =  1/8

        cb = f.colorbar(im, ax=axs.ravel().tolist(),fraction=0.047 * ratiov, pad=0.01)  # ,size='small')#,shrink=0.6)
        cb.set_label(r'$||\vec{v}||~[\mathrm{km\,s^{-1}}]$', size=fs)  # , size='large'
        cb.ax.tick_params(labelsize=fs)

        #plt.grid(b=None)

        plt.show()

        sys.exit()

    # show_map(vx,"total x",v_comp_x,"compressive x",v_sol_x,"solenoidal x","comp")

    # print("files saved")

    #show_map_nx3(v, "Total", v_comp, "Compressive", v_sol, "Solenoidal", "tot")

    #show_map_1x3_rectangle(v, "Total", v_comp, "Compressive", v_sol, "Solenoidal", "comp")

    show_map_nx3(v, "total x", v_comp, "compressive x", v_sol, "solenoidal x", "tot")

#show_decomposed_vel_fields()

def vproj_maps_hist():
    def gauss(x, A, mu, sigma):
        return A * norm.pdf(x, mu, sigma)
    def gaussian_fit(data, count, bins):

        # Bin centers for fitting
        bin_centers = (bins[:-1] + bins[1:]) / 2

        # Define Gaussian function

        model = GaussianModel()

        params = model.make_params(amplitude=1, center=np.mean(data), sigma=np.std(data))
        result = model.fit(count, params, x=bin_centers)

        print(result.fit_report())

        # print("best values",result.best_values)
        # print("best fit",result.best_fit)
        # print("red chi2",result.redchi)

        amplitude_best = result.params['amplitude'].value
        amplitude_error = result.params['amplitude'].stderr

        center_best = result.params['center'].value
        center_error = result.params['center'].stderr

        sigma_best = result.params['sigma'].value
        sigma_error = result.params['sigma'].stderr

        height_best = result.params['height'].value
        height_error = result.params['height'].stderr

        chi2 = result.chisqr
        red_chi2 = result.redchi

        popt = [amplitude_best, center_best, sigma_best, height_best]
        errors = [amplitude_error, center_error, sigma_error, height_error]

        # sys.exit()

        # Initial guesses: amplitude, mean, std
        # p0 = [0.015, np.mean(data), np.std(data)]

        # Perform the fit
        # popt, pcov = curve_fit(gauss, bin_centers, count, p0=p0)

        # Extract fitted parameters and their errors from the covariance matrix
        # A_fit, mu_fit, sigma_fit = popt
        # errors = np.sqrt(np.diag(pcov))

        # print(f'Fit results: A = {A_fit:.3g} ± {errors[0]:.3g}, '
        #      f'mean = {mu_fit:.3g} ± {errors[1]:.3g}, '
        #      f'std = {sigma_fit:.3g} ± {errors[2]:.3g}')

        # expected = gauss(bin_centers, *popt)
        # chi2 = np.sum((count - expected) ** 2 / expected)
        # dof = len(count) - len(popt)

        # print(f'Chi-square: {chi2:.3g}')
        # print(f'Reduced Chi-square: {chi2 / dof:.3g}')

        # red_chi2= chi2/dof

        # normalisation = A_fit/(np.sqrt(2*np.pi)*sigma_fit)
        # normalisation_error = (1/(np.sqrt(2*np.pi)*sigma_fit))*np.sqrt(errors[0]**2+((-A_fit*errors[2])/sigma_fit)**2)

        # print("normalisation",normalisation,"±",normalisation_error)

        return popt, errors, chi2, red_chi2, bin_centers

    def load_map(file,proj):
        h = FortranFile(file, 'r')
        nx, ny, nz = h.read_ints()
        cen_x, cen_y, cen_z = h.read_reals()
        print("nx,ny,nz", nx, ny, nz)
        # print("cen_x,cen_y,cen_z",cen_x,cen_y,cen_z)
        # sys.exit()
        # ncell = 0
        # print(ncell.dtype)
        # sys.exit()
        # ncell = np
        # sys.exit()
        if proj == "x":
            ncell = nz * ny
        elif proj == "y":
            ncell = nx * nz
        elif proj == "z":
            # ncell = np.uint64(nx * ny)
            ncell = nx * ny
        print("ncell", ncell)

        map = np.zeros(ncell)

        map = ftp.f90_to_py.read_map_file(ncell, file, 0)
        # map=h.read_reals()

        print("len nan",len(map[np.isnan(map)==True]),'ratio, ',len(map[np.isnan(map)==True])/len(map))
        #sys.exit()

        if proj == "x":
            map = np.reshape(map, (nz, ny))
            #map2 = np.reshape(map2, (nz, ny))
            #map3 = np.reshape(map3, (nz, ny))
            nl = nx
        elif proj == "y":
            map = np.reshape(map, (nx, nz))
            #map2 = np.reshape(map2, (nx, nz))
            #map3 = np.reshape(map3, (nx, nz))
            nl = ny
        elif proj == "z":
            map = np.reshape(map, (ny, nx))
            #map2 = np.reshape(map2, (ny, nx))
            #map3 = np.reshape(map3, (ny, nx))
            nl = nz

        print('file loaded')
        # sys.exit()

        # for i in range(len(map[0, :])):
        #    for j in range(len(map[:, 0])):
        #        if (map[i, j] == 0):
        #            map[i, j] = np.nan

        # map2 = np.loadtxt("map_2000px_z_ne_los_cleanall+T1e5_highr.txt")
        # for i in range(len(map2[0, :])):
        #    for j in range(len(map2[:, 0])):
        # if (map2[i, j] == 0):
        #    map2[i, j] = 1e-50

        # map-=map2
        print(map)
        print("min=", np.min(map))
        print("max=", np.max(map))

        return map

    map_vx = load_map('./maps/high_res/velocity/map_high_16f16_x_map_vx_ew_5Mpc2.bin','x')
    map_vy = load_map('./maps/high_res/velocity/map_high_16f16_y_map_vy_ew_5Mpc2.bin','y')
    map_vz = load_map('./maps/high_res/velocity/map_high_16f16_z_map_vz_ew_5Mpc2.bin','z')

    min_x = 0
    max_x = len(map_vx[:, 0])
    print("max", max_x)
    min_y = 0
    max_y = len(map_vx[0, :])

    cen_x = int(max_x / 2)
    cen_y = int(max_y / 2)

    print(cen_x, cen_y)

    #d = np.array([[np.sqrt((cen_x - (i + 0.5)) ** 2 + (cen_y - (j + 0.5)) ** 2) * 11.2524 for i in range(max_x)] for j in range(max_y)]) #still useful for toy model tests

    #np.save("map_distances_lvl16_proj_vel",d)

    d = np.load("map_distances_lvl16_proj_vel.npy")

    #print(np.shape(d))
    #print(np.min(d),np.max(d))
    #print(n)

    #plt.hist(d,bins=100)
    #plt.show()

    map_vx_rvir = map_vx[d<2147]
    map_vx_r500 = map_vx[d<1087]
    map_vx_500 = map_vx[d<500]
    map_vx_100 = map_vx[d<100]

    map_vy_rvir = map_vy[d<2147]
    map_vy_r500 = map_vy[d<1087]
    map_vy_500 = map_vy[d<500]
    map_vy_100 = map_vy[d<100]

    map_vz_rvir = map_vz[d<21147]
    map_vz_r500 = map_vz[d<1087]
    map_vz_500 = map_vz[d<500]
    map_vz_100 = map_vz[d<100]

    f, axs = plt.subplots(5, 4, figsize=(20, 8), constrained_layout=True)

    large = 1

    if large == 1:
        r1 = 2147
        r2 = 1087
        plt.suptitle("Velocity components distribution within $R_{vir}$(=2.15Mpc) and $R_{500}$(=1.1Mpc)")
        pdfrvir = "PDF \n ($R<R_{vir}$)"
        pdfr500 = "PDF \n ($R<R_{500}$)"
        height = 0.0017
        # height = 0.0035

    if large == 0:
        r1 = 500
        r2 = 100
        plt.suptitle("Velocity components distribution within spheres of 500kpc and 100kpc radii")
        pdfrvir = "PDF \n ($R<500 kpc$)"
        pdfr500 = "PDF \n ($R<100 kpc$)"
        height = 0.0035

    #vx_r1 = vx[r < r1]
    #vy_r1 = vy[r < r1]
    #vz_r1 = vz[r < r1]

    # print("nbr of cells with T<10^7K", len(vx_r1[t[r < r1] < 1e7]), "ratio",len(vx_r1[t[r < r1] < 1e7]) / len(vx_r1))

    # vx_r1 = vx_r1[t[r < r1] < 1e7]

    # meanvx,stdvx = norm.fit(vx_r1)
    # print(f"Vx fit:mean={meanvx:.4g},"f"std={stdvx:.4g}")

    print(
        f"map_vx_rvir stat: mean={np.mean(map_vx_rvir):.4g}, "f"median={np.median(map_vx_rvir):.4g}, "f"std={np.std(map_vx_rvir):.4g}, "f"skewness={skew(map_vx_rvir):.4g}, "f"kurtosis={kurtosis(map_vx_rvir):.4g}")

    # sys.exit()

    print(
        f"map_vy_rvir stat: mean={np.mean(map_vy_rvir):.4g}, "f"median={np.median(map_vy_rvir):.4g}, "f"std={np.std(map_vy_rvir):.4g}, "f"skewness={skew(map_vy_rvir):.4g}, "f"kurtosis={kurtosis(map_vy_rvir):.4g}")
    print(
        f"map_vz_rvir stat: mean={np.mean(map_vz_rvir):.4g}, "f"median={np.median(map_vz_rvir):.4g}, "f"std={np.std(map_vz_rvir):.4g}, "f"skewness={skew(map_vz_rvir):.4g}, "f"kurtosis={kurtosis(map_vz_rvir):.4g}")

    #map_vx_r500 = vx[r < r2]
    #map_vy_r500 = vy[r < r2]
    #map_vz_r500 = vz[r < r2]

    print(
        f"map_vx_r500 stat: mean={np.mean(map_vx_r500):.4g}, "f"median={np.median(map_vx_r500):.4g}, "f"std={np.std(map_vx_r500):.4g}, "f"skewness={skew(map_vx_r500):.4g}, "f"kurtosis={kurtosis(map_vx_r500):.4g}")
    print(
        f"map_vy_r500 stat: mean={np.mean(map_vy_r500):.4g}, "f"median={np.median(map_vy_r500):.4g}, "f"std={np.std(map_vy_r500):.4g}, "f"skewness={skew(map_vy_r500):.4g}, "f"kurtosis={kurtosis(map_vy_r500):.4g}")
    print(
        f"map_vz_r500 stat: mean={np.mean(map_vz_r500):.4g}, "f"median={np.median(map_vz_r500):.4g}, "f"std={np.std(map_vz_r500):.4g}, "f"skewness={skew(map_vz_r500):.4g}, "f"kurtosis={kurtosis(map_vz_r500):.4g}")

    count, bins, _ = axs[0, 0].hist(map_vx_rvir, bins=20, alpha=0.6, color='darkred', density=True)

    # print("count",count)
    # sys.exit()

    popt, errors, chi2, red_chi2, bin_cen = gaussian_fit(map_vx_rvir, count, bins)
    axs[0, 0].plot(bin_cen, gauss(bin_cen, *popt[0:3]), color='black',
                   label='best fit : \n $\mu$={:.3g}±{:.3g}, $\sigma$={:.3g}±{:.3g} \n $Norm$={:.3g}±{:.3g} \n $\chi^2$={:.3g}, $\chi^2/dof$={:.3g}'.format(
                       popt[1], errors[1], popt[2], errors[2], popt[3], errors[3], chi2, red_chi2))

    # axs[0,0].yscale('log')
    axs[0, 0].set_xlim(-1900, 1900)
    axs[0, 0].set_ylim(0, height)
    # axs[0, 0].set_ylim(0, 1.2e6)
    # axs[0,0].legend(fontsize=10)
    # axs[0, 0].set_title("$V_x$")
    # axs[0,0].set_xlabel("Vx [km/s]")
    # axs[0, 0].set_ylabel("PDF ($R<R_{vir}$)")
    axs[0, 0].set_ylabel(pdfrvir)
    axs[0, 0].set_xticks([])
    # axs[0, 0].legend(fontsize=10)

    count, bins, _ = axs[0, 1].hist(map_vy_rvir, bins=20, alpha=0.6, color='darkblue', density=True)

    popt, errors, chi2, red_chi2, bin_cen = gaussian_fit(map_vy_rvir, count, bins)
    axs[0, 1].plot(bin_cen, gauss(bin_cen, *popt[0:3]), color='black',
                   label='best fit : \n $\mu$={:.3g}±{:.3g}, $\sigma$={:.3g}±{:.3g} \n $Norm$={:.3g}±{:.3g} \n $\chi^2$={:.3g}, $\chi^2/dof$={:.3g}'.format(
                       popt[1], errors[1], popt[2], errors[2], popt[3], errors[3], chi2, red_chi2))

    # xs[0,1].yscale('log')
    axs[0, 1].set_xlim(-1900, 1900)
    axs[0, 1].set_ylim(0, height)
    # axs[0, 1].set_ylim(0, 1.2e6)
    # axs[0,1].legend(fontsize=10)
    # axs[0, 1].set_title("y component")
    axs[0, 1].set_yticks([])
    # axs[0,1].set_xlabel("Vy [km/s]")
    axs[0, 1].set_xticks([])
    # axs[0, 1].legend(fontsize=10)

    count, bins, _ = axs[0, 2].hist(map_vz_rvir, bins=20, alpha=0.6, color='darkgreen', density=True)

    popt, errors, chi2, red_chi2, bin_cen = gaussian_fit(map_vz_rvir, count, bins)
    axs[0, 2].plot(bin_cen, gauss(bin_cen, *popt[0:3]), color='black',
                   label='best fit : \n $\mu$={:.3g}±{:.3g}, $\sigma$={:.3g}±{:.3g} \n $Norm$={:.3g}±{:.3g} \n $\chi^2$={:.3g}, $\chi^2/dof$={:.3g}'.format(
                       popt[1], errors[1], popt[2], errors[2], popt[3], errors[3], chi2, red_chi2))

    # axs[0,2].yscale('log')
    axs[0, 2].set_xlim(-1900, 1900)
    axs[0, 2].set_ylim(0, height)
    # axs[0, 2].set_ylim(0, 1.2e6)
    # axs[0,2].legend(fontsize=10)
    # axs[0, 2].set_title("z projection")
    # axs[0,2].set_xlabel("Vz [km/s]")
    axs[0, 2].set_yticks([])
    axs[0, 2].set_xticks([])
    # axs[0, 2].legend(fontsize=10)

    axs[0, 3].hist(map_vx_rvir, bins=20, alpha=0.6, color='darkred', label='$V_x$', density=True)
    axs[0, 3].hist(map_vy_rvir, bins=20, alpha=0.6, color='darkblue', label='$V_y$', density=True)
    axs[0, 3].hist(map_vz_rvir, bins=20, alpha=0.6, color='darkgreen', label='$V_z$', density=True)
    axs[0, 3].set_xlim(-1900, 1900)
    axs[0, 3].set_ylim(0, height)
    # axs[0, 3].set_ylim(0, 1.2e6)
    axs[0, 3].set_yticks([])
    axs[0, 3].set_xticks([])
    axs[0, 3].legend(fontsize=10)

    # plt.text(1.05, 0.5, '$R<R_{vir}$', transform=axs[0, 3].transAxes, verticalalignment='center',
    #         rotation=270, fontsize=14)

    count, bins, _ = axs[1, 0].hist(map_vx_r500, bins=20, alpha=0.6, color='red', density=True)

    popt, errors, chi2, red_chi2, bin_cen = gaussian_fit(map_vx_r500, count, bins)
    axs[1, 0].plot(bin_cen, gauss(bin_cen, *popt[0:3]), color='black',
                   label='best fit : \n $\mu$={:.3g}±{:.3g}, $\sigma$={:.3g}±{:.3g} \n $Norm$={:.3g}±{:.3g} \n $\chi^2$={:.3g}, $\chi^2/dof$={:.3g}'.format(
                       popt[1], errors[1], popt[2], errors[2], popt[3], errors[3], chi2, red_chi2))

    # axs[1,0].yscale('log')
    axs[1, 0].set_xlim(-1900, 1900)
    axs[1, 0].set_ylim(0, height)
    # axs[1, 0].set_ylim(0, 1.3e4)
    # axs[1,0].legend(fontsize=10)
    # axs[1,0].set_title("x projection")
    # axs[1,0].set_xlabel("Vx [km/s]")
    # axs[1, 0].set_ylabel("PDF ($R<R_{500}$)")
    axs[1, 0].set_ylabel(pdfr500)
    axs[1, 0].set_xticks([])
    # axs[1, 0].legend(fontsize=10)

    count, bins, _ = axs[1, 1].hist(map_vy_r500, bins=20, alpha=0.6, color='blue', density=True)

    popt, errors, chi2, red_chi2, bin_cen = gaussian_fit(map_vy_r500, count, bins)
    axs[1, 1].plot(bin_cen, gauss(bin_cen, *popt[0:3]), color='black',
                   label='best fit : \n $\mu$={:.3g}±{:.3g}, $\sigma$={:.3g}±{:.3g} \n $Norm$={:.3g}±{:.3g} \n $\chi^2$={:.3g}, $\chi^2/dof$={:.3g}'.format(
                       popt[1], errors[1], popt[2], errors[2], popt[3], errors[3], chi2, red_chi2))

    # xs[1,1].yscale('log')
    axs[1, 1].set_xlim(-1900, 1900)
    axs[1, 1].set_ylim(0, height)
    # axs[1, 1].set_ylim(0, 1.3e4)
    # axs[1,1].legend(fontsize=10)
    # axs[1,1].set_title("y projection")
    axs[1, 1].set_yticks([])
    # axs[1,1].set_xlabel("Vy [km/s]")
    axs[1, 1].set_xticks([])
    # axs[1, 1].legend(fontsize=10)

    count, bins, _ = axs[1, 2].hist(map_vz_r500, bins=20, alpha=0.6, color='green', density=True)

    popt, errors, chi2, red_chi2, bin_cen = gaussian_fit(map_vz_r500, count, bins)
    axs[1, 2].plot(bin_cen, gauss(bin_cen, *popt[0:3]), color='black',
                   label='best fit : \n $\mu$={:.3g}±{:.3g}, $\sigma$={:.3g}±{:.3g} \n $Norm$={:.3g}±{:.3g} \n $\chi^2$={:.3g}, $\chi^2/dof$={:.3g}'.format(
                       popt[1], errors[1], popt[2], errors[2], popt[3], errors[3], chi2, red_chi2))

    # axs[1,2].yscale('log')
    axs[1, 2].set_xlim(-1900, 1900)
    axs[1, 2].set_ylim(0, height)
    # axs[1, 2].set_ylim(0, 1.3e4)
    # axs[1,2].legend(fontsize=10)
    # axs[1,2].set_title("z projection")
    # axs[1,2].set_xlabel("Vz [km/s]")
    axs[1, 2].set_yticks([])
    axs[1, 2].set_xticks([])
    # axs[1, 2].legend(fontsize=10)

    axs[1, 3].hist(map_vx_r500, bins=20, alpha=0.6, color='red', label='$V_x$', density=True)
    axs[1, 3].hist(map_vy_r500, bins=20, alpha=0.6, color='blue', label='$V_y$', density=True)
    axs[1, 3].hist(map_vz_r500, bins=20, alpha=0.6, color='green', label='$V_z$', density=True)
    axs[1, 3].set_xlim(-1900, 1900)
    axs[1, 3].set_ylim(0, height)
    # axs[1, 3].set_ylim(0, 1.3e4)
    axs[1, 3].set_yticks([])
    axs[1, 3].set_xticks([])
    axs[1, 3].legend(fontsize=10)
    # plt.text(1.05, 0.5, '$R<R_{500}$', transform=axs[1, 3].transAxes, verticalalignment='center',
    #         rotation=270, fontsize=14)

    ########################################
    ########################################

    large = 0

    if large == 0:
        r1 = 500
        r2 = 100
        plt.suptitle("Velocity components PDFs within spheres of 500kpc and 100kpc radii")
        pdfrvir = "PDF \n ($R<500 kpc$)"
        pdfr500 = "PDF \n ($R<100 kpc$)"
        height1 = 0.002
        height = 0.0035

    #vx_r1 = vx[r < r1]
    #vy_r1 = vy[r < r1]
    #vz_r1 = vz[r < r1]

    #vx_r2 = vx[r < r2]
    #vy_r2 = vy[r < r2]
    #vz_r2 = vz[r < r2]

    count, bins, _ = axs[2, 0].hist(map_vx_500, bins=20, alpha=0.6, color='orange', density=True)

    # print("count",count)
    # sys.exit()

    popt, errors, chi2, red_chi2, bin_cen = gaussian_fit(map_vx_500, count, bins)
    axs[2, 0].plot(bin_cen, gauss(bin_cen, *popt[0:3]), color='black',
                   label='best fit : \n $\mu$={:.3g}±{:.3g}, $\sigma$={:.3g}±{:.3g} \n $Norm$={:.3g}±{:.3g} \n $\chi^2$={:.3g}, $\chi^2/dof$={:.3g}'.format(
                       popt[1], errors[1], popt[2], errors[2], popt[3], errors[3], chi2, red_chi2))

    # axs[2,0].yscale('log')
    axs[2, 0].set_xlim(-1900, 1900)
    axs[2, 0].set_ylim(0, height1)
    # axs[2, 0].set_ylim(0, 1.2e6)
    # axs[2,0].legend(fontsize=10)
    # axs[2, 0].set_title("$V_x$")
    # axs[2,0].set_xlabel("Vx [km/s]")
    # axs[2, 0].set_ylabel("PDF ($R<R_{vir}$)")
    axs[2, 0].set_ylabel(pdfrvir)
    axs[2, 0].set_xticks([])
    # axs[2, 0].legend(fontsize=10)

    count, bins, _ = axs[2, 1].hist(map_vy_500, bins=20, alpha=0.6, color='purple', density=True)

    popt, errors, chi2, red_chi2, bin_cen = gaussian_fit(map_vy_500, count, bins)
    axs[2, 1].plot(bin_cen, gauss(bin_cen, *popt[0:3]), color='black',
                   label='best fit : \n $\mu$={:.3g}±{:.3g}, $\sigma$={:.3g}±{:.3g} \n $Norm$={:.3g}±{:.3g} \n $\chi^2$={:.3g}, $\chi^2/dof$={:.3g}'.format(
                       popt[1], errors[1], popt[2], errors[2], popt[3], errors[3], chi2, red_chi2))

    # xs[2,1].yscale('log')
    axs[2, 1].set_xlim(-1900, 1900)
    axs[2, 1].set_ylim(0, height1)
    # axs[2, 1].set_ylim(0, 1.2e6)
    # axs[2,1].legend(fontsize=10)
    # axs[2, 1].set_title("y component")
    axs[2, 1].set_yticks([])
    # axs[2,1].set_xlabel("Vy [km/s]")
    axs[2, 1].set_xticks([])
    # axs[2, 1].legend(fontsize=10)

    count, bins, _ = axs[2, 2].hist(map_vz_500, bins=20, alpha=0.6, color='turquoise', density=True)

    popt, errors, chi2, red_chi2, bin_cen = gaussian_fit(map_vz_500, count, bins)
    axs[2, 2].plot(bin_cen, gauss(bin_cen, *popt[0:3]), color='black',
                   label='best fit : \n $\mu$={:.3g}±{:.3g}, $\sigma$={:.3g}±{:.3g} \n $Norm$={:.3g}±{:.3g} \n $\chi^2$={:.3g}, $\chi^2/dof$={:.3g}'.format(
                       popt[1], errors[1], popt[2], errors[2], popt[3], errors[3], chi2, red_chi2))

    # axs[2,2].yscale('log')
    axs[2, 2].set_xlim(-1900, 1900)
    axs[2, 2].set_ylim(0, height1)
    # axs[2, 2].set_ylim(0, 1.2e6)
    # axs[2,2].legend(fontsize=10)
    # axs[2, 2].set_title("z projection")
    # axs[2,2].set_xlabel("Vz [km/s]")
    axs[2, 2].set_yticks([])
    axs[2, 2].set_xticks([])
    # axs[2, 2].legend(fontsize=10)

    axs[2, 3].hist(map_vx_500, bins=20, alpha=0.6, color='orange', label='$V_x$', density=True)
    axs[2, 3].hist(map_vy_500, bins=20, alpha=0.6, color='purple', label='$V_y$', density=True)
    axs[2, 3].hist(map_vz_500, bins=20, alpha=0.6, color='turquoise', label='$V_z$', density=True)
    axs[2, 3].set_xlim(-1900, 1900)
    axs[2, 3].set_ylim(0, height1)
    # axs[2, 3].set_ylim(0, 1.2e6)
    axs[2, 3].set_yticks([])
    axs[2, 3].set_xticks([])
    axs[2, 3].legend(fontsize=10)

    # plt.text(1.05, 0.5, '$R<R_{vir}$', transform=axs[2, 3].transAxes, verticalalignment='center',
    #         rotation=270, fontsize=14)

    count, bins, _ = axs[3, 0].hist(map_vx_100, bins=20, alpha=0.6, color='gold', density=True)

    popt, errors, chi2, red_chi2, bin_cen = gaussian_fit(map_vx_100, count, bins)
    axs[3, 0].plot(bin_cen, gauss(bin_cen, *popt[0:3]), color='black',
                   label='best fit : \n $\mu$={:.3g}±{:.3g}, $\sigma$={:.3g}±{:.3g} \n $Norm$={:.3g}±{:.3g} \n $\chi^2$={:.3g}, $\chi^2/dof$={:.3g}'.format(
                       popt[1], errors[1], popt[2], errors[2], popt[3], errors[3], chi2, red_chi2))

    # axs[3,0].yscale('log')
    axs[3, 0].set_xlim(-1900, 1900)
    axs[3, 0].set_ylim(0, height)
    # axs[3, 0].set_ylim(0, 1.3e4)
    # axs[3,0].legend(fontsize=10)
    # axs[3,0].set_title("x projection")
    # axs[3,0].set_xlabel("Vx [km/s]")
    # axs[3, 0].set_ylabel("PDF ($R<R_{500}$)")
    axs[3, 0].set_ylabel(pdfr500)
    axs[3, 0].set_xticks([])
    # axs[3, 0].legend(fontsize=10)

    count, bins, _ = axs[3, 1].hist(map_vy_100, bins=20, alpha=0.6, color='pink', density=True)

    popt, errors, chi2, red_chi2, bin_cen = gaussian_fit(map_vy_100, count, bins)
    axs[3, 1].plot(bin_cen, gauss(bin_cen, *popt[0:3]), color='black',
                   label='best fit : \n $\mu$={:.3g}±{:.3g}, $\sigma$={:.3g}±{:.3g} \n $Norm$={:.3g}±{:.3g} \n $\chi^2$={:.3g}, $\chi^2/dof$={:.3g}'.format(
                       popt[1], errors[1], popt[2], errors[2], popt[3], errors[3], chi2, red_chi2))

    # xs[3,1].yscale('log')
    axs[3, 1].set_xlim(-1900, 1900)
    axs[3, 1].set_ylim(0, height)
    # axs[3, 1].set_ylim(0, 1.3e4)
    # axs[3,1].legend(fontsize=10)
    # axs[3,1].set_title("y projection")
    axs[3, 1].set_yticks([])
    # axs[3,1].set_xlabel("Vy [km/s]")
    axs[3, 1].set_xticks([])
    # axs[3, 1].legend(fontsize=10)

    count, bins, _ = axs[3, 2].hist(map_vz_100, bins=20, alpha=0.6, color='cyan', density=True)

    popt, errors, chi2, red_chi2, bin_cen = gaussian_fit(map_vz_100, count, bins)
    axs[3, 2].plot(bin_cen, gauss(bin_cen, *popt[0:3]), color='black',
                   label='best fit : \n $\mu$={:.3g}±{:.3g}, $\sigma$={:.3g}±{:.3g} \n $Norm$={:.3g}±{:.3g} \n $\chi^2$={:.3g}, $\chi^2/dof$={:.3g}'.format(
                       popt[1], errors[1], popt[2], errors[2], popt[3], errors[3], chi2, red_chi2))

    # axs[3,2].yscale('log')
    axs[3, 2].set_xlim(-1900, 1900)
    axs[3, 2].set_ylim(0, height)
    # axs[3, 2].set_ylim(0, 1.3e4)
    # axs[3,2].legend(fontsize=10)
    # axs[3,2].set_title("z projection")
    # axs[3,2].set_xlabel("Vz [km/s]")
    axs[3, 2].set_yticks([])
    axs[3, 2].set_xticks([])
    # axs[3, 2].legend(fontsize=10)

    axs[3, 3].hist(map_vx_100, bins=20, alpha=0.6, color='gold', label='$V_x$', density=True)
    axs[3, 3].hist(map_vy_100, bins=20, alpha=0.6, color='pink', label='$V_y$', density=True)
    axs[3, 3].hist(map_vz_100, bins=20, alpha=0.6, color='cyan', label='$V_z$', density=True)
    axs[3, 3].set_xlim(-1900, 1900)
    axs[3, 3].set_ylim(0, height)
    # axs[3, 3].set_ylim(0, 1.3e4)
    axs[3, 3].set_yticks([])
    axs[3, 3].legend(fontsize=10)
    # plt.text(1.05, 0.5, '$R<R_{500}$', transform=axs[3, 3].transAxes, verticalalignment='center',
    #         rotation=270, fontsize=14)

    ########################################
    ########################################

    rvir = 2147
    r500 = 1087

    axs[4, 0].hist(map_vx_rvir, bins=20, alpha=0.6, color='darkred', label='$R<R_{vir}$',
                   density=True)  # ,histtype='step',linewidth=2)
    axs[4, 0].hist(map_vx_r500, bins=20, alpha=0.6, color='red', label='$R<R_{500}$',
                   density=True)  # ,histtype='step',linewidth=2)
    axs[4, 0].hist(map_vx_500, bins=20, alpha=0.6, color='orange', label='$R<500kpc$',
                   density=True)  # ,histtype='step',linewidth=2)
    axs[4, 0].hist(map_vx_100, bins=20, alpha=0.6, color='gold', label='$R<100kpc$',
                   density=True)  # ,histtype='step',linewidth=2)
    # axs[4,0].yscale('log')
    axs[4, 0].set_xlim(-1900, 1900)
    axs[4, 0].set_ylim(0, height)
    # axs[4, 0].set_ylim(0, 700)
    # axs[4,0].legend(fontsize=10)
    # axs[4,0].set_title("x projection")
    axs[4, 0].set_xlabel("Vx [km/s]")
    axs[4, 0].set_ylabel("PDF")
    axs[4, 0].legend(fontsize=10)
    # axs[4, 0].set_yscale('log')

    axs[4, 1].hist(map_vy_rvir, bins=20, alpha=0.6, color='darkblue', label='$R<R_{vir}$',
                   density=True)  # ,histtype='step',linewidth=2)
    axs[4, 1].hist(map_vy_r500, bins=20, alpha=0.6, color='blue', label='$R<R_{500}$',
                   density=True)  # ,histtype='step',linewidth=2)
    axs[4, 1].hist(map_vy_500, bins=20, alpha=0.6, color='purple', label='$R<500kpc$',
                   density=True)  # ,histtype='step',linewidth=2)
    axs[4, 1].hist(map_vy_100, bins=20, alpha=0.6, color='pink', label='$R<100kpc$',
                   density=True)  # ,histtype='step',linewidth=2)
    # xs[4,1].yscale('log')
    axs[4, 1].set_xlim(-1900, 1900)
    axs[4, 1].set_ylim(0, height)
    # axs[4, 1].set_ylim(0, 700)
    # axs[4,1].legend(fontsize=10)
    # axs[4,1].set_title("y projection")
    axs[4, 1].set_yticks([])
    axs[4, 1].set_xlabel("Vy [km/s]")
    axs[4, 1].legend(fontsize=10)
    # axs[4, 1].set_yscale('log')

    axs[4, 2].hist(map_vz_rvir, bins=20, alpha=0.6, color='darkgreen', label='$R<R_{vir}$',
                   density=True)  # ,histtype='step',linewidth=2)
    axs[4, 2].hist(map_vz_r500, bins=20, alpha=0.6, color='green', label='$R<R_{500}$',
                   density=True)  # ,histtype='step',linewidth=2)
    axs[4, 2].hist(map_vz_500, bins=20, alpha=0.6, color='turquoise', label='$R<500kpc$',
                   density=True)  # ,histtype='step',linewidth=2)
    axs[4, 2].hist(map_vz_100, bins=20, alpha=0.6, color='cyan', label='$R<100kpc$',
                   density=True)  # ,histtype='step',linewidth=2)
    # axs[4,2].yscale('log')
    axs[4, 2].set_xlim(-1900, 1900)
    axs[4, 2].set_ylim(0, height)
    # axs[4, 2].set_ylim(0, 700)
    # axs[4,2].legend(fontsize=10)
    # axs[4,2].set_title("z projection")
    axs[4, 2].set_xlabel("Vz [km/s]")
    axs[4, 2].set_yticks([])
    axs[4, 2].legend(fontsize=10)
    # axs[4, 2].set_yscale('log')

    # axs[3].yscale('log')

    # plt.hist(v_Trange_list_100kpc, bins=20, alpha=0.6, color='green', label='100kpc**2, XRISM T range')
    # plt.grid(b=None)
    plt.legend()
    # plt.ylabel("PDF")
    axs[3, 3].set_xlabel("V [km/s]")

    axs[4, 3].axis("off")

    # plt.yscale('log')
    # plt.title("Vz distribution on ew sightline velocity along z axis")
    # plt.title("Velocity distribution along sightlines, 100kpc**2 maps, XRISM T range")

    plt.suptitle("Velocity components distribution within $R_{vir}$, $R_{500}$, 500kpc and 100kpc radii spheres")

    plt.show()

    sys.exit()


    sys.exit()

#vproj_maps_hist()

def vproj_maps_1x4_plots():

    def show_radii():
        xpx = np.linspace(-10,10,4000)
        ypx = np.linspace(-10,10,4000)

        Xpx, Ypx = np.meshgrid(xpx, ypx)

        xpx = np.linspace(-2.5, 2.5, 4000)
        ypx = np.linspace(-2.5, 2.5, 4000)

        Xpx, Ypx = np.meshgrid(xpx, ypx)

        F = Xpx ** 2 + Ypx ** 2 - (1087 * (10 / 1e4)) ** 2  # r500
        r500 = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=2, alpha=0.9, linestyles='dashed')
        pos = [(-0.75, 2)]
        plt.clabel(r500, r500.levels, inline=True, fmt="$R_{500}$", fontsize=20, manual=pos)
        F = Xpx ** 2 + Ypx ** 2 - (2147 * (10 / 1e4)) ** 2  # rvir
        rvir = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=2, alpha=0.9, linestyles='solid')
        plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{vir}$", fontsize=20, manual=pos)
        F = Xpx ** 2 + Ypx ** 2 - (515 * (10 / 1e4)) ** 2  # r2500
        r2500 = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=2, alpha=0.9, linestyles='dotted')
        plt.clabel(r2500, r2500.levels, inline=True, fmt="$R_{2500}$", fontsize=20, manual=pos)

        #F = Xpx ** 2 + Ypx ** 2 - (8000 * (10 / 1e4)) ** 2  # background lower limit
        #rbgl = plt.contour(Xpx, Ypx, F, [0], colors='grey', linewidths=0.6, alpha=0.9)
        #plt.clabel(rbgl, rbgl.levels, inline=True, fmt="$Background lower limit}$", fontsize=10)

        #F = Xpx ** 2 + Ypx ** 2 - (10000 * (10 / 1e4)) ** 2  # background upper limit
        #rbgh = plt.contour(Xpx, Ypx, F, [0], colors='grey', linewidths=0.6, alpha=0.9)
        #plt.clabel(rbgh, rbgh.levels, inline=True, fmt="$Background upper limit}$", fontsize=10)

        #plt.fill(rbgl,rbgh)

        #F = Xpx ** 2 + Ypx ** 2 - (1648 * (10 / 1e4)) ** 2 #r200
        #r200=plt.contour(Xpx, Ypx, F, [0], colors='blue', linewidths=0.6, alpha=0.9)
        #plt.clabel(r200, r200.levels, inline=True, fmt="$R_{200}$", fontsize=10)
        #F = Xpx ** 2 + Ypx ** 2 - (1085 * (10 / 1e4)) ** 2 #r500
        #r500=plt.contour(Xpx, Ypx, F, [0], colors='blue',linestyles="dashed", linewidths=0.6, alpha=0.9)
        #plt.clabel(r500, r500.levels, inline=True, fmt="$R_{500}$", fontsize=10)
        #F = Xpx ** 2 + Ypx ** 2 - (900 * (10 / 1e4)) ** 2  # bump radius
        #rbump = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
        #plt.clabel(rbump, rbump.levels, inline=True, fmt="$R_{bump}$", fontsize=10)
        #F = Xpx ** 2 + Ypx ** 2 - (850 * (10 / 1e4)) ** 2  # bump radius
        #rbump = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
        #plt.clabel(rbump, rbump.levels, inline=True, fmt="$R_{bump}$", fontsize=10)

        #F = Xpx ** 2 + Ypx ** 2 - (400 * (10 / 1e4)) ** 2  # low mw deproj radius
        #rdec = plt.contour(Xpx, Ypx, F, [0], colors='black', linestyles="dashed", linewidths=0.6, alpha=0.9)
        #plt.clabel(rdec, rdec.levels, inline=True, fmt="$R_{decrease in deproj prof}$", fontsize=10)
        #plt.scatter(x[i], y[i], alpha=0.9, s=1, c='white')
        #plt.show()
        #sys.exit()

    def show_quadrants():
        # xpx = np.linspace(-2.5, 2.5, 4000)
        # ypx = np.linspace(-2.5, 2.5, 4000)

        xpx = np.linspace(-2.5, 2.5, 4000)
        ypx = np.linspace(-2.5, 2.5, 4000)

        rvir = 2.147
        xrvir = rvir * 0.5 * np.sqrt(2)
        xpx_vir = np.linspace(-xrvir, xrvir, 4000)

        # Xpx, Ypx = np.meshgrid(xpx, ypx)

        # F = Xpx**2 + Ypx**2 - (1087*(10/1e4))**2  #r500
        # r500=plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=2, alpha=0.9,linestyles='dashed')
        # plt.clabel(r500,r500.levels,inline=True,fmt="$R_{500}$",fontsize=20)

        # F = Xpx ** 2 + Ypx ** 2 - (2147 * (10 / 1e4)) ** 2  # rvir
        # rvir = plt.contour(Xpx, Ypx, F, [0], colors='black', linewidths=2, alpha=0.9,linestyles='dashdot')
        # plt.clabel(rvir, rvir.levels, inline=True, fmt="$R_{vir}$", fontsize=20)

        plt.plot(xpx_vir, xpx_vir, color='black', lw=2, alpha=0.9)
        plt.plot(xpx_vir, -xpx_vir, color='black', lw=2, alpha=0.9)
        plt.axvline(x=0, ymin=(2.5 - rvir) / 5, ymax=(2.5 + rvir) / 5, color='black', lw=2, alpha=0.9)
        plt.axhline(y=0, xmin=(2.5 - rvir) / 5, xmax=(2.5 + rvir) / 5, color='black', lw=2, alpha=0.9)


    def load_map(file, proj):
        h = FortranFile(file, 'r')
        nx, ny, nz = h.read_ints()
        cen_x, cen_y, cen_z = h.read_reals()
        print("nx,ny,nz", nx, ny, nz)

        if proj == "x":
            ncell = nz * ny
        elif proj == "y":
            ncell = nx * nz
        elif proj == "z":
            ncell = nx * ny
        print("ncell", ncell)

        map = np.zeros(ncell)

        map = ftp.f90_to_py.read_map_file(ncell, file, 0)

        if proj == "x":
            map = np.reshape(map, (nz, ny))
            nl = nx
        elif proj == "y":
            map = np.reshape(map, (nx, nz))
            nl = ny
        elif proj == "z":
            map = np.reshape(map, (ny, nx))
            nl = nz

        return map,proj,nx,ny,nz

    def EW_vlos_1x4_plot():

        vx,proj,nx,ny,nz = load_map('./maps/high_res/velocity/15f15_analysis/map_high_15f15_x_map_vx_ew_Tsup7_5Mpc2.bin','x')
        vy,proj,nx,ny,nz = load_map('./maps/high_res/velocity/15f15_analysis/map_high_15f15_y_map_vy_ew_Tsup7_5Mpc2.bin','y')
        vz,proj,nx,ny,nz = load_map('./maps/high_res/velocity/15f15_analysis/map_high_15f15_z_map_vz_ew_Tsup7_5Mpc2.bin','z')
        vcen,proj,nx,ny,nz = load_map('./maps/high_res/velocity/15f15_analysis/map_high_15f15_cen_map_vz_ew_Tsup7_5Mpc2.bin','z')

        size = 22
        lvl = 15

        if proj == "x":
            dimx = (ny / 2) * (737.441 / 2 ** lvl)
            dimy = (nz / 2) * (737.441 / 2 ** lvl)
        elif proj == "y":
            dimx = (nz / 2) * (737.441 / 2 ** lvl)
            dimy = (nx / 2) * (737.441 / 2 ** lvl)
        elif proj == "z":
            dimx = (nx / 2) * (737.441 / 2 ** lvl)
            dimy = (ny / 2) * (737.441 / 2 ** lvl)


        dim = [-dimx, dimx, -dimy, dimy]
        print("dim", dim)

        fig, ax = plt.subplots(1, 4, figsize=(20, 5), constrained_layout=True,facecolor='white')#,gridspec_kw = {'wspace':0, 'hspace':0})

        plt.suptitle("Emission-weighted",size=size)

        plt.sca(ax[0])
        im = plt.imshow(vx, cmap="BR", origin='lower', alpha=1, extent=dim, vmin=-600, vmax=600)  # vmin=-2000, vmax=2000
        show_radii()
        #show_quadrants()
        #cb = fig.colorbar(im, ax=ax)

        #cb.ax.tick_params(labelsize=size)
        plt.xlabel("Mpc", size=size)
        plt.ylabel("Mpc", size=size)
        plt.yticks(fontsize=size)
        plt.xticks(fontsize=size)
        plt.grid(b=None)
        plt.title("$v_x$",size=size)

        plt.sca(ax[1])
        im = plt.imshow(vy, cmap="BR", origin='lower', alpha=1, extent=dim, vmin=-600, vmax=600)  # vmin=-2000, vmax=2000
        show_radii()
        #show_quadrants()
        # cb = fig.colorbar(im, ax=ax)

        # cb.ax.tick_params(labelsize=size)
        plt.xlabel("Mpc", size=size)
        #ax.set_ylabel("Mpc", size=size)
        plt.yticks([])
        plt.xticks(fontsize=size)
        plt.grid(b=None)
        plt.title("$v_y$", size=size)

        plt.sca(ax[2])
        im = plt.imshow(vz, cmap="BR", origin='lower', alpha=1, extent=dim, vmin=-600, vmax=600)  # vmin=-2000, vmax=2000
        show_radii()
        # show_quadrants()
        # cb = fig.colorbar(im, ax=ax)

        # cb.ax.tick_params(labelsize=size)
        plt.xlabel("Mpc", size=size)
        # ax.set_ylabel("Mpc", size=size)
        plt.yticks([])
        plt.xticks(fontsize=size)
        plt.grid(b=None)
        plt.title("$v_z$", size=size)

        plt.sca(ax[3])
        im = plt.imshow(vcen, cmap="BR", origin='lower', alpha=1, extent=dim, vmin=-600, vmax=600)  # vmin=-2000, vmax=2000
        show_radii()
        # show_quadrants()
        # cb = fig.colorbar(im, ax=ax)

        # cb.ax.tick_params(labelsize=size)
        plt.xlabel("Mpc", size=size)
        # ax.set_ylabel("Mpc", size=size)
        plt.yticks([])
        plt.xticks(fontsize=size)
        plt.grid(b=None)
        plt.title("$v_{cen}$", size=size)

        plt.show()
        
    def MW_vlos_1x4_plot():

        vx,proj,nx,ny,nz = load_map('./maps/high_res/velocity/15f15_analysis/map_high_15f15_x_map_vx_mw_Tsup7_5Mpc2.bin','x')
        vy,proj,nx,ny,nz = load_map('./maps/high_res/velocity/15f15_analysis/map_high_15f15_y_map_vy_mw_Tsup7_5Mpc2.bin','y')
        vz,proj,nx,ny,nz = load_map('./maps/high_res/velocity/15f15_analysis/map_high_15f15_z_map_vz_mw_Tsup7_5Mpc2.bin','z')
        vcen,proj,nx,ny,nz = load_map('./maps/high_res/velocity/15f15_analysis/map_high_15f15_cen_map_vz_mw_Tsup7_5Mpc2.bin','z')

        size = 22
        lvl = 15

        if proj == "x":
            dimx = (ny / 2) * (737.441 / 2 ** lvl)
            dimy = (nz / 2) * (737.441 / 2 ** lvl)
        elif proj == "y":
            dimx = (nz / 2) * (737.441 / 2 ** lvl)
            dimy = (nx / 2) * (737.441 / 2 ** lvl)
        elif proj == "z":
            dimx = (nx / 2) * (737.441 / 2 ** lvl)
            dimy = (ny / 2) * (737.441 / 2 ** lvl)


        dim = [-dimx, dimx, -dimy, dimy]
        print("dim", dim)

        fig, ax = plt.subplots(1, 4, figsize=(20, 5), constrained_layout=True,facecolor='white')#,gridspec_kw = {'wspace':0, 'hspace':0})

        plt.suptitle("Mass-weighted",size=size)

        plt.sca(ax[0])
        im = plt.imshow(vx, cmap="BR", origin='lower', alpha=1, extent=dim, vmin=-600, vmax=600)  # vmin=-2000, vmax=2000
        show_radii()
        #show_quadrants()
        #cb = fig.colorbar(im, ax=ax)

        #cb.ax.tick_params(labelsize=size)
        plt.xlabel("Mpc", size=size)
        plt.ylabel("Mpc", size=size)
        plt.yticks(fontsize=size)
        plt.xticks(fontsize=size)
        plt.grid(b=None)
        plt.title("$v_x$",size=size)

        plt.sca(ax[1])
        im = plt.imshow(vy, cmap="BR", origin='lower', alpha=1, extent=dim, vmin=-600, vmax=600)  # vmin=-2000, vmax=2000
        show_radii()
        #show_quadrants()
        # cb = fig.colorbar(im, ax=ax)

        # cb.ax.tick_params(labelsize=size)
        plt.xlabel("Mpc", size=size)
        #ax.set_ylabel("Mpc", size=size)
        plt.yticks([])
        plt.xticks(fontsize=size)
        plt.grid(b=None)
        plt.title("$v_y$", size=size)

        plt.sca(ax[2])
        im = plt.imshow(vz, cmap="BR", origin='lower', alpha=1, extent=dim, vmin=-600, vmax=600)  # vmin=-2000, vmax=2000
        show_radii()
        # show_quadrants()
        # cb = fig.colorbar(im, ax=ax)

        # cb.ax.tick_params(labelsize=size)
        plt.xlabel("Mpc", size=size)
        # ax.set_ylabel("Mpc", size=size)
        plt.yticks([])
        plt.xticks(fontsize=size)
        plt.grid(b=None)
        plt.title("$v_z$", size=size)

        plt.sca(ax[3])
        im = plt.imshow(vcen, cmap="BR", origin='lower', alpha=1, extent=dim, vmin=-600, vmax=600)  # vmin=-2000, vmax=2000
        show_radii()
        # show_quadrants()
        # cb = fig.colorbar(im, ax=ax)

        # cb.ax.tick_params(labelsize=size)
        plt.xlabel("Mpc", size=size)
        # ax.set_ylabel("Mpc", size=size)
        plt.yticks([])
        plt.xticks(fontsize=size)
        plt.grid(b=None)
        plt.title("$v_{cen}$", size=size)

        plt.show()

    def EW_vdisp_1x4_plot():

        vdx,proj,nx,ny,nz = load_map('./maps/high_res/velocity/15f15_analysis/map_high_15f15_x_map_vdx_ew_Tsup7_5Mpc2.bin','x')
        vdy,proj,nx,ny,nz = load_map('./maps/high_res/velocity/15f15_analysis/map_high_15f15_y_map_vdy_ew_Tsup7_5Mpc2.bin','y')
        vdz,proj,nx,ny,nz = load_map('./maps/high_res/velocity/15f15_analysis/map_high_15f15_z_map_vdz_ew_Tsup7_5Mpc2.bin','z')
        vdcen,proj,nx,ny,nz = load_map('./maps/high_res/velocity/15f15_analysis/map_high_15f15_cen_map_vdz_ew_Tsup7_5Mpc2.bin','z')

        size = 22
        lvl = 15

        if proj == "x":
            dimx = (ny / 2) * (737.441 / 2 ** lvl)
            dimy = (nz / 2) * (737.441 / 2 ** lvl)
        elif proj == "y":
            dimx = (nz / 2) * (737.441 / 2 ** lvl)
            dimy = (nx / 2) * (737.441 / 2 ** lvl)
        elif proj == "z":
            dimx = (nx / 2) * (737.441 / 2 ** lvl)
            dimy = (ny / 2) * (737.441 / 2 ** lvl)


        dim = [-dimx, dimx, -dimy, dimy]
        print("dim", dim)

        fig, ax = plt.subplots(1, 4, figsize=(20, 5), constrained_layout=True,facecolor='white')#,gridspec_kw = {'wspace':0, 'hspace':0})

        plt.suptitle("Emission-weighted",size=size)

        plt.sca(ax[0])
        im = plt.imshow(vdx, cmap="Sunset_r", origin='lower', alpha=1, extent=dim, vmin=0, vmax=800)
        show_radii()
        #show_quadrants()
        #cb = fig.colorbar(im, ax=ax)

        #cb.ax.tick_params(labelsize=size)
        plt.xlabel("Mpc", size=size)
        plt.ylabel("Mpc", size=size)
        plt.yticks(fontsize=size)
        plt.xticks(fontsize=size)
        plt.grid(b=None)
        plt.title("$\sigma_x$",size=size)

        plt.sca(ax[1])
        im = plt.imshow(vdy, cmap="Sunset_r", origin='lower', alpha=1, extent=dim, vmin=0, vmax=800)
        show_radii()
        #show_quadrants()
        # cb = fig.colorbar(im, ax=ax)

        # cb.ax.tick_params(labelsize=size)
        plt.xlabel("Mpc", size=size)
        #ax.set_ylabel("Mpc", size=size)
        plt.yticks([])
        plt.xticks(fontsize=size)
        plt.grid(b=None)
        plt.title("$\sigma_y$", size=size)

        plt.sca(ax[2])
        im = plt.imshow(vdz, cmap="Sunset_r", origin='lower', alpha=1, extent=dim, vmin=0, vmax=800)
        show_radii()
        # show_quadrants()
        # cb = fig.colorbar(im, ax=ax)

        # cb.ax.tick_params(labelsize=size)
        plt.xlabel("Mpc", size=size)
        # ax.set_ylabel("Mpc", size=size)
        plt.yticks([])
        plt.xticks(fontsize=size)
        plt.grid(b=None)
        plt.title("$\sigma_z$", size=size)

        plt.sca(ax[3])
        im = plt.imshow(vdcen, cmap="Sunset_r", origin='lower', alpha=1, extent=dim, vmin=0, vmax=800)
        show_radii()
        # show_quadrants()
        # cb = fig.colorbar(im, ax=ax)

        # cb.ax.tick_params(labelsize=size)
        plt.xlabel("Mpc", size=size)
        # ax.set_ylabel("Mpc", size=size)
        plt.yticks([])
        plt.xticks(fontsize=size)
        plt.grid(b=None)
        plt.title("$\sigma_{cen}$", size=size)

        plt.show()
        
    def MW_vdisp_1x4_plot():

        vdx,proj,nx,ny,nz = load_map('./maps/high_res/velocity/15f15_analysis/map_high_15f15_x_map_vdx_mw_Tsup7_5Mpc2.bin','x')
        vdy,proj,nx,ny,nz = load_map('./maps/high_res/velocity/15f15_analysis/map_high_15f15_y_map_vdy_mw_Tsup7_5Mpc2.bin','y')
        vdz,proj,nx,ny,nz = load_map('./maps/high_res/velocity/15f15_analysis/map_high_15f15_z_map_vdz_mw_Tsup7_5Mpc2.bin','z')
        vdcen,proj,nx,ny,nz = load_map('./maps/high_res/velocity/15f15_analysis/map_high_15f15_cen_map_vdz_mw_Tsup7_5Mpc2.bin','z')
        
        #vdx,proj,nx,ny,nz = load_map('./maps/high_res/velocity/16f16_analysis/map_high_16f16_x_map_vdx_mw_Tsup7_5Mpc2.bin','x')
        #vdy,proj,nx,ny,nz = load_map('./maps/high_res/velocity/16f16_analysis/map_high_16f16_y_map_vdy_mw_Tsup7_5Mpc2.bin','y')
        #vdz,proj,nx,ny,nz = load_map('./maps/high_res/velocity/16f16_analysis/map_high_16f16_z_map_vdz_mw_Tsup7_5Mpc2.bin','z')
        #vdcen,proj,nx,ny,nz = load_map('./maps/high_res/velocity/16f16_analysis/map_high_16f16_cen_map_vdz_mw_Tsup7_5Mpc2.bin','z')

        size = 22
        lvl = 15

        if proj == "x":
            dimx = (ny / 2) * (737.441 / 2 ** lvl)
            dimy = (nz / 2) * (737.441 / 2 ** lvl)
        elif proj == "y":
            dimx = (nz / 2) * (737.441 / 2 ** lvl)
            dimy = (nx / 2) * (737.441 / 2 ** lvl)
        elif proj == "z":
            dimx = (nx / 2) * (737.441 / 2 ** lvl)
            dimy = (ny / 2) * (737.441 / 2 ** lvl)


        dim = [-dimx, dimx, -dimy, dimy]
        print("dim", dim)

        fig, ax = plt.subplots(1, 4, figsize=(20, 5), constrained_layout=True,facecolor='white')#,gridspec_kw = {'wspace':0, 'hspace':0})

        plt.suptitle("Mass-weighted",size=size)

        plt.sca(ax[0])
        im = plt.imshow(vdx, cmap="Sunset_r", origin='lower', alpha=1, extent=dim, vmin=0, vmax=800)
        show_radii()
        #show_quadrants()
        #cb = fig.colorbar(im, ax=ax)

        #cb.ax.tick_params(labelsize=size)
        plt.xlabel("Mpc", size=size)
        plt.ylabel("Mpc", size=size)
        plt.yticks(fontsize=size)
        plt.xticks(fontsize=size)
        plt.grid(b=None)
        plt.title("$\sigma_x$",size=size)

        plt.sca(ax[1])
        im = plt.imshow(vdy, cmap="Sunset_r", origin='lower', alpha=1, extent=dim, vmin=0, vmax=800)
        show_radii()
        #show_quadrants()
        # cb = fig.colorbar(im, ax=ax)

        # cb.ax.tick_params(labelsize=size)
        plt.xlabel("Mpc", size=size)
        #ax.set_ylabel("Mpc", size=size)
        plt.yticks([])
        plt.xticks(fontsize=size)
        plt.grid(b=None)
        plt.title("$\sigma_y$", size=size)

        plt.sca(ax[2])
        im = plt.imshow(vdz, cmap="Sunset_r", origin='lower', alpha=1, extent=dim, vmin=0, vmax=800)
        show_radii()
        # show_quadrants()
        # cb = fig.colorbar(im, ax=ax)

        # cb.ax.tick_params(labelsize=size)
        plt.xlabel("Mpc", size=size)
        # ax.set_ylabel("Mpc", size=size)
        plt.yticks([])
        plt.xticks(fontsize=size)
        plt.grid(b=None)
        plt.title("$\sigma_z$", size=size)

        plt.sca(ax[3])
        im = plt.imshow(vdcen, cmap="Sunset_r", origin='lower', alpha=1, extent=dim, vmin=0, vmax=800)
        show_radii()
        # show_quadrants()
        # cb = fig.colorbar(im, ax=ax)

        # cb.ax.tick_params(labelsize=size)
        plt.xlabel("Mpc", size=size)
        # ax.set_ylabel("Mpc", size=size)
        plt.yticks([])
        plt.xticks(fontsize=size)
        plt.grid(b=None)
        plt.title("$\sigma_{cen}$", size=size)

        plt.show()


    MW_vlos_1x4_plot()

vproj_maps_1x4_plots()




