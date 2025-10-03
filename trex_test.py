import numpy as np

#Add sources to path
import sys
#sys.path.insert(1, '../src/')

import trex_python.networkAlgo as na          #Graph management
import trex_python.regularization_gmm as rg   #RegGMM algorithm
import matplotlib.pyplot as plt
import trex_python.utility as utility         #For plots and internal functions

from scipy.io import FortranFile

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
import pickle
from mpl_toolkits import mplot3d
from ipywidgets import interact
import sys



unit_kpc = 2.2755E27 / 3.08567758128E21
unit_l = 0.227550518265197E+28
unit_d = 0.264088847033921E-29
unit_t = 0.455724529603677E+18
unit_p = unit_d * (unit_l / unit_t) ** 2

# Plots parameters
axSize = [0.1,0.1,0.8,0.8]
#cbarPos = [0.9, 0.1, 0.025, 0.8]
figSize = (6, 5)

#test=np.zeros(10)
#print(test)
#np.save('test.npy',test)

#data = np.load('test.npy')
#print(data)

#sys.exit()

#Loading data
#data = np.load('./trex/toy_dataset.npy')

#print("len data",len(data))

#sys.exit()

def load_dm():
    file="virgo_xyz_dm_low_res.dat"

    d = FortranFile(file, 'r')

    print('fichier dm ouvert')

    ncell = d.read_ints()
    print("n dm cell", ncell)

    x = d.read_reals()
    y = d.read_reals()
    z = d.read_reals()
    vx = d.read_reals()
    vy = d.read_reals()
    vz = d.read_reals()
    m = d.read_reals()

    step = 1
    n = int(len(x) / step)

    xplot = np.array([x[step * i] for i in range(n)])
    yplot = np.array([y[step * i] for i in range(n)])
    zplot = np.array([z[step * i] for i in range(n)])

    data=np.array([xplot,yplot,zplot])
    #data=np.array([xplot,yplot])

    print("data len",len(data.T))

def load_gal_list():
    gal=np.loadtxt("/data/cluster/byopic/SIMS/VirgoClone/HighRes/list_gal_251.dat_js_nocontam")
    #print("gal pos",gal[0,3:6])

    cond = np.logical_and(gal[:, 3] > 0.46, np.logical_and(gal[:, 3] < 0.51, np.logical_and(gal[:, 4] > 0.48,np.logical_and(gal[:, 4] < 0.53, np.logical_and(gal[:, 5] > 0.47, gal[:, 5] < 0.52)))))
    x_gal = gal[:,3][cond]
    #print("xgal",x_gal)
    #sys.exit()
    y_gal = gal[:,4][cond]
    z_gal = gal[:,5][cond]

    x_gal = (x_gal - 0.5) * (unit_l / 3.08567758128E21)
    y_gal = (y_gal - 0.5) * (unit_l / 3.08567758128E21)
    z_gal = (z_gal - 0.5) * (unit_l / 3.08567758128E21)

    data=np.array([x_gal,y_gal,z_gal])
    #data=data.T

    return data


    #sys.exit()

data=load_gal_list()

#fig = plt.figure()
#ax = plt.axes(projection='3d')
#cond=np.logical_and(xplot>xl,np.logical_and(xplot<xu,np.logical_and(yplot>yl,np.logical_and(yplot<yu,np.logical_and(zplot>zl,zplot<zu)))))
#ax.scatter3D(xplot, yplot, zplot, s=2, c='black')
#ax.set_xlim(0.4711, 0.4981)
#ax.set_ylim(0.49459, 0.52159)
#ax.set_zlim(0.4833, 0.5103)

#plt.legend()
#plt.show()

#sys.exit()


#print('test')

#Plot
#fig, ax = plt.subplots(figsize=figSize)
#ax.scatter(data[0], data[1], s=5, color='k')
#plt.show()

#sys.exit()

data=data.T

# Set parameters of the T-ReX algorithm
param = rg.param_()
A0=0.05
#A0=7
param.set(data=data, topology=na.buildMST_fast, verbose=0, covariance_type='spherical',
          covariance_update='fixed', lam=5, lam_sig=5, A0=A0, lam_pi=1, maxIter=300, eps=1e-5, denoisingParameter=0)

# Set parameters for the background
param.background_noise = False
if param.background_noise:
    param.alpha = 0.25  # Initilisation of the noise level (true is 0.25)
else:
    param.alpha = 0
param.domain_area = np.prod(np.max(data, axis=0) - np.min(data, axis=0))  # Estimate of the volume containing the data

#param.sig = 0.01  #Initialisation of the variance of Gaussian components

X = data      # Datapoints used for the computation

# Initialise the graph with 100 random datapoints
np.random.seed(10)
M = X[np.random.choice(np.arange(len(X)), 100, replace=False)]

#Get the regularised graph
regGraph, _, _, _ = rg.regGMM(X, M, param, computeObj=0, display_progress_bar=1)
opt_edges = regGraph.edges
F = regGraph.F

#Compute XYZ positions of edges for plots
Xe_fil, Ye_fil, Ze_fil = na.compute_XYZ_of_edges(opt_edges, F)


def plot_reg_graph():
    #Plot the regularised graph
    fig, ax = utility.getFigure(figSize, axSize)
    plt.scatter(data.T[0], data.T[1], color = 'k', s=3, alpha=0.5)

    plt.plot(Xe_fil, Ye_fil, lw=3, c='crimson', alpha=0.8, zorder=1)
    plt.title('Regularized MST ($\lambda_\mu$={}, $\sigma_0^2$={:.2f}, $A_0$={})' .format(param.lam, param.sig, A0))

    #Plot covariance ellipses
    for j in range(len(F)):
        if isinstance(regGraph.sig, np.ndarray):
            if param.covariance_update == 'adaptative':
                covMat = np.eye(2)*regGraph.sig[j]
            else:
                covMat = np.eye(2)*regGraph.sig[j]
        else:
            covMat = np.eye(2)*regGraph.sig
        ell = utility.get_cov_ellipse(covMat, F[j][0:2], 0.395, alpha=1, zorder=-1)   #0.395 to see 1-sigma circles
        ell.set_facecolor('silver')
        ax.add_artist(ell)
    plt.show()
    sys.exit()

def plot_reg_graph_3d():

    #print("Xe_fil",Xe_fil)
    #sys.exit()
    for i in range(len(Xe_fil)) :
        if (Xe_fil[i]==None) :
            Xe_fil[i]=np.nan

    for i in range(len(Ye_fil)) :
        if (Ye_fil[i]==None) :
            Ye_fil[i]=np.nan

    for i in range(len(Ze_fil)) :
        if (Ze_fil[i]==None) :
            Ze_fil[i]=np.nan

    #step = 20
    #n = int(len(x) / step)

    #xplot = np.array([x[step * i] for i in range(n)])
    #yplot = np.array([y[step * i] for i in range(n)])
    #zplot = np.array([z[step * i] for i in range(n)])

    Virgo = np.array([0.48461068, 0.50809848, 0.49687076])
    MW = np.array([0.5, 0.5, 0.5])
    fil = np.array([0.497, 0.4984, 0.4996])
    filx = np.array([0.4839, 0.5087, 0.51])
    fily = np.array([0.476, 0.497, 0.4953])
    opposfil = 2 * Virgo - fil
    opposMW = 2 * Virgo - MW
    opposfilx = 2 * Virgo - filx
    opposfily = 2 * Virgo - fily
    mainfil = np.array([opposfil, fil])
    vmw = np.array([opposMW, MW])
    vfilx = np.array([opposfilx, filx])
    vfily = np.array([opposfily, fily])

    mainfil = (mainfil - 0.5) * (unit_l / 3.08567758128E21)
    vfilx = (vfilx - 0.5) * (unit_l / 3.08567758128E21)
    vfily = (vfily - 0.5) * (unit_l / 3.08567758128E21)
    vmw = (vmw - 0.5) * (unit_l / 3.08567758128E21)

    fig = plt.figure()
    ax = plt.axes(projection='3d')#,computed_zorder=False)
    # cond=np.logical_and(xplot>xl,np.logical_and(xplot<xu,np.logical_and(yplot>yl,np.logical_and(yplot<yu,np.logical_and(zplot>zl,zplot<zu)))))
    ax.scatter3D(data.T[0], data.T[1], data.T[2], s=1, c='black',alpha=0.5,zorder=1)
    #ax.scatter3D(xplot, yplot, zplot, s=1, c='black',alpha=0.15,zorder=1)
    ax.plot(Xe_fil,Ye_fil,Ze_fil,lw=5, c='crimson', alpha=1,zorder=2)
    #ax.plot(mainfil[:, 0], mainfil[:, 1], mainfil[:, 2], color='royalblue', lw=8, label="Fil")
    #ax.plot(vfilx[:, 0], vfilx[:, 1], vfilx[:, 2], color='red', lw=8, label="Filx")
    #ax.plot(vfily[:, 0], vfily[:, 1], vfily[:, 2], color='violet', lw=8, label="Fily")
    #ax.plot(vmw[:, 0], vmw[:, 1], vmw[:, 2], color='green', lw=8, label="Cen")
    plt.title('Regularized MST ($\lambda_\mu$={}, $\sigma_0^2$={:.2f}, $A_0$={})'.format(param.lam, param.sig, A0))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")


    #for j in range(len(F)):
        #if isinstance(regGraph.sig, np.ndarray):
        #    if param.covariance_update == 'adaptative':
        #        covMat = np.eye(2)*regGraph.sig[j]
        #    else:
        #        covMat = np.eye(2)*regGraph.sig[j]
        #else:
        #    covMat = np.eye(2)*regGraph.sig
        #ell = utility.get_cov_ellipse(covMat, F[j][0:2], 0.395, alpha=1, zorder=-1)   #0.395 to see 1-sigma circles
        #ell.set_facecolor('silver')
        #ax.add_artist(ell)

    plt.show()
    sys.exit()

plot_reg_graph_3d()
#plot_reg_graph()



filaments, dataSeg, affil_clusters, res = rg.filamentFinding(data, regGraph, param)

 #Plots
d = na.getDegree(np.array(opt_edges), len(F))  #Function take numpy array as input, not list

#Plot each branch individually
fig, ax = utility.getFigure(figSize, axSize)
ax.scatter(data.T[0], data.T[1], s=1, color='k')
length_filaments = []
for f in filaments:
    if f.type==1:
        clr = 'red'
    elif f.type==2:
        clr = 'green'
    else:
        clr = 'blue'
    length_filaments.append(f.length)
    ax.plot(F[f.nodes].T[0], F[f.nodes].T[1], '-o', markersize=2) #, color=clr)

unit="kpc"

ax.set_title('Filaments from regularized graph')
ax.set_xlabel('x ' + unit); ax.set_ylabel('y ' + unit)


#Plot the segregation of input points according detected filaments
fig, ax = utility.getFigure(figSize, axSize)
plt.scatter(F.T[0][(d>=3)], F.T[1][(d >= 3)], color='green', alpha=1, s=100,zorder=3)
plt.scatter(F.T[0][(d==1)], F.T[1][(d==1)], color='red', alpha=1, s=100,zorder=4)
plt.plot(Xe_fil, Ye_fil, lw = 1.5, c = 'black', alpha = 0.8, zorder=2)
plt.xlabel("x " + unit)
plt.ylabel("y " + unit)
plt.legend(('Bifurcation nodes', 'Extremity nodes', 'Edges', 'Junction nodes', 'Used data points'))

vals = np.linspace(0,1,256)
np.random.shuffle(vals)
random_cmap = plt.cm.colors.ListedColormap(plt.cm.jet(vals))

ax.scatter(data.T[0], data.T[1], c='gray', s=2, alpha=0.5)
ax.scatter(data[dataSeg>-1].T[0], data[dataSeg>-1].T[1], c=dataSeg[dataSeg>-1], s=1, alpha=0.8, cmap=random_cmap, zorder=1)
ax.scatter(data[dataSeg==-1].T[0], data[dataSeg==-1].T[1], c='grey', s=1, alpha=0.8, cmap=random_cmap)

plt.show()