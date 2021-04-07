from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.tri import triangulation
from scipy.spatial import Delaunay
from scipy import stats
import numpy as np
import math
import time
import struct
import sys
import os
import subprocess
from sklearn.cluster import KMeans
#from scipy.fftpack import fft
import scipy.fftpack as fft
import Queue
from threading import Thread
import argparse
from scipy.spatial import Delaunay
np.set_printoptions(threshold=np.inf)

def plot(data,r,z,filename):
    points = np.transpose(np.array([z,r]))
    Delaunay_t = Delaunay(points)
    conn=Delaunay_t.simplices
    fig,ax=plt.subplots(figsize=(8,8))
    plt.rc('xtick', labelsize=26)          # fontsize of the tick labels
    plt.rc('ytick', labelsize=26)

    axis_font = {'fontname':'Arial', 'size':'38'}

    #plt.xlabel('R', **axis_font)
    #plt.ylabel('Z',**axis_font )
    plt.tricontourf(r, z, conn, data,cmap=plt.cm.jet, levels=np.linspace(np.min(data),np.max(data),num=25));
    #plt.colorbar();
    plt.xticks([])
    plt.yticks([])
    for key, spine in ax.spines.items():
        # 'left', 'right', 'bottom', 'top'
        if key == 'right' or key == 'top' or key == 'left' or key == 'bottom':
            spine.set_visible(False)
    plt.savefig(filename, format='png')

full_len = 44928785
filename = "/ssd/full_data_xgc.bin"
f = open(filename, "rb")
dpot_L1_str=f.read(full_len*8)
r_L1_str=f.read(full_len*8)
z_L1_str=f.read(full_len*8)
f.close()
end_ssd = time.time()
dpot_L1=struct.unpack(str(full_len)+'d',dpot_L1_str)
r_L1=struct.unpack(str(full_len)+'d',r_L1_str)
z_L1=struct.unpack(str(full_len)+'d',z_L1_str)
start = time.time()
plot(dpot_L1,r_L1,z_L1,"xgc_original.png")
end = time.time()
print "Plot time = ", end -start

