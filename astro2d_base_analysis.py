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
import cv2
from scipy.spatial import Delaunay
from skimage.measure import compare_ssim
np.set_printoptions(threshold=np.inf)

#image.append(cv2.imread("test.png"))
#image_original = cv2.imread("astro2d/test2/original.png")
#image_reduced = cv2.imread("astro2d/test2/reduced.png")
#image_reduced = cv2.imread("xgc/original.png")
def calculate_dice(original_imag_name, reduced_imag_name):
    image_original = cv2.imread(original_imag_name,cv2.IMREAD_GRAYSCALE)
    image_reduced = cv2.imread(reduced_imag_name,cv2.IMREAD_GRAYSCALE)

    #image_original = cv2.cvtColor(image_ori, cv2.COLOR_BGR2GRAY)
    #image_reduced = cv2.cvtColor(image_red, cv2.COLOR_BGR2GRAY)
    cnt=0
    for i in range(np.shape(image_original)[0]):
        for j in range(np.shape(image_original)[1]):
            if image_original[i][j] == image_reduced[i][j]:
                cnt+=1
    return 2*cnt/(np.shape(image_original)[0]*np.shape(image_original)[1] + np.shape(image_reduced)[0]*np.shape(image_reduced)[1])

    #image_original = image_original/255
    #image_reduced = image_reduced/255
    #union = image_original * image_reduced
    #dice = 2*np.sum(union)/(np.sum(image_original)+np.sum(image_reduced))
    #return dice



def SSIM(original_imag_name, reduced_imag_name):
    image_original = cv2.imread(original_imag_name)
    image_reduced = cv2.imread(reduced_imag_name)
    
    #Convert the images to grayscale
    gray_original = cv2.cvtColor(image_original, cv2.COLOR_BGR2GRAY)
    gray_reduced = cv2.cvtColor(image_reduced, cv2.COLOR_BGR2GRAY)

    (score, diff) = compare_ssim(gray_original, gray_reduced, full=True)
    diff = (diff * 255).astype("uint8")

    print("SSIM: {}".format(score))
    return score

def psnr_c(original_data, base, leveldata_len, deci_ratio, level_id):
    for i in range(len(leveldata_len)-level_id,len(leveldata_len)):
        leveldata=np.zeros(leveldata_len[i])
        for j in range(leveldata_len[i]):
            index1=j//deci_ratio
            index2=j%deci_ratio
            if index1!= len(base)-1:
                if index2 != 0:
                    leveldata[j]=(base[index1]+base[index1+1])*index2/deci_ratio
                else:
                    leveldata[j]=base[index1]
            else:
                if index2 != 0:
                    leveldata[j]=(base[index1]*2)*index2/deci_ratio
                else:
                    leveldata[j]=base[index1]
        base=leveldata
    if len(base) !=len(original_data):
        print "len(leveldata) !=len(original_data)"

    MSE = 0.0
    for i in range(len(original_data)):
        #print i, original_data[i]-base[i]
        MSE=(original_data[i]-base[i])*(original_data[i]-base[i])+MSE
    MSE=MSE/len(original_data)
    if MSE < 1e-6:
        MSE=0.0
    if MSE ==0.0:
        print "Couldn't get PSNR of two identical data."
        return 0
    else:
        psnr=10*math.log(np.max(original_data)**2/MSE,10)
    #print "psnr=",psnr
    return psnr

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
    print len(r), len(z), len(conn), len(data)
    plt.tricontourf(r, z, conn, data,cmap=plt.cm.jet, levels=np.linspace(np.min(data),np.max(data),num=25));
    #plt.colorbar();
    plt.xticks([])
    plt.yticks([])
    for key, spine in ax.spines.items():
        # 'left', 'right', 'bottom', 'top'
        if key == 'right' or key == 'top' or key == 'left' or key == 'bottom':
            spine.set_visible(False)
    plt.savefig(filename, format='png')

deci_ratio = 8192
time_interval = 60
reduced_len = 5787
finer_len = 47403736
fp = np.load("/ssd/astro2d_30_120_base_plot.npz")
data = fp['finer']
data_r = fp['finer_r']
data_z = fp['finer_z']

fp = np.load("/ssd/astro2d_30_120_base_psnr.npz")
finer = fp['psnr_finer']
finer_r = fp['psnr_finer_r']
finer_z = fp['psnr_finer_z']
print "len(finer_plot_data)=", len(data), len(data_r), len(data_z)
print "len(psnr_finer) = ", len(finer)

filename = "/ssd/reduced_data_astro2d.bin"
f = open(filename, "rb")
dpot_L1_str=f.read(reduced_len*8)
r_L1_str=f.read(reduced_len*8)
z_L1_str=f.read(reduced_len*8)
f.close()
dpot_L1=struct.unpack(str(reduced_len)+'d',dpot_L1_str)
r_L1=struct.unpack(str(reduced_len)+'d',r_L1_str)
z_L1=struct.unpack(str(reduced_len)+'d',z_L1_str)

filename = "/ssd/full_data_astro2d.bin"
f = open(filename, "rb")
dpot_str = f.read(finer_len*8)
r_str = f.read(finer_len*8)
z_str = f.read(finer_len*8)
f.close()
number_of_original_elements = str(finer_len)
dpot=struct.unpack(number_of_original_elements+'d',dpot_str)
r=struct.unpack(number_of_original_elements+'d',r_str)
z=struct.unpack(number_of_original_elements+'d',z_str)
data_len=[len(dpot_L1),len(dpot)]
print data_len
if len(finer)!=len(dpot):
    print "This timestep doesn't read delta\n"
    if len(finer) != len(dpot_L1):
        print "ERROR: length of finer data for psnr is wrong\n"
        print "len(finer)=%d len(dpot)=%d\n"%(len(finer), len(dpot))
if len(finer) == len(dpot):
    print "No augment and PSNR"
    psnr_finer=psnr_c(dpot, finer, data_len, deci_ratio, 0)
else:
    print "Augment and PSNR"
    psnr_finer = psnr_c(dpot,dpot_L1,data_len,deci_ratio, 1)
print "finer PSNR=",psnr_finer

plot_name = "/hdd/astro2d_30_120_macro_base.png"
plot(data,data_r,data_z,plot_name) 
score = SSIM("astro2d_original.png", plot_name)
dice_co = calculate_dice("astro2d_original.png", plot_name)
print "SSIM score = ",score
print "Dice similarity coefficient",dice_co
