# -*- coding: utf-8 -*-
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
np.set_printoptions(threshold=np.inf)


def calc_area(p1, p2, p3):
    (x1, y1), (x2, y2), (x3, y3) = p1,p2,p3
    #return 0.5 * abs(x2 * y3 + x1 * y2 + x3 * y1 - x3 * y2 - x2 * y1 - x1 * y3)
    area = abs((x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))*0.5)
    return area

def high_pressure_area(dpot,R,Z,thre):
    points = np.transpose(np.array([Z,R]))
    Delaunay_o = Delaunay(points)
    conn=Delaunay_o.simplices
    area=0.0
    total_force=0.0
    #cont=0
    #conn_p=[]
    #rz_p=[]
    #cont_l=0
    #cont_s=0
    for i in conn:
        index1=i[0]
        index2=i[1]
        index3=i[2]
        #if (dpot[index1]>thre and dpot[index2]>thre) or (dpot[index1]>thre and dpot[index3]>thre) or (dpot[index2]>thre and dpot[index3]>thre):
        if dpot[index1]>thre or dpot[index2]>thre or dpot[index3]>thre:
        #if dpot[index1]>thre and 
        #if (dpot[index1] + dpot[index2] + dpot[index3])/3 >135 :
            #if abs(dpot[index1]-thre) >33 or abs(dpot[index2]-thre) >33 or abs(dpot[index3]-thre) >33:
                #print (dpot[index1]-thre,dpot[index2]-thre,dpot[index3]-thre)
            each_area = calc_area((R[index1],Z[index1]),(R[index2],Z[index2]),(R[index3],Z[index3]))
            area += each_area
            total_force += (dpot[index1] + dpot[index2] + dpot[index3])*each_area/3.0
            #cont=cont+1
    return area,total_force

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

thre = 150
deci_ratio = 8192
reduced_len = 3756
finer_len = 30764603
initial_step = 2460 - 1800
tt_a = [1980, 2100, 2160, 2340, 2400, 2520, 2700, 2880, 3000, 3060, 3240, 3300, 3420, 3600]
tt_a = [2460, 2640, 2760, 2940]
tt_b = [i-1800 for i in tt_a]
#tt_b = [180]
for i in tt_b:
#for i in range(60, 1860, 60):
    print "timestep = ",1800+i
    fname_head ="/ssd/aug_PSNR/cfd_30_120_"+str(i)+"_1_plot"
    fname_results = fname_head + ".npz"
    fp = np.load(fname_results)
    data = fp['finer']
    data_r = fp['finer_r']
    data_z = fp['finer_z']

    fname_head ="/ssd/aug_PSNR/cfd_30_120_"+str(i)+"_1_psnr"
    fname_results = fname_head + ".npz"
    fp = np.load(fname_results)    
    finer = fp['psnr_finer']
    finer_r = fp['psnr_finer_r']
    finer_z = fp['psnr_finer_z']
    print "len(finer_plot_data)=", len(data), len(data_r), len(data_z)
    print "len(psnr_finer) = ", len(finer)

    filename = "/ssd/reduced_data_cfd.bin"
    f = open(filename, "rb")
    dpot_L1_str=f.read(reduced_len*8)
    r_L1_str=f.read(reduced_len*8)
    z_L1_str=f.read(reduced_len*8)
    f.close()
    dpot_L1=struct.unpack(str(reduced_len)+'d',dpot_L1_str)
    r_L1=struct.unpack(str(reduced_len)+'d',r_L1_str)
    z_L1=struct.unpack(str(reduced_len)+'d',z_L1_str)

    psnr_start = time.time()
    filename = "/ssd/full_data_cfd.bin"
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
        psnr_finer=psnr_c(dpot, finer, data_len, deci_ratio, 0)
    else:
        psnr_finer = psnr_c(dpot,dpot_L1,data_len,deci_ratio, 1)
    print "finer PSNR=",psnr_finer
    fname_head = "/ssd/aug_PSNR/cfd_30_120_" + str(deci_ratio)+"X"
    fname_psnr = fname_head + "_psnr.npz"
    if i == initial_step:
        np.savez(fname_psnr, psnr = [psnr_finer])
        print "Total PSNR =", [psnr_finer]
    else:
        fpp = np.load(fname_psnr)
        s_psnr = fpp["psnr"]
        s_psnr = s_psnr.tolist()
        s_psnr.append(psnr_finer)
        np.savez(fname_psnr, psnr = s_psnr)
        print "Total PSNR =", s_psnr
    psnr_end = time.time()
    print "Time for calculate PSNR = ",psnr_end-psnr_start

    start = time.time()
    area, force = high_pressure_area(data, data_r, data_z, thre)
    end = time.time()
    print "Analysis time = ", end - start
    macro_name = fname_head + "_macro.npz"
    if i!= initial_step:
        fpf = np.load(macro_name)
        total_area=fpf['total_area']
        total_force=fpf['total_force']
        total_area = total_area.tolist()
        total_force = total_force.tolist()
    else:
        total_area=[]
        total_force=[]
    total_area.append(area)
    total_force.append(force)
    np.savez(macro_name, total_area = total_area, total_force = total_force)
    print "total_area = ", total_area
    print "total_force = ", total_force
    print "*****************************************"
    



