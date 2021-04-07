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
import cv2
np.set_printoptions(threshold=np.inf)
sys.path.append('/usr/local/lib/python2.7/site-packages')
def pdist(pt1, pt2):
    x = pt1[0] - pt2[0]
    y = pt1[1] - pt2[1]
    return math.sqrt(math.pow(x, 2) + math.pow(y, 2))

def blob_detection(fname,timestep,initial_timestep):
    # construct the argument parse and parse the arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--image", help = "path to the image")
    args = vars(ap.parse_args())

    # load the image
    #image = cv2.imread(args["image"])
    image = []
    #image.append(cv2.imread('/home/qliu/Dropbox/DataReductionEvaluate-2018-restored/figures/mergeCompress/dpot-1x.png'))
    name_hat = "/media/sf_shared/"
    name_tail1 = "original"
    #name_tail2 = "original"
    name_tail2 = "finer_0"
    print (name_tail2)
    name1 = name_hat+name_tail1+".png"
    name2 = name_hat+name_tail2+".png"
    image.append(cv2.imread("xgc_original.png"))
    image.append(cv2.imread(fname))
    #image.append(cv2.imread("test.png"))
    #image.append(cv2.imread("xgca_256_stats_preserving.png"))
    #image.append(cv2.imread("xgca_8.png"))
    #image.append(cv2.imread("dpot-16.png"))
    #image.append(cv2.imread("dpot-32.png"))
    #image.append(cv2.imread("dpot-64.png"))

    #image.append(cv2.imread("/home/qliu/Downloads/compressed-images/c4-zfp.png"))
    #image.append(cv2.imread("/home/qliu/Downloads/compressed-images/c4-sz.png"))
    #image.append(cv2.imread("/home/qliu/Downloads/compressed-images/c5-zfp.png"))
    #image.append(cv2.imread("/home/qliu/Downloads/compressed-images/c5-sz.png"))
    #image.append(cv2.imread("/home/qliu/Downloads/compressed-images/c6-zfp.png"))
    #image.append(cv2.imread("/home/qliu/Downloads/compressed-images/c6-sz.png"))

    # define the list of boundaries
    boundaries = [
        ([0, 0, 100], [204, 204, 255]), #red 
        ([86, 31, 4], [220, 88, 50]),
        ([25, 146, 190], [62, 174, 250]),
        ([103, 86, 65], [145, 133, 128])
    ]

    (lower, upper) = boundaries[0]
    print (boundaries[0])

    # create NumPy arrays from the boundaries
    lower = np.array(lower, dtype = "uint8")
    upper = np.array(upper, dtype = "uint8")

        # find the colors within the specified boundaries and apply
        # the mask


    # show the images
    #plt.imshow(output)
    #plt.show()

    # Setup SimpleBlobDetector parameters.
    params = cv2.SimpleBlobDetector_Params()

    # Change thresholds
    params.minThreshold = 10;
    params.maxThreshold = 200;
    print (params.thresholdStep)
    # Filter by Area.
    params.filterByArea = 1
    params.minArea = 120

    # Filter by Circularity
    #params.filterByCircularity = True
    #params.minCircularity = 0.1

    # Filter by Convexity
    params.filterByConvexity = 1
    params.minConvexity = 0.3
    #params.maxConvexity = 1


    # Filter by Inertia
    params.filterByInertia = 1
    params.minInertiaRatio = 0.1
    #params.maxInertiaRatio = 1

    # Set up the detector with default parameters.
    detector = cv2.SimpleBlobDetector_create(params)


    mask = []
    output = []
    gray_image = []
    keypoints = []
    abc = cv2.inRange(image[0], lower, upper)
    print(len(image))
    for i in range(len(image)):
        print i
        time1 = time.time()
        mask.append(cv2.inRange(image[i], lower, upper))
        output.append(cv2.bitwise_and(image[i], image[i], mask = mask[i]))
        gray_image.append(cv2.cvtColor(output[i], cv2.COLOR_BGR2GRAY))
        keypoints.append(detector.detect(gray_image[i]))
        time2 = time.time()
        print ('blob detection time', (time2 - time1))
        print ('image %d blob # %d' %(i, len(keypoints[i])))


        total_diameter = 0
        total_blob_area  = 0
        for k in keypoints[i]:
            total_diameter = total_diameter + k.size
            total_blob_area = total_blob_area + 3.14 * math.pow(k.size/2, 2)
        if len(keypoints[i]):
            print ('avg diameter', total_diameter / len(keypoints[i]))
        else:
            print ('ERROR: avg diameter', 0)
        print ('aggregate blob area', total_blob_area)

        if i > 0:
            overlap = 0
            for k in keypoints[i]:
                for p in keypoints[0]:
                    if pdist(k.pt, p.pt) < (k.size + p.size) * 1.0 / 2.0:
                        overlap = overlap + 1.0
                        break
            if len(keypoints[i]):
                print ('overlap ratio', overlap / len(keypoints[i]))
            else:
                print ('ERROR: overlap ratio', 0)
            
        im_with_keypoints = cv2.drawKeypoints(cv2.cvtColor(image[i], cv2.COLOR_BGR2RGB), keypoints[i], np.array([]), (0,0,255), cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
        plt.imshow(im_with_keypoints)
        plt.axis('off')
        if i ==1:
            if timestep!=initial_timestep:
                fp = np.load("/ssd/8192X_PSNR_30/xgc_blob.npz")
                blob_number = fp["blob_number"]
                blob_diameter = fp["blob_diameter"]
                blob_area = fp["blob_area"]
                overlap_ratio = fp["overlap_ratio"]
                blob_number = blob_number.tolist()
                blob_diameter = blob_diameter.tolist()
                blob_area = blob_area.tolist()
                overlap_ratio = overlap_ratio.tolist()
            else:
                blob_number=[]
                blob_diameter=[]
                blob_area=[]
                overlap_ratio=[]
            blob_number.append(len(keypoints[i]))
            blob_diameter.append(total_diameter / len(keypoints[i]))
            blob_area.append(total_blob_area)
            overlap_ratio.append(overlap / len(keypoints[i]))
            np.savez("/ssd/8192X_PSNR_30/xgc_blob.npz", blob_number = blob_number, blob_diameter = blob_diameter, blob_area = blob_area, overlap_ratio = overlap_ratio)
            print "blob_number = ", blob_number
            print "blob_diameter = ", blob_diameter
            print "blob_area = ", blob_area
            print "overlap_ratio = ", overlap_ratio 
        plt.savefig("/ssd/8192X_PSNR_30/xgc_30_120_"+str(timestep)+'_1_blobed.pdf', dpi=600,format='pdf')
        #plt.savefig('/media/sf_shared/blobed/test.pdf', dpi = 600, format='pdf')

        #plt.savefig('xgca_256_after_decimation_blobed.pdf', format='pdf')
        #plt.savefig('xgca_256_stats_preserving_blobed.pdf', format='pdf')
        #plt.show()





    #for kp in keypoints:
    #    print '(x, y)', int(kp.pt[0]), int(kp.pt[1])
    #    print 'diameter', int(kp.size)
    #    print 'strength', kp.response

    #print 'keypoint type', type(keypoints[0])
    # Show keypoints

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

deci_ratio = 8192
initial_timestep = 60
time_interval = 60
reduced_len = 5485
finer_len = 44928785
for i in range(480, 1801, time_interval):
    print "Timestep = ", i+1800
    fname_head ="/ssd/8192X_PSNR_30/xgc_30_120_"+str(i)+"_1_plot"
    fname_results = fname_head + ".npz" 
    fp = np.load(fname_results)
    data = fp['finer']
    data_r = fp['finer_r']
    data_z = fp['finer_z']

    fname_head ="/ssd/8192X_PSNR_30/xgc_30_120_"+str(i)+"_1_psnr"
    fname_results = fname_head + ".npz"
    fp = np.load(fname_results)

    finer = fp['psnr_finer']
    finer_r = fp['psnr_finer_r']
    finer_z = fp['psnr_finer_z']  
    print "len(finer_plot_data)=", len(data), len(data_r), len(data_z)
    print "len(psnr_finer) = ", len(finer)

    filename = "/ssd/reduced_data_xgc.bin"
    f = open(filename, "rb")
    dpot_L1_str=f.read(reduced_len*8)
    #dpot_L1=zfpy._decompress(dpot_L1_compressed, 4, [2496111], tolerance=0.01)
    r_L1_str=f.read(reduced_len*8)
    #r_L1=zfpy._decompress(r_L1_compressed, 4, [2496111], tolerance=0.01)
    z_L1_str=f.read(reduced_len*8)
    #z_L1=zfpy._decompress(z_L1_compressed, 4, [2496111], tolerance=0.01)
    f.close()
    dpot_L1=struct.unpack(str(reduced_len)+'d',dpot_L1_str)
    r_L1=struct.unpack(str(reduced_len)+'d',r_L1_str)
    z_L1=struct.unpack(str(reduced_len)+'d',z_L1_str)

    psnr_start = time.time()
    filename = "/ssd/full_data_xgc.bin"
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
    fname_psnr = "/ssd/8192X_PSNR_30/xgc_30_120_psnr.npz"
    if i == initial_timestep:
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
    fname_head = "/ssd/8192X_PSNR_30/xgc_30_120_macro_"+str(i)
    plot_name = fname_head + ".png"
    plot(data,data_r,data_z,plot_name) 
    blob_detection(plot_name,i,initial_timestep)
    end = time.time()
    print "Analysis time = ", end - start
    print "*********************************************************"

#fp = np.load("xgc/xgc_60_120_75_1.npz")
#data = fp['finer']
#r = fp['finer_r']
#z = fp['finer_z']    
#plot(data,r,z,"12345.png") 

