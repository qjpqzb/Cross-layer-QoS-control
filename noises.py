# -*- coding: utf-8 -*-
from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.tri import triangulation
import numpy as np
import math
import time
import struct
import sys
import os
import argparse
import datetime
from pytz import timezone
np.set_printoptions(threshold=np.inf)

def start_analysis(read_size, interval, error_target_index, docker_path_id):
    docker_cgroup_path = "/sys/fs/cgroup/blkio/docker/" + docker_path_id + "/blkio.weight"
    if sys.argv[2] == "on":
        os.system("echo %d > %s"%(weight, docker_cgroup_path))
    else:
        os.system("echo 100 > %s"%(docker_cgroup_path))
    os.system("cat %s"%(docker_cgroup_path))
    filename_d = "/hdd/noise1.bin"
    x = []
    y = []
    io_time = []
    bandwidth = []
    for i in range(37):
        start1 = time.time()
        print "%s s\n"%(i*interval)
        x.append(i*interval)
        y.append(0)
        x.append(i*interval)
        start = time.time()
        f = open(filename_d, "rb")
        delta_str = f.read(read_size*1024*1024)
        f.close()
        end = time.time()
        io_time.append(end - start)
        x.append(i*interval+end-start)
        x.append(i*interval+end-start)
        #print "Read time = %f s\n"%(end - start)
        bw = read_size/(end-start)
        bandwidth.append(bw)
        print "Percievd bandwidth = %f MB/s"%(bw)
        y.append(bw)
        y.append(bw)
        y.append(0)
        end1 = time.time()
        print "Analysis time = ", end1- start1
        if end1- start1 > interval:
            print "Analysis time is larger than interval!\n"
        time.sleep(interval - (end1 - start1))

    print "x2 = ",x
    print "app2 = ",y
    print "x2_io_time = ", io_time
    print "x2_bandwidth = ", bandwidth

def main():
    tag = 1
    read_size = 1024
    interval = 100
    docker_path_id = "f6af9c2810d799932818a10153fa8a1f7f3e7955f9d1a7c26af45844113cea05"
    while (tag):
        est = timezone('EST')
        now_time = datetime.datetime.now(est)
        hour = now_time.hour
        minutes = now_time.minute
        seconds = now_time.second
        #print "hour = %d, minutes = %d, seconds = %d"%(hour, minutes, seconds)
        #print "input:",sys.argv[1], sys.argv[2], sys.argv[3]
        if sys.argv[3] == 'now':
            start_analysis(read_size, interval,int(sys.argv[1]), docker_path_id)
            tag = 0
            break
        if hour == int(sys.argv[3]) and minutes == int(sys.argv[4]) and seconds == int(sys.argv[5]):
            start_analysis(read_size, interval,int(sys.argv[1]), docker_path_id)
            tag = 0
            break

if __name__ == "__main__":
    main()


