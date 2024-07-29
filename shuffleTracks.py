import pandas as pd
import numpy as np
import sys,os
import cv2 
from tqdm import tqdm 

from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.collections import LineCollection

from processTracks import get_stats,get_meta,get_tracks,get_distance

if __name__ == "__main__":
    #csv_file = sys.argv[1]
    file_dir = sys.argv[1]
    pid = sys.argv[2]
    d_dist = []
    a_dist = []
    v_dist = []
    count = 0
    for f in sorted(os.listdir(file_dir)):
        if pid not in f:
            continue
        count += 1
        if count > 10:
            pass
            #break
        csv_file = '/'.join([file_dir,f])
        track_array,track_polar,(n_frames,n_fish,fishIDs) = get_tracks(csv_file)

        first_frame = np.argmax(~np.isnan(track_array[:,0]))
        clean_array = track_array[first_frame:]
        clean_polar = track_polar[first_frame:]

        shifted_array = np.empty_like(clean_array)
        shifted_polar = np.empty_like(clean_array)
       
        iterations = 30
        dists = np.empty(iterations)
        angs = np.empty(iterations)
        vels = np.empty(iterations)

        for i in tqdm(range(iterations)):
            for f in range(4):
                shift = np.random.randint(3600,len(shifted_array) - 3600)
                shifted_array[:,f] = np.roll(clean_array[:,f],shift)
                shifted_polar[:,f] = np.roll(clean_polar[:,f],shift)
                shifted_stats,_ = get_stats(shifted_array,shifted_polar)
                dists[i] = shifted_stats['dist_mean']
                angs[i] = shifted_stats['angleC_mean']
                vels[i] = shifted_stats['velC_mean']
        real_stats,_ = get_stats(clean_array,clean_polar)
        real_dist = real_stats['dist_mean']
        real_ang = real_stats['angleC_mean']
        real_vel = real_stats['velC_mean']

        d_dist.append(np.mean(dists) - real_dist)
        a_dist.append(real_ang - np.mean(angs))
        v_dist.append(real_vel - np.mean(vels))
        if False: ## Just always 0th percentile
            percentile = sum(dists < real_dist) / iterations
            percentileA = sum(angs < real_ang) / iterations
            percentileV = sum(vels < real_vel) / iterations
            print(percentile,percentileA,percentileV)
        if False:
            fig,axes = plt.subplots(1,3)
            axes[0].hist(dists)
            axes[0].axvline(real_dist)

            axes[1].hist(angs)
            axes[1].axvline(real_ang)

            axes[2].hist(vels)
            axes[2].axvline(real_vel)

            plt.show()
        #import pdb;pdb.set_trace()
fig,axes = plt.subplots(3,1)
axes[0].plot(d_dist)
axes[1].plot(a_dist)
axes[2].plot(v_dist)
plt.show()
print('Done!')
