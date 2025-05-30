import pandas as pd
import numpy as np
import sys,os
import cv2 
from tqdm import tqdm 

from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.collections import LineCollection

from processTracks import get_stats,get_meta,get_tracks,get_distance
from joblib import Parallel, delayed

iterations = 100
current_pid = 0
d_dist = []
all_pis = []

def shuffle_stats(clean_array,clean_polar):
    shifted_array = np.full_like(clean_array,np.nan)
    shifted_polar = np.full_like(clean_polar,np.nan)
    for f in range(4):
        shift = np.random.randint(3600,len(shifted_array) - 3600)
        shifted_array[:,f] = np.roll(clean_array[:,f],shift,axis=0)
        shifted_polar[:,f] = np.roll(clean_polar[:,f],shift,axis=0)
    shifted_stats,_ = get_stats(shifted_array,shifted_polar)
    return shifted_stats

if __name__ == "__main__":
    #csv_file = sys.argv[1]
    file_dir = sys.argv[1]
    all_ds = np.full([14,100],np.nan)
    all_as = np.full([14,100],np.nan)

    all_vs = np.full([14,100],np.nan)
    all_rvs = np.full([14,100],np.nan)

    all_ps = np.full([14,100],np.nan)
    all_rps = np.full([14,100],np.nan)
    pi_count = -1
    for f in sorted(os.listdir(file_dir)):
        pid = f.split('_')[1]
        if 'pi' not in pid:
            continue
        if pid != current_pid:
            print('starting pid:',pid)
            if pid != 0:
                all_pis.append(current_pid)
            pi_count += 1
            count = 0
            #d_dist = []
            #a_dist = []
            #v_dist = []
            current_pid = pid
        count += 1

        csv_file = '/'.join([file_dir,f])
        track_array,track_polar,(n_frames,n_fish,fishIDs) = get_tracks(csv_file)
        first_frame = np.argmax(~np.isnan(track_array[:,0]))
        clean_array = track_array[first_frame:]
        clean_polar = track_polar[first_frame:]

        shifted_array = np.empty_like(clean_array)
        shifted_polar = np.empty_like(clean_array)
       
        dists = np.empty(iterations)
        angs = np.empty(iterations)
        vels = np.empty(iterations)
        pDists = np.empty(iterations)

        raw_vels = np.empty(iterations)
        raw_pDists = np.empty(iterations)

        if True:
            all_stats = Parallel(n_jobs=10)(delayed(shuffle_stats)(clean_array,clean_polar) for i in tqdm(range(iterations)))
            for i in range(iterations):
                dists[i] = all_stats[i]['dist_mean']
                angs[i] = all_stats[i]['angleC_mean']
                vels[i] = all_stats[i]['velC_mean']
                raw_vels[i] = all_stats[i]['vel_mean']
                raw_pDists[i] = all_stats[i]['pDist_mean']
                pDists[i] = all_stats[i]['pDistC_mean']
        else:
            for i in tqdm(range(iterations)):
                for f in range(4):
                    shift = np.random.randint(3600,len(shifted_array) - 3600)
                    shifted_array[:,f] = np.roll(clean_array[:,f],shift)
                    shifted_polar[:,f] = np.roll(clean_polar[:,f],shift)
                shifted_stats,_ = get_stats(shifted_array,shifted_polar)
                dists[i] = shifted_stats['dist_mean']
                angs[i] = shifted_stats['angleC_mean']
                vels[i] = shifted_stats['velC_mean']
                raw_vels[i] = shifted_stats['vel_mean']

                pDists[i] = shifted_stats[i]['pDistC_mean']
                raw_pDists[i] = shifted_stats[i]['pDist_mean']
        real_stats,_ = get_stats(clean_array,clean_polar)
        real_dist = real_stats['dist_mean']
        real_ang = real_stats['angleC_mean']

        real_vel = real_stats['velC_mean']
        real_rvel = real_stats['vel_mean']

        real_pDist = real_stats['pDistC_mean']
        real_rpDist = real_stats['pDist_mean']

        #d_dist.append(np.mean(dists) - real_dist)
        #a_dist.append(real_ang - np.mean(angs))
        #v_dist.append(real_vel - np.mean(vels))
        all_ds[pi_count,count] = np.mean(dists) - real_dist
        all_as[pi_count,count] = real_ang - np.mean(angs)
        all_vs[pi_count,count] = real_vel - np.mean(vels)
        all_rvs[pi_count,count] = real_rvel - np.mean(raw_vels)

        all_ps[pi_count,count] = real_pDist - np.mean(pDists)
        all_rps[pi_count,count] = real_rpDist - np.mean(raw_pDists)



        if False: ## Just always 0th percentile
            percentile = sum(dists < real_dist) / iterations
            percentileA = sum(angs < real_ang) / iterations
            percentileV = sum(vels < real_vel) / iterations
            percentileP = sum(pDists < real_pDists) / iterations
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
    import pdb;pdb.set_trace()
    fig,axes = plt.subplots(2,3)
    for i in range(len(all_pis)):
        d_dist,a_dist,v_dist,p_dist = all_ds[i],all_as[i],all_vs[i],all_ps[i]
        rp_dist = all_rps[i]
        rv_dist = all_rvs[i]

        axes[0,0].plot(d_dist,alpha=0.1,color='black')
        axes[0,1].plot(rv_dist,alpha=0.1,color='black')
        axes[0,2].plot(v_dist,alpha=0.1,color='black')

        axes[1,0].plot(rp_dist,alpha=0.1,color='black')
        axes[1,1].plot(p_dist,alpha=0.1,color='black')
        axes[1,2].plot(a_dist,alpha=0.1,color='black')

    mean_ds = np.nanmean(all_ds,axis=0)
    std_ds = np.nanstd(all_ds,axis=0)

    mean_as = np.nanmean(all_as,axis=0)
    std_as = np.nanstd(all_as,axis=0)

    mean_vs = np.nanmean(all_vs,axis=0)
    std_vs = np.nanstd(all_vs,axis=0)

    mean_rvs = np.nanmean(all_rvs,axis=0)
    std_rvs = np.nanstd(all_rvs,axis=0)

    mean_ps = np.nanmean(all_ps,axis=0)
    std_ps = np.nanstd(all_ps,axis=0)

    mean_rps = np.nanmean(all_rps,axis=0)
    std_rps = np.nanstd(all_rps,axis=0)

    xs = np.arange(100)
    axes[0,0].set_xlim([0,54])
    axes[0,0].plot(mean_ds,color='black')
    axes[0,0].fill_between(xs,mean_ds - std_ds,mean_ds + std_ds,alpha=0.3,color='gray')
    axes[0,0].axhline(0,linestyle=':',color='black')
    axes[0,0].set_ylabel('Distance:Coh (pixels)')

    axes[0,1].set_xlim([0,54])
    axes[0,1].plot(mean_rvs,color='black')
    axes[0,1].fill_between(xs,mean_rvs - std_rvs,mean_rvs + std_rvs,alpha=0.3,color='gray')
    axes[0,1].axhline(0,linestyle=':',color='black')
    axes[0,1].set_ylabel('Distance:vel (pixels)')

    axes[0,2].set_xlim([0,54])
    axes[0,2].plot(mean_vs,color='black')
    axes[0,2].fill_between(xs,mean_vs - std_vs,mean_vs + std_vs,alpha=0.3,color='gray')
    axes[0,2].axhline(0,linestyle=':',color='black')
    axes[0,2].set_ylabel('Pearson r:vel')

    axes[1,0].set_xlim([0,54])
    axes[1,0].plot(mean_rps,color='black')
    axes[1,0].fill_between(xs,mean_rps - std_rps,mean_rps + std_rps,alpha=0.3,color='gray')
    axes[1,0].axhline(0,linestyle=':',color='black')
    axes[1,0].set_ylabel('Distance:pDist (pixels)')

    axes[1,1].set_xlim([0,54])
    axes[1,1].plot(mean_ps,color='black')
    axes[1,1].fill_between(xs,mean_ps - std_ps,mean_ps + std_ps,alpha=0.3,color='gray')
    axes[1,1].axhline(0,linestyle=':',color='black')
    axes[1,1].set_ylabel('Pearson r:pDist')

    axes[1,2].set_xlim([0,54])
    axes[1,2].plot(mean_as,color='black')
    axes[1,2].fill_between(xs,mean_as - std_as,mean_as + std_as,alpha=0.3,color='gray')
    axes[1,2].axhline(0,linestyle=':',color='black')
    axes[1,2].set_ylabel('Pearson r:angle')


    fig.tight_layout()
    fig.set_size_inches([11,6.5])
    fig.savefig('fig_test.png',dpi=300)
    plt.show()
    print('Done!')
