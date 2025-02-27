import pandas as pd
import sys,os
import numpy as np
from matplotlib import pyplot as plt
import itertools
from scipy.stats import pearsonr
from scipy import ndimage

import statsmodels.api as sm
from statsmodels.formula.api import ols,mixedlm

import warnings
from tqdm import tqdm
from datetime import datetime

from processTracks import deep_clean_track

import cv2

warnings.filterwarnings('ignore')

#in_file = sys.argv[1]
#vid_file = sys.argv[2]

MAX_VEL = 200
MAX_THETA =  np.pi/4
SPLIT_DIST = 30 ## pixels, picks confident tracks. 
## Reads in array of xy coordints (frames,2)
## spirts out array of polar coordints (frames,2)
def xy_to_polar(a,center=(400,400)):
    a_polar = np.empty(a.shape)
    dx = a[:,0] - center[0]
    dy = a[:,1] - center[1]
    a_polar[:,0] = np.sqrt(dx**2 + dy**2)
    a_polar[:,1] = np.arctan2(dy,dx)
    return a_polar

## Takes an array of xy coordinates (frames,2) 
## Returns an vector of distances of length (frames - 1)
def get_distance(a):
    a0 = a[:-1]
    a1 = a[1:]
    dist = np.linalg.norm(a1-a0,axis=1)
    return dist

## Gets the index of the next matching value (by default nan)

## Should use numba if I want this fast.
def get_next(vec,val=np.nan):
    for i in range(len(vec)):
        if np.isnan(vec[i]):
            return i
    return i+1

## read in the csv, convert it into a numpy array
## Outputs array of xy coordinates, polar coordinates, meta_info
###  Out arrays are n_frames x n_fish x 2
###  meta_info = n_frames,n_fish,fishIDs
def get_tracks(in_file):
    track_df = pd.read_csv(in_file)

    try:
        track = np.array(track_df[['x','y']])
    except:
        import pdb;pdb.set_trace()

    fishIDs = track_df.id.unique()
    n_fish = len(fishIDs)
    n_frames = max(track_df.frame) + 1
    track_array = np.full([n_frames,n_fish,2],np.nan)
    track_polar = np.array(track_array)

    for f in range(n_fish):
        f_df = track_df[track_df.id == fishIDs[f]]
        
        xy_f = np.array(f_df[['x','y']]).astype(float)
        xy_f[xy_f == -1] = np.nan

        indices = np.array(f_df.frame)
        track_array[indices,f] = xy_f

        
        track_polar[indices,f] = xy_to_polar(xy_f)


    #track_array[track_array == -1] = np.nan
    return track_array,track_polar,[n_frames,n_fish,fishIDs]

def bin_tracks(track_array,bin_size=3600):
    n_arrays = int(track_array.shape[0] / 3600) ## note, this will drop possibe trailing incomplete hours
    sub_tracks = []
    for n in range(n_arrays):
        i0 = n*bin_size
        i1 = (n+1)*bin_size
        sub_tracks.append(track_array[i0:i1])
    return sub_tracks


##### METRICS AND STUFF #####
""" 
We will be calculating a few metrics
0. Overall space use
1. interfish distance
2. Distance to the wall (which is maybe already in the csv?)
    if not, I can get distance from center
3. Motion
    a) Angular Speed 
    b) Angular Direction
    c) normal speed

"""
def get_histogram(track_polar,plot_me = False):
    hists = []
    if plot_me:
        fig,(ax,ax1) = plt.subplots(1,2)
    for f in range(n_fish):
        subtrack = track_polar[:,f]
        subtrack = subtrack[~np.isnan(subtrack[:,0])]
        subtrack = subtrack[~np.isnan(subtrack[:,1])]
        hist_theta,edges = np.histogram(subtrack[:,1])
        hist_r,edges_ = np.histogram(subtrack[:,0])
        #hists.append(r)
        r_width = edges_[1] - edges_[0]
        theta_width = edges[1] - edges[0]
        if plot_me:
            ax.bar(edges_[:-1],hist_r,alpha=0.5,width=r_width)
            ax1.bar(edges[:-1],hist_theta,alpha=0.5,width=theta_width)

    if plot_me:
        plt.show()
    return hists

def clean_vel(track_array,cleanup=True):
    velocity_array = np.full([n_frames,n_fish],np.nan)
    for n in range(n_fish):
        n_vel = get_distance(track_array[:,n])
        if cleanup:
            n_vel[n_vel > MAX_VEL] = np.nan
            n_vel[n_vel == 0] = np.nan
        velocity_array[1:,n] = n_vel
    smooth_velocity = ndimage.gaussian_filter1d(velocity_array,5,0,radius=1)
    return velocity_array,smooth_velocity

def clean_dist(track_array):
    distance_array = np.full([n_frames,n_fish,n_fish],np.nan)
    for i,j in itertools.combinations(np.arange(n_fish),2):
## Calculate correlation of distance to wall
        xs = track_polar[:,i,0]
        ys = track_polar[:,j,0]
        good_indices = (~np.isnan(xs)) & (~np.isnan(ys))

### Calculate mean distance between fish
        track_i = track_array[:,i]
        track_j = track_array[:,j]
        good_indices = (~np.isnan(track_i[:,0])) & (~np.isnan(track_j[:,0]))
        distance_array[good_indices,i,j] = np.linalg.norm(track_i[good_indices] - track_j[good_indices],axis=1)
        distance_array[good_indices,j,i] = distance_array[good_indices,i,j]
    return distance_array

## FUNCTION to calculate lots of motion statistics
def get_stats(track_array,track_polar):
    n_frames,n_fish,_ = track_array.shape
    n_stats = 5
    stat_array = np.zeros([n_stats,n_fish,n_fish])

    velocity_array = np.full([n_frames,n_fish],np.nan)
    angMom_array = np.array(velocity_array)
    distance_array = np.full([n_frames,n_fish,n_fish],np.nan)

    MAX_VEL = 200
    MAX_THETA =  np.pi/4

## I might be able to speed these up:
    for n in range(n_fish):
        velocity_array[1:,n] = get_distance(track_array[:,n])
        angMom_array[1:,n] = np.diff(track_polar[:,n,1])
## Clean up tracking errors
    velocity_array[velocity_array > MAX_VEL] = np.nan
    velocity_array[velocity_array == 0] = np.nan

    angMom_array[np.abs(angMom_array) > MAX_THETA] = np.nan
    angMom_array[angMom_array == 0] = np.nan

    diff_array = np.diff(track_array,axis=0,prepend=np.nan)
    diff_array[np.isnan(velocity_array)] = np.nan
    #diff_array[velocity_array > MAX_VEL] = np.nan

    angle_array = np.arctan2(diff_array[:,:,0],diff_array[:,:,1])/np.pi*180

    #angle_array[angle_array == 0] = np.nan ##solves arctan(0,0)
    angMom_array_N = np.array(angMom_array) / (2*np.pi)
    angle_array_N = angle_array / 180
    angle_medians = np.nanmedian(angle_array_N,axis=1)
    angMom_medians = np.nanmedian(angMom_array_N,axis=1)
    for n in range(n_fish):
        angMom_array_N[:,n] = 1 - np.abs(angMom_array[:,n] - angMom_medians)
        angle_array_N[:,n] = 1 - np.abs(angle_array_N[:,n] - angle_medians)

        angle_array_N[:,n][np.isnan(velocity_array[:,n])] = np.nan
        angMom_array_N[:,n][np.isnan(angMom_array[:,n])] = np.nan
    
    polarity_array = np.nanmean(angle_array_N,axis=1)
    polarity_array[polarity_array == 1] = np.nan

    rotation_array = np.nanmean(angMom_array_N,axis=1)
    rotation_array[rotation_array == 1] = np.nan

    all_good_indices = np.zeros(n_frames)
    for i,j in itertools.combinations(np.arange(n_fish),2):
## Calculate correlation of distance to wall
        xs = track_polar[:,i,0]
        ys = track_polar[:,j,0]
        good_indices = (~np.isnan(xs)) & (~np.isnan(ys))
        if sum(good_indices) > 10:
            xs = xs[good_indices]
            ys = ys[good_indices]
            r_stat,p = pearsonr(xs,ys)
        else:
            r_stat = np.nan
        #print(i,j,r_stat,p)
        stat_array[0,i,j] = r_stat

### Calculate mean distance between fish
        track_i = track_array[:,i]
        track_j = track_array[:,j]
        good_indices = (~np.isnan(track_i[:,0])) & (~np.isnan(track_j[:,0]))
        distance_array[good_indices,i,j] = np.linalg.norm(track_i[good_indices] - track_j[good_indices],axis=1)
        distance_array[good_indices,j,i] = distance_array[good_indices,i,j]
        all_good_indices[good_indices] = True
        mean_distance = np.nanmean(np.linalg.norm(track_i[good_indices] - track_j[good_indices],axis=1))
        mean_distance = np.nanmedian(np.linalg.norm(track_i[good_indices] - track_j[good_indices],axis=1))
        stat_array[1,i,j] = mean_distance

### Calculate speed correlation
        #vel_i = get_distance(track_array[:,i])
        #vel_j = get_distance(track_array[:,j])

## Will probably need to add smoothing, in addition to cutting out unreasonable values.
        #vel_i[vel_i > 200] = np.nan
        #vel_j[vel_j > 200] = np.nan
        vel_i = velocity_array[:,i]
        vel_j = velocity_array[:,j]

        good_indices = (~np.isnan(vel_i)) & (~np.isnan(vel_j))
        if sum(good_indices) > 100:
            r_stat,p = pearsonr(vel_i[good_indices],vel_j[good_indices])
            #good_indices1 = good_indices[1:]
            r_stat2,p2 = pearsonr(angle_array[good_indices,i],angle_array[good_indices,j])
        else:
            r_stat,p = np.nan,np.nan
            r_stat2,p2 = np.nan,np.nan
        stat_array[2,i,j] = r_stat
        stat_array[4,i,j] = r_stat2
## Get heading correlations
        #head_i = np.diff(track_polar[:,i,1])
        #head_j = np.diff(track_polar[:,j,1])
        #head_i[np.abs(head_i) > np.pi/8] = np.nan
        #head_j[np.abs(head_j) > np.pi/8] = np.nan
        head_i = angMom_array[:,i]
        head_j = angMom_array[:,j]

        good_indices = ~np.isnan(head_i) & ~np.isnan(head_j)
        if sum(good_indices) > 10:
            r_stat,p = pearsonr(head_i[good_indices],head_j[good_indices])
        else:
            r_stat = np.nan
        stat_array[3,i,j] = r_stat
    if False:
        print('dist from center:')
        print(stat_array[0]) ## This is cohesion
        print('Inter-fish Distance')
        print(stat_array[1]) ## This is mean distance
        print('velocity r')
        print(stat_array[2]) ## Velocity cohesion 
        print('Angular Mom. r')
        print(stat_array[3]) ## Heading Cohesion
        print('Angle r')
        print(stat_array[4]) ## Angle Cohesion
## What do I actually want here? 
## Speed correlation (and std)
## Mean distance (and std) STD is across fish


    #distance_array = distance_array[all_good_indices == True] ## NOTE: this is bad, do I need it?
    dist_mean = np.nanmean(distance_array) ## This is just mean distance
    dist_std = np.nanmean(np.nanstd(distance_array,axis=(1,2)))

    nn_dist = np.nanmin(distance_array,axis=1) ## Mean nearest neighbord distance
    nn_mean = np.nanmean(nn_dist)
    nn_std = np.nanstd(nn_dist)

    pDist_mean = np.nanmean(track_polar[:,:,0])
    pDist_std = np.nanstd(track_polar[:,:,0])
    pDistC_mean = np.nanmean(stat_array[0])
    pDistC_std = np.nanstd(stat_array[0])

    vel_mean = np.nanmean(velocity_array)
    vel_std = np.nanstd(np.nanstd(velocity_array,axis=1))
    angMomC_mean = np.nanmean(stat_array[3]) ## this is the mean cohesion
    angMomC_std = np.nanstd(stat_array[3])
    velC_mean = np.nanmean(stat_array[2])
    velC_std = np.nanstd(stat_array[2])

    angleC_mean = np.nanmean(stat_array[4])
    angleC_std = np.nanstd(stat_array[4])

    polarity_mean = np.nanmean(polarity_array)
    polarity_std = np.nanstd(polarity_array)

    rotation_mean = np.nanmean(rotation_array)
    rotation_std = np.nanstd(rotation_array)

    stat_arrays = [velocity_array,angMom_array,distance_array]
    track_stats = {
        'dist_mean':dist_mean,
        'dist_std':dist_std,
        'vel_mean':vel_mean,
        'vel_std':vel_std,
        'velC_mean':velC_mean,
        'velC_std':velC_std,
        'angMC_mean':angMomC_mean,
        'angMC_std':angMomC_std,
        'polarity_mean':polarity_mean,
        'polarity_std':polarity_std,
        'rotation_mean':rotation_mean,
        'rotation_std':rotation_std,
        'pDist_mean':pDist_mean,
        'pDist_std':pDist_std,
        'pDistC_mean':pDistC_mean,
        'pDistC_std':pDistC_std,
        'NearN_mean':nn_mean,
        'NearN_std':nn_std,
        'angleC_mean':angleC_mean,
        'angleC_std':angleC_std
    }
    return track_stats,stat_arrays

### Function that takes a list of lists stats for two pis (2:n_vids:mean/std
### Returns an array of that list, along with an axis with a the stats plotted

def plot_lines(stat_list,ax=None,labels=[None,None]):
    max_days = max([len(stat_list[0]),len(stat_list[1])])
    stat_array = np.full([2,max_days,2],np.nan)
    if ax is None:
        fig,ax = plt.subplots()
    for i in range(2):
        stat_array[i,:len(stat_list[i])] = stat_list[i]

    xs_i = np.arange(len(stat_list[0]))
    xs_j = np.arange(len(stat_list[1]))

    ys_i = stat_array[0,xs_i,0]
    std_i = stat_array[0,xs_i,1]
    ys_j = stat_array[1,xs_j,0]
    std_j = stat_array[1,xs_j,1]

    ax.plot(xs_i,ys_i,label=labels[0])
    ax.fill_between(xs_i,ys_i - std_i,ys_i + std_i,color='gray',alpha=0.2)

    ax.plot(xs_j,ys_j,label=labels[1])
    ax.fill_between(xs_j,ys_j - std_j,ys_j + std_j,color='gray',alpha=0.5)
    return ax,stat_array


## optionally overlay the tracks on a video to show how they work
def plot_video(track_array,vid_file,viz=True,trail=0):
    if viz:
        i = 0
        cors = [(255,0,0),(0,255,0),(0,0,255),(255,0,255)]
        cap = cv2.VideoCapture(vid_file)
        n_fish = np.shape(track_array)[1]
        n_points = trail 
        wt = 1
        while True:
            ret, frame = cap.read()
            if not ret:
                break
            for f in range(n_fish):
                if np.isnan(track_array[i,f,0]):
                    continue
                x,y = track_array[i,f].astype(int)
                cor = cors[f]
                cv2.circle(frame,(x,y),9,cor,-1)
                for n in range(1,n_points+1):
                    if np.isnan(track_array[i-n,f,0]):
                        break
                    x0,y0 = track_array[i-n,f].astype(int)
                    x1,y1 = track_array[i-n+1,f].astype(int)
                    cv2.circle(frame,(x0,y0),max(9-n,1),cor,-1)
                    cv2.line(frame,(x0,y0),(x1,y1),cor,2)

            cv2.imshow('overlay',frame)
            if cv2.waitKey(1) & 0xFF == ord('v'):
                wt = 500
            if cv2.waitKey(wt) & 0xFF == ord('q'):
                break
            i += 1
        cv2.destroyAllWindows()
        cap.release()

## Extract exp day and pi from video
def get_meta(csv_file,csv_dir):
    pi = csv_file.split('_')[1]
    ds = csv_file.split('_')[-1].replace('.csv','')
    all_pi_csvs = []
    for c in os.listdir(csv_dir):
        if pi in c:
            all_pi_csvs.append(c)
    all_pi_csvs = sorted(all_pi_csvs)
    first_ds = all_pi_csvs[0].split('_')[-1].replace('.csv','')
    first_datetime = datetime.strptime(first_ds, "%y%m%d")
    csv_datetime = datetime.strptime(ds, "%y%m%d")
    delta_t = csv_datetime - first_datetime
    delta_days = delta_t.days

    return pi,delta_days

if __name__ == '__main__':
    out_file = open('test_rpt2.csv','w')
    csv_dir = sys.argv[1]
    pi_dict = {'pi13':0,'pi14':1} #,'pi54':2}

    df_list = []
    long_list = []
    dist_list = [[] for p in pi_dict.keys()]
    vel_list = [[] for p in pi_dict.keys()]
    head_list = [[] for p in pi_dict.keys()]

    median_arrays = {}
    median_arrays2 = {}
    median_arrays3 = {}

    for csv_file in tqdm(sorted(os.listdir(csv_dir))):
        if 'GR' in csv_file:
            continue
        in_list = False
        for p in pi_dict.keys():
            if p in csv_file:
                p_ = pi_dict[p]
                in_list = True
        if not in_list:
            pass
            continue

        csv_file = '/'.join([csv_dir,csv_file])
        track_array,track_polar,[n_frames,n_fish,fishIDs] = get_tracks(csv_file)
        clean_array,velocity_array,distance_array = deep_clean_track(track_array)
        track_polar[np.isnan(clean_array)] = np.nan
        track_array = clean_array
        pi,day = get_meta(csv_file,csv_dir)
        if pi not in median_arrays.keys():
            median_arrays[pi] = []
            median_arrays2[pi] = []
            median_arrays3[pi] = []

        pDist_array = track_polar[:,:,0]
        #distance_array = clean_dist(track_array)
        #velocity_array,smooth_vel = clean_vel(track_array)
        if False:
            smooth_vel = ndimage.gaussian_filter1d(velocity_array,5,0,radius=1)
        else:
            smooth_vel = velocity_array
        velocity_confident = np.array(velocity_array)
        min_dist = np.nanmin(distance_array,axis=1)

        #velocity_confident[min_dist < SPLIT_DIST]
        #distanceM_array = np.nanmean(distance_array,axis=1)
        distanceM_array = np.nanmedian(distance_array,axis=1)

## Split all these tracks based on where we are confident with veloicty
        distance_confident = np.array(distanceM_array)
        distance_confident[np.isnan(velocity_array)] = np.nan

        pDist_confident = np.array(pDist_array)
        pDist_confident[np.isnan(velocity_array)] = np.nan

        if False:
            velocity_array = smooth_vel
        smooth_confident = np.array(smooth_vel)
        smooth_confident[np.isnan(velocity_array)] = np.nan

        if False:
            velocity_array = smooth_confident

        ranked_velocity = np.argsort(velocity_array,axis=1)
        ranked_pDist = np.argsort(pDist_array,axis=1)
        ranked_dist = np.argsort(distanceM_array,axis=1)

        drop_count = np.sum(np.isnan(velocity_array),axis=1)
        #drop_count = np.sum(np.isnan(velocity_confident),axis=1)

## Get unique spots where all 4 tracks are present
        good_indices = np.arange(len(drop_count))[drop_count == 0]
        diff_indices = np.diff(good_indices)
        good_indices = good_indices[1:][diff_indices > 1]
        
        #velocity_array = velocity_confident ## Definitely some question on what I should do this on. 
        ranked_array = np.full(velocity_array.shape,np.nan)
        ranked_pDistarray = np.full(velocity_array.shape,np.nan)
        ranked_distarray = np.full(velocity_array.shape,np.nan)
## This is a tricky little bit
## We're using these good indices, working backwards and forwards to define tracklets
## Then finding the means, ranking the tracklets, and assigning them to the ranked array.
        if True:
            len_list = []
            means = [0] * 4
            means_dist = [0] * 4
            means_Pdist = [0] * 4
            n_indices = np.zeros([n_fish,2]).astype(int)
            for i in good_indices:
                i_list = []
                n_list = []
                v_list = []
                next_nans = []
                last_nans = []
                for n in range(n_fish):
                    next_nan = get_next(velocity_array[i:,n])
                    next_nans.append(next_nan)
                    last_nan = get_next(velocity_array[i::-1,n])
                    last_nans.append(last_nan)
                
                next_nan = min(next_nans)
                last_nan = min(last_nans)
                if next_nan - last_nan < 10:
                    continue
                len_list.append(next_nan - last_nan)
                for n in range(n_fish):
                    #next_nan = get_next(velocity_array[i:,n])
                    #last_nan = get_next(velocity_array[i::-1,n])
                    tracklet_n = velocity_array[i-last_nan+1:i+next_nan,n]
                    if np.sum(np.isnan(tracklet_n)) > 1:
                        import pdb;pdb.set_trace()
                    tracklet_dist = distanceM_array[i-last_nan+1:i+next_nan,n]
                    tracklet_Pdist = pDist_array[i-last_nan+1:i+next_nan,n]
                    i_list.extend(np.arange(i-last_nan+1,i+next_nan))
                    n_list.extend([n] * len(tracklet_n))
                    v_list.extend(tracklet_n)                    

                    means[n] = np.nanmean(tracklet_n)
                    means_dist[n] = np.nanmean(tracklet_dist)
                    means_Pdist[n] = np.nanmean(tracklet_Pdist)

                    n_indices[n,0] = i-last_nan + 1
                    n_indices[n,1] = i+next_nan
## Calculate repeatability for this tracklet
                #print(pi,i)
## Store repeatability
                #import pdb;pdb.set_trace()
                tracklet_rank = np.argsort(means)
                tracklet_rank_dist = np.argsort(means_dist)
                tracklet_rank_Pdist = np.argsort(means_Pdist)
                #tracklet_rank = tracklet_rank_dist
## Do you want to sort by one of them, or by all of them? 
                if len(i_list) <= 10:
                    #print('##### SKIPPING ####',i)
                    continue
                for n_ in range(n_fish):
                    n_v = tracklet_rank[n_]
                    n_c = tracklet_rank_dist[n_]
                    n_p = tracklet_rank_Pdist[n_]
                    if False:
                        n_c = n_v
                        n_p = n_v
                    #n = tracklet_rank[n_]

                    i0_v,i1_v = n_indices[n_v]
                    i0_c,i1_c = n_indices[n_c]
                    i0_p,i1_p = n_indices[n_p]

                    ranked_array[i0_v:i1_v,n_] = velocity_array[i0_v:i1_v,n_v] ## This is lowest to highest
                    ranked_pDistarray[i0_p:i1_p,n_] = pDist_array[i0_p:i1_p,n_p] ## This is lowest to highest
                    ranked_distarray[i0_c:i1_c,n_] = distanceM_array[i0_c:i1_c,n_c] ## This is lowest to highest

                if True:
                    tmp_df = pd.DataFrame(zip(*[i_list,n_list,v_list]),columns=['i','id','vel'])
                    #cw_lm=ols('vel ~ i + C(id)',data=tmp_df).fit()
                    #import pdb;pdb.set_trace()
                    model=mixedlm('vel ~ i',data=tmp_df,groups=tmp_df["id"]).fit()
                    var_w = model.scale
                    var_a = model.cov_re.iloc[0,0]
                    #res = sm.stats.anova_lm(cw_lm,typ=2)
                    #var_w = res['sum_sq']['C(id)']
                    #var_a = res['sum_sq']['Residual']
                    rpt = var_a / (var_a+var_w)
                    rpt_output = [pi,str(day),str(i),str(rpt),str(var_w),str(var_a),csv_file,'\n']
                    out_file.write(','.join(rpt_output)) 
            #print(np.unique(len_list,return_counts=True))
        else:  ## alternatively, you can just go frame by frame, but the above is more reliable and more data.
            for i in range(4):
                ranked_array[:,0] = velocity_array[np.arange(n_frames),ranked_velocity[:,0]]
                ranked_array[:,1] = velocity_array[np.arange(n_frames),ranked_velocity[:,1]]
                ranked_array[:,2] = velocity_array[np.arange(n_frames),ranked_velocity[:,2]]
                ranked_array[:,3] = velocity_array[np.arange(n_frames),ranked_velocity[:,3]]
            ranked_array[drop_count > 0] = np.nan
## Quick and dirty. I might need to speed this up if it's too slow. 
        if False:
            fig,ax = plt.subplots()
            ax.hist(ranked_array[:,0],alpha=0.5)
            ax.hist(ranked_array[:,1],alpha=0.5)
            ax.hist(ranked_array[:,2],alpha=0.5)
            ax.hist(ranked_array[:,3],alpha=0.5)
            plt.show()

        median_speeds = np.nanmedian(ranked_array,axis=0)
        median_pdists = np.nanmedian(ranked_pDistarray,axis=0)
        median_dists = np.nanmedian(ranked_distarray,axis=0)

        median_arrays[pi].append(median_speeds)
        median_arrays2[pi].append(median_pdists)
        median_arrays3[pi].append(median_dists)
        
        if False:
            fig,ax = plt.subplots()
            ax.plot(velocity_array[:,0],color='black',alpha=0.5)
            ax.plot(smooth_vel[:,0],color='red',linestyle=':')
            plt.show()


    import pdb;pdb.set_trace()
    out_file.close()

    #pi_speeds11 = median_arrays['pi11']
    #pi_speeds13 = median_arrays['pi13']
    n_plots = len(median_arrays.keys())
    if True:
        fig,axes = plt.subplots(n_plots,3,sharex=True)
        for k_ in range(len(median_arrays.keys())):
            k = list(median_arrays.keys())[k_]
            xs = range(len(median_arrays[k]))
            axes[k_,0].plot(xs,median_arrays[k])
            axes[k_,1].plot(xs,median_arrays2[k])
            axes[k_,2].plot(xs,median_arrays3[k])
        axes[0,0].set_xlim([0,55])
        axes[0,0].set_xlabel('Days since birth')
        axes[0,0].set_ylabel('Speed')
        axes[0,1].set_ylabel('Center Distance')
        axes[0,2].set_ylabel('Median IID by fish')
        axes[1,0].set_ylabel('Speed')
        axes[1,1].set_ylabel('Center Distance')
        axes[1,2].set_ylabel('Median IID by fish')
    else:
        n_pis = len(median_arrays.keys())
        all_std_v = np.full([n_pis,100],np.nan)
        all_std_p = np.full([n_pis,100],np.nan)
        all_std_c = np.full([n_pis,100],np.nan)
        fig,(ax0,ax1,ax2) = plt.subplots(1,3,sharex=True)
        for k_ in range(len(median_arrays.keys())):
            k = list(median_arrays.keys())[k_]
            xs = np.arange(len(median_arrays[k]))
            pi_std_v = np.nanstd(np.array(median_arrays[k]),axis=1)
            pi_std_p = np.nanstd(np.array(median_arrays2[k]),axis=1)
            pi_std_c = np.nanstd(np.array(median_arrays3[k]),axis=1)
            all_std_v[k_,:len(xs)] = pi_std_v
            all_std_p[k_,:len(xs)] = pi_std_p
            all_std_c[k_,:len(xs)] = pi_std_c

            ax0.plot(pi_std_v,color='black',alpha=0.1)
            ax1.plot(pi_std_p,color='black',alpha=0.1)
            ax2.plot(pi_std_c,color='black',alpha=0.1)

        ax0.plot(np.nanmean(all_std_v,axis=0),color='black')
        ax1.plot(np.nanmean(all_std_p,axis=0),color='black')
        ax2.plot(np.nanmean(all_std_c,axis=0),color='black')

        ax0.set_xlabel('Days since birth')
        ax1.set_xlabel('Days since birth')
        ax2.set_xlabel('Days since birth')

        ax0.set_ylabel('Std Speed')
        ax1.set_ylabel('Std Center Distance')
        ax2.set_ylabel('Std IID')
    plt.show()
