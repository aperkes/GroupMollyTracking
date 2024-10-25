import pandas as pd
import sys,os
import numpy as np
from matplotlib import pyplot as plt
import itertools
from scipy.stats import pearsonr
import warnings
from tqdm import tqdm
from datetime import datetime

import cv2

warnings.filterwarnings('ignore')

#in_file = sys.argv[1]
#vid_file = sys.argv[2]

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

def clean_track(track_array,bins=80,thresh=1000):
    clean_array = np.array(track_array)
    a_x = clean_array[:,:,0]
    a_y = clean_array[:,:,1]
    flat_array = np.reshape(track_array,[-1,2])
    flat_array = flat_array[~np.isnan(flat_array[:,0])]

    hist,xedges,yedges = np.histogram2d(flat_array[:,0],flat_array[:,1],bins=bins)
    xedges_ = [0]
    xedges_.extend(xedges)
    yedges_ = [0]
    yedges_.extend(yedges)
    clips = np.argwhere(hist > thresh)
    for c in clips:
        x0,x1 = xedges[c[0]],xedges[c[0]+1]
        y0,y1 = yedges[c[1]],yedges[c[1]+1]
        bad_spots = np.argwhere((a_y >= y0) & (a_y <= y1) & (a_x >= x0) & (a_x <= x1))
        clean_array[bad_spots[:,0],bad_spots[:,1]] = np.nan
    #import pdb;pdb.set_trace()
    return clean_array

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
def plot_video(track_array,video_file):
    visualize = False
    if visualize:
        i = 0
        cors = [(255,0,0),(0,255,0),(0,0,255),(255,0,255)]
        cap = cv2.VideoCapture(vid_file)
        while True:
            ret, frame = cap.read()
            if not ret:
                break
            for f in range(n_fish):
                x,y = track_array[i,f].astype(int)

                cor = cors[f]
                if ~np.isnan(x):
                    cv2.circle(frame,(x,y),5,cor,4)
            cv2.imshow('overlay',frame)
            if cv2.waitKey(5) & 0xFF == ord('q'):
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
    csv_dir = sys.argv[1]
    pi_dict = {'pi11':0,'pi13':1}

    df_list = []
    long_list = []
    dist_list = [[] for p in pi_dict.keys()]
    vel_list = [[] for p in pi_dict.keys()]
    head_list = [[] for p in pi_dict.keys()]

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
            #continue

        csv_file = '/'.join([csv_dir,csv_file])
        track_array,track_polar,_ = get_tracks(csv_file)
        clean_array = clean_track(track_array)
        clean_polar = clean_track(track_polar)
        track_array,track_polar = clean_array,clean_polar
        pi,day = get_meta(csv_file,csv_dir)
        #import pdb;pdb.set_trace()
        track_stats,_ = get_stats(track_array,track_polar)
## Build hourly dicts
        binned_tracks = bin_tracks(track_array)
        binned_polars = bin_tracks(track_polar)
        for t in range(len(binned_tracks)):
            binned_stats,_ = get_stats(binned_tracks[t],binned_polars[t])

            keys = list(binned_stats.keys())
            for k in keys:
                binned_stats[k+'_'] = binned_stats[k]
                binned_stats[k] = track_stats[k]
            binned_stats['CSV'] = csv_file
            binned_stats['ExpDay'] = day
            binned_stats['Pi'] = pi
            binned_stats['Hour'] = t

            long_list.append(binned_stats)

        track_stats['CSV'] = csv_file
        track_stats['ExpDay'] = day
        track_stats['Pi'] = pi

        #import pdb;pdb.set_trace()
        df_list.append(track_stats)

        if False:
            dist_mean = track_stats['dist_mean']
            dist_std = track_stats['dist_std']
            vel_mean = track_stats['vel_mean']
            vel_std = track_stats['vel_std']
            head_mean = track_stats['head_mean']
            head_std = track_stats['head_std']
            
            dist_list[p_].append([dist_mean,dist_std])
            vel_list[p_].append([vel_mean,vel_std])
            head_list[p_].append([head_mean,head_std])

## Build df
    csv_df = pd.DataFrame(df_list)
    long_df = pd.DataFrame(long_list)
    cols = csv_df.columns.to_list()
    long_cols = long_df.columns.to_list()
## Move some columns around 
    new_cols = cols[-3:]
    new_cols.extend(cols[:-3])
    csv_df = csv_df[new_cols]
    new_cols_ = long_cols[-3:]
    new_cols_.extend(long_cols[:-3])
    long_df = long_df[new_cols_]

    if True:
        csv_df.to_csv('./JolleTracksAll_3.csv',index=False)
        long_df.to_csv('./JolleTracksHourly_3.csv',index=False)

    if False:
        fig,(ax,ax2,ax3) = plt.subplots(1,3)
        ax,dist_array = plot_lines(dist_list,ax,labels=['pi11','pi13'])
        ax2,vel_array = plot_lines(vel_list,ax2,labels=['pi11','pi13'])
        ax3,head_array = plot_lines(head_list,ax3,labels=['pi11','pi13'])

        ax.set_xlabel('Days since birth')
        ax.set_ylabel('Pairwise mean distance')

        ax2.set_xlabel('Days since birth')
        ax2.set_ylabel('Pairwise mean velocity')

        ax3.set_xlabel('Days since birth')
        ax3.set_ylabel('Mean heading correlation')

        ax.legend()
        ax2.legend()
        ax3.legend()

        plt.show()
