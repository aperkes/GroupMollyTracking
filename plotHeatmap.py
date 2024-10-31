
import pandas as pd
import numpy as np
import sys,os
import cv2 

from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.collections import LineCollection

from processTracks import get_stats,get_meta,get_tracks,get_distance,clean_track


## Builds the string for a given fish
# Relies on global variables
def build_fstr(i,f):
    x,y = track_array[i,f]
    r,thet = track_polar[i,f]
    vel = velocity_array[i,f]
    heading = heading_array[i,f]
    nn_dist = nn_array[i,f]
    mean_dist = fishMeanDist_array[i,f]

    text_block = ('F' + str(f) + '\n' +  
                ','.join([str((x,y)),str((r,thet))]) + '\n' +
                ','.join([str(vel),str(heading)]) + '\n' +
                ','.join([str(mean_dist),str(nn_dist)]))
    return text_block

## Normalizes and scales array for a given min and max.
def norm_array(arr,vmin,vmax):
    arr2 = arr - np.nanmin(arr) 
    arr2 = arr / np.nanmax(arr2) * (vmax - vmin) + vmin
    return arr2
  
if __name__ == "__main__":

    csv_file = sys.argv[1]
    TEXT = False   
    #csv_file = sys.argv[1]
    if len(sys.argv) > 2:
        vid_file = sys.argv[2]
    else:
        vid_file = None
    #csv_path = './groupCSVs/' + csv_file
    csv_path = csv_file

    track_array,track_polar,(n_frames,n_fish,fishIDs) = get_tracks(csv_path)
    clean_array = clean_track(track_array)
    #track_array = clean_array
## Calculate lots of stat arrays for plotting

    print(track_array.shape)
    flat_tracks = np.reshape(track_array,[-1,2])
    
    print(flat_tracks.shape)
    flat_tracks = flat_tracks[~np.isnan(flat_tracks[:,0])]
    print(flat_tracks.shape)
    print(np.min(flat_tracks),np.max(flat_tracks))

    x_edges = np.arange(0,801,10)
    y_edges = x_edges
    #foo = np.histogram2d(flat_tracks[:,0],flat_tracks[:,1],bins=80)
    foo = np.histogram2d(flat_tracks[:,0],flat_tracks[:,1],bins=x_edges)
    import pdb;pdb.set_trace()
    fig,ax=plt.subplots()
    if vid_file is not None:
        cap = cv2.VideoCapture(vid_file)
        cap.set(cv2.CAP_PROP_POS_FRAMES,3600) ## skip initial darkness
        res,frame = cap.read()
        ax.imshow(frame)
        alpha = 0.5
    else:
        alpha = None
    bar = np.flipud(foo[0])
    #bar = np.fliplr(bar)
    #ax.imshow(np.flipud(foo[0]),vmin=0,vmax=250,alpha=alpha,cmap='viridis',extent=[0,800,0,800])
    #ax.imshow(bar,vmin=0,vmax=250,alpha=alpha,cmap='viridis',extent=[0,799,0,799])
    ax.imshow(bar,vmin=0,vmax=250,alpha=alpha,cmap='viridis')
    plt.show()

    fig.savefig('./imgs2/' + os.path.basename(csv_file).replace('.csv','png'),dpi=300)
