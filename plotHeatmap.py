
import pandas as pd
import numpy as np
import sys,os
import cv2 

from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.collections import LineCollection
from matplotlib.widgets import Slider

from processTracks import get_stats,get_meta,get_tracks,get_distance,clean_track, deep_clean_track


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
    #import pdb;pdb.set_trace()
    #clean_array = clean_track(track_array)
    clean_array,velocity_array,distance_array = deep_clean_track(track_array)
    track_array = clean_array
## Calculate lots of stat arrays for plotting

    print(track_array.shape)

    n_tracks = track_array.shape[1]
    window_size = 20
    half_win = window_size // 2 
    n_frames = track_array.shape[0]

    fig,ax = plt.subplots()
    lines = []
    points = []
    from matplotlib import colormaps

    cmap = colormaps.get_cmap('tab10')
    cors = [cmap(i) for i in range(n_tracks)] 
    for i in range(n_tracks):

        line, = ax.plot([],[],label=f'Trajectory {i+1}',c=cors[i])
        point, = ax.plot([],[],'o',markersize=8,c=cors[i])
        points.append(point)
        lines.append(line)

    ax.set_xlim([0,800])
    ax.set_ylim([0,800])
    
    ax_slider = fig.add_axes([0.25,0.1,0.65,0.03])
    slider = Slider(ax=ax_slider,label='Frame Index',valmin=0,valmax = n_frames,valinit=100,valstep=1)

    def update(frame_index):
        start = max(0,frame_index - half_win)
        end = min(n_frames,frame_index + half_win+1)

        for i,line in enumerate(lines):
            line.set_data(clean_array[start:end,i,0],clean_array[start:end,i,1])
            points[i].set_data(clean_array[end-1,i,0],clean_array[end-1,i,1])

        fig.canvas.draw_idle()
        return lines

    slider.on_changed(update)
    plt.axis('off')
    fig.set_size_inches(6,6)
    #plt.show()

    flat_tracks = np.reshape(track_array,[-1,2])
    
    print(flat_tracks.shape)
    flat_tracks = flat_tracks[~np.isnan(flat_tracks[:,0])]
    print(flat_tracks.shape)
    print(np.min(flat_tracks),np.max(flat_tracks))

    x_edges = np.arange(0,801,10)
    y_edges = x_edges
    #foo = np.histogram2d(flat_tracks[:,0],flat_tracks[:,1],bins=80)
    foo = np.histogram2d(flat_tracks[:,0],flat_tracks[:,1],bins=x_edges)
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
    print('bar sum:',np.nansum(bar))
    low_pass = np.nansum(bar) / 10000
    low_pass = 0
    print(low_pass)
    #bar[bar <= low_pass] = np.nan
    #bar = np.fliplr(bar)
    #ax.imshow(np.flipud(foo[0]),vmin=0,vmax=250,alpha=alpha,cmap='viridis',extent=[0,800,0,800])
    #ax.imshow(bar,vmin=0,vmax=250,alpha=alpha,cmap='viridis',extent=[0,799,0,799])
    cmap = 'viridis'
    #cmap = 'Oranges'
    ax.imshow(bar,vmin=None,vmax=150,alpha=alpha,cmap=cmap)
    plt.show()

    fig.savefig('./figs/' + os.path.basename(csv_file).replace('.csv','.png'),dpi=300)
    fig.savefig('./figs/' + os.path.basename(csv_file).replace('.csv','.svg'))
