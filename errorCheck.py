
import pandas as pd
import numpy as np
import sys
import cv2 

from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.collections import LineCollection

from processTracks import get_stats,get_meta,get_tracks,get_distance


## Builds the string for a given fish
# Relies on global variables
def build_fstr(i,f):
    x,y = track_array[i,f]
    r,thet = np.round(track_polar[i,f],3)
    vel = np.round(velocity_array[i,f],3)
    heading = np.round(heading_array[i,f],3)
    nn_dist = np.round(nn_array[i,f],3)
    mean_dist = np.round(fishMeanDist_array[i,f],3)

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

    TEXT = True
    csv_file = sys.argv[1]

    track_array,track_polar,(n_frames,n_fish,fishIDs) = get_tracks(csv_file)
## Calculate lots of stat arrays for plotting
    meanR_array = np.nanmean(track_polar[:,:,0],axis=1)
    meanThet_array = np.nanmean(track_polar[:,:,1],axis=1)

    stdR_array = np.nanstd(track_polar[:,:,0],axis=1)
    stdThet_array = np.nanstd(track_polar[:,:,1],axis=1)
    track_stats,[velocity_array,heading_array,distance_array] = get_stats(track_array,track_polar)
    
    meanVel_array = np.nanmean(velocity_array,axis=1)
    stdVel_array = np.nanstd(velocity_array,axis=1)

    xyVel_array = np.full([n_frames,n_fish,2],np.nan)
    for n in range(n_fish):
        xyVel_array[1:,n] = np.diff(track_array[:,n],axis=0)
    stdHead_array = np.nanstd(heading_array,axis=1)
    meanDist_array = np.nanmean(distance_array,axis=(1,2))
    fishMeanDist_array = np.nanmean(distance_array,axis=1)
    #import pdb;pdb.set_trace()
    nn_array = np.nanmin(distance_array,axis=1)
    meanNn_array = np.nanmean(nn_array,axis=1)

## NOTE: start here
## Something seems off with nn-distance still
## Fish are often stationary, which feels suspicious
    first_frame = np.argmax(~np.isnan(track_array[:,0,0])) + 500
    print(first_frame)
    if len(sys.argv) > 2:
        vid_file = sys.argv[2]
        cap = cv2.VideoCapture(vid_file)
        cap.set(cv2.CAP_PROP_POS_FRAMES, first_frame-1)
        res,frame = cap.read()
        gray0 = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    else:
        vid_file = None

    #import pdb;pdb.set_trace()
    fig,ax = plt.subplots()
    ax.axis('off')
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
    #scat = ax.plot(track_array[0,:,0],track_array[0,:,1],marker='.',linestyle='None')

    if vid_file is not None:
        im=plt.imshow(gray0,cmap='gray',vmin=0,vmax=255)

    sub_xs = np.arange(n_frames - first_frame)
    sub_xs = sub_xs * 800 / len(sub_xs)
    alph = 0.2
    min_vel = str(np.round(np.nanmin(meanVel_array),3))
    max_vel = str(np.round(np.nanmax(meanVel_array),3))


    min_vstd = str(np.round(np.nanmin(stdVel_array),3))
    max_vstd = str(np.round(np.nanmax(stdVel_array),3))

    min_Tstd = str(np.round(np.nanmin(stdThet_array),3))
    max_Tstd = str(np.round(np.nanmax(stdThet_array),3))

    min_Rstd = str(np.round(np.nanmin(stdR_array),3))
    max_Rstd = str(np.round(np.nanmax(stdR_array),3))

    min_dist = str(np.round(np.nanmin(meanDist_array),3))
    max_dist = str(np.round(np.nanmax(meanDist_array),3))

    min_Ndist = str(np.round(np.nanmin(meanNn_array),3))
    max_Ndist = str(np.round(np.nanmax(meanNn_array),3))

    PLOT = False
    if PLOT:
        ax.plot(sub_xs,norm_array(meanVel_array[first_frame:],0,100),alpha=alph)
        ax.plot(sub_xs,norm_array(stdVel_array[first_frame:],100,200),alpha=alph)
        ax.plot(sub_xs,norm_array(stdThet_array[first_frame:],200,300),alpha=alph)
        ax.plot(sub_xs,norm_array(stdR_array[first_frame:],300,400),alpha=alph)
        ax.plot(sub_xs,norm_array(meanDist_array[first_frame:],400,500),alpha=alph)
        ax.plot(sub_xs,norm_array(meanNn_array[first_frame:],500,600),alpha=alph)

    vid_text = str(track_stats).replace(',','\n')
    vid_text = " " + vid_text[1:]
    vid_text = csv_file + '\n' + vid_text[:-1]
    if TEXT:
        ax.text(805,5,vid_text)
        ax.text(0,550,"MeanNnDistance: " + ",".join([min_Ndist,max_Ndist]))
        ax.text(0,450,"MeanDistance: " + ",".join([min_dist,max_dist]))
        ax.text(0,50,"MeanVelocity: " + ",".join([min_vel,max_vel]))
        ax.text(0,150,"StdVelocity: " + ",".join([min_vstd,max_vstd]))
        ax.text(0,250,"StdTheta: " + ",".join([min_Tstd,max_Tstd]))
        ax.text(0,350,"StdR: " + ",".join([min_Rstd,max_Rstd]))

        f_texts = []
        for f in range(n_fish):
            x,y = track_array[first_frame,f] + 5
            f_texts.append(ax.text(x,y,build_fstr(first_frame,f)))

    key = "Fish ID \n" \
            "(x,y); (r,theta)) \n" \
            "velocity; heading \n" \
            "mean_dist; nn_dist."

    ax.text(800,700,key)
    cors = ['red','green','blue','purple']
    scat = ax.scatter([0,0,0,0],[0,0,0,0],c=cors)
    #frame_plot, = ax.plot((0,0),(0,600),color='black',alpha=alph)
    frame_id = ax.text(5,790,"0")
    plt.xlim(0,800)
    plt.ylim(0,800)

    lines = [[(0,0),(0,0)] for n in range(n_fish)]
    vel_lc = LineCollection(lines,colors=cors,lw=4)
    #import pdb;pdb.set_trace()
    ax.add_collection(vel_lc)

        
## This could probably be defined separately, but easier here
## This animates the plot, given a bunch of global variables
    def animate(i):
        xs = track_array[i+first_frame,:,0]
        ys = track_array[i+first_frame,:,1]

        if vid_file is not None:
            cap.set(cv2.CAP_PROP_POS_FRAMES, i + first_frame -1)
            res,frame=cap.read()
            gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
            #if cv2.waitKey(5) & 0xFF == ord('p'):
            #    while 0xFF is not ord('c'):
            #        cv2.waitKey(0)
            #import pdb;pdb.set_trace()
            im.set_array(gray)

        #print(i,x,y)
        array = np.c_[xs,ys]
        scat.set_offsets(array)
        #frame_plot.set_data([(i,i),(0,600)])
        frame_id.set_text(str(i))
        lines = []
        for j in range(len(xs)):
## Get fish lines
            line0x = xs[j]
            line0y = ys[j]
            line1x = xs[j] - xyVel_array[i+first_frame,j,0]
            line1y = ys[j] - xyVel_array[i+first_frame,j,1]
            line = [(line0x,line0y),(line1x,line1y)]
            lines.append(line)

## Update fish text too
            if TEXT:
                fish_str = build_fstr(i+first_frame,j)
                fish_text = f_texts[j]
                fish_text.set_text(fish_str)
                fish_text.set_x(line0x + 5)
                fish_text.set_y(line0y + 5)
        vel_lc.set_segments(lines) ## Lines
        return 0
    def onClick(event):
        global ani_running
        if ani_running:
            ani.pause()
            ani_running = False
        else:
            ani.resume()
            ani_running = True
    #ani = animation.FuncAnimation(fig, animate, frames = n_frames-first_frame,interval=1000)
    ani_running = True
    fig.canvas.mpl_connect('button_press_event',onClick)
    ani = animation.FuncAnimation(fig, animate, frames = 200,interval=1)
    #writervideo = animation.FFMpegWriter(fps=10)
    #ani.save('test_day37.mp4', writer=writervideo)
    plt.show()

