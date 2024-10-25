
import cv2

from scipy import stats
import numpy as np

from matplotlib import pyplot as plt
import sys

from processTracks import get_tracks,clean_track

### This code will quantify algae and compare to space use

### Quantifying algae:
# Take the mode of the background for every video (ugh) 
# Bin the background, excluding the edges 
def get_algae_locs(vid_file,mask_file = None,down_ratio=0.1,n_frames=100):
    cap = cv2.VideoCapture(vid_file)
    total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    n_skip = total_frames / n_frames

    width  = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

    f_loc = 0
    count = 0
    #all_frames = np.empty([n_frames,width,height])
    small_frames = np.empty([n_frames,int(width*down_ratio),int(height*down_ratio)])

    print('reading video...')
    while cap.isOpened():
        #print(count,f_loc)
        ret,frame = cap.read()
        if not ret:
            break
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        small_gray = cv2.resize(gray,(0,0),
                    fx=down_ratio,fy=down_ratio,
                    interpolation=cv2.INTER_NEAREST)
        #all_frames[count] = gray
        small_frames[count] = small_gray

        cap.set(cv2.CAP_PROP_POS_FRAMES,f_loc + n_skip)
        f_loc += n_skip
        count += 1

    print('calculating mode')

    br,_ = stats.mode(small_frames,axis=0)

    if mask is not None:
        print('masking...')
        br[~mask] = np.nan

    algae_flat = br.flatten()
    algae_flat = algae_flat[~np.isnan(algae_flat)]
    return algae_flat,br

## Grab fish locations
def get_fish_locs(csv_path,mask=None,down_ratio=0.1):
    track_array,track_polar,(n_frames,n_fish,fishIDs) = get_tracks(csv_path)
    clean_array = clean_track(track_array)
    flat_tracks = np.reshape(track_array,[-1,2])
    flat_tracks = flat_tracks[~np.isnan(flat_tracks[:,0])]

    x_edges = np.arange(0,801,int(1/down_ratio))
    y_edges = x_edges
    counts,x_edges,y_edges = np.histogram2d(flat_tracks[:,0],flat_tracks[:,1],bins=x_edges)

    if mask is not None:
        counts[~mask] = np.nan
    flat_counts = counts.flatten()
    flat_counts = flat_counts[~np.isnan(flat_counts)]
    return flat_counts,counts

def get_mask(mask_path='./br-mask.png',down_ratio=0.1):
    mask = cv2.imread(mask_path)
    mask = cv2.resize(mask,(0,0),
                    fx=down_ratio,fy=down_ratio,
                    interpolation=cv2.INTER_NEAREST)
    mask = mask[:,:,0]
    mask = (mask / 255).astype(bool)

    return mask

if __name__ == '__main__':
    
    vid_path = sys.argv[1]
    csv_path = sys.argv[2]
    mask = get_mask()
    algae_flat,br = get_algae_locs(vid_path,mask)
    flat_counts,counts = get_fish_locs(csv_path,mask)
    print(stats.pearsonr(algae_flat,flat_counts))

    fig,ax=plt.subplots()

    ax.imshow(br,cmap='gray')

#cv2.imwrite('./br.png',br)
    plt.show()


    import pdb;pdb.set_trace()

## At some point I need to move to R probably

# Binwise correlation of algae to space use

# Calculate repeatability of algae (sliding scale?) 
