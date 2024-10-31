
import cv2

from scipy import stats
import numpy as np

from matplotlib import pyplot as plt
import sys,os

from processTracks import get_tracks,clean_track

### This code will quantify algae and compare to space use

### Quantifying algae:
# Take the mode of the background for every video (ugh) 
# Bin the background, excluding the edges 
def get_algae_locs(vid_file,mask = None,down_ratio=0.1,n_frames=100):
    cap = cv2.VideoCapture(vid_file)
    total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    if total_frames < n_frames:
        #import pdb;pdb.set_trace()
        raise Exception("Vid file too short")
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

    br_masked = np.array(br)
    #import pdb;pdb.set_trace()
    if mask is not None:
        #print('masking...')
        br_masked[~mask] = np.nan

    algae_flat = br_masked.T.flatten() ## need to get it into xy coord space
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

## opencv uses [row,column] (i.e., y,x). I forget this once per year. 
def get_pixIDs(w=800,h=800,down_ratio=0.1):
    if w <= 1000:
        y_scale = int(1000 * down_ratio)
    else:
        y_scale = int(w * down_ratio)
    
    w = int(w * down_ratio)
    h = int(h * down_ratio)

    id_array = np.empty([w,h])
    xy_array = np.empty([w,h,2])
    yx_array = np.empty([w,h,2])

    for x in np.arange(w): ## There's almost certainly a more clever way to do this
        for y in np.arange(h):
            id_array[x,y] = x + y*y_scale
            xy_array[x,y] = (x,y)
            yx_array[x,y] = (y,x)
    return id_array,xy_array,yx_array

if __name__ == '__main__':
    
## This could be arg parsed...
    try:
        vid_path = sys.argv[1]
        csv_path = sys.argv[2]
        out_csv = sys.argv[3]
    except:
        print("needs 3 arguments, vid path, csv path, and out_csv")
    mask = get_mask()
    maskT = mask.T ## opencv imgs are in y,x!!!
    id_array,xy_array,yx_array = get_pixIDs()
    flat_ids = id_array[maskT != 0].flatten()
    flat_xys = xy_array[maskT != 0].reshape([-1,2])
    #mask=None

    algae_flat,br = get_algae_locs(vid_path,mask)
    flat_counts,counts = get_fish_locs(csv_path,mask)
    
    basename = os.path.basename(vid_path)
    basename = basename.replace('.mp4','')
    pi,day = basename.split('_')

    if True:
        with open(out_csv,'a') as f:
            for i in range(len(algae_flat)):
                x,y = flat_xys[i]
                line_list = [flat_ids[i],x,y,pi,day,algae_flat[i],flat_counts[i]]
                line = ','.join([str(v) for v in line_list]) + '\n'
                #print(line)
                f.write(line)
    if True:
        print(stats.pearsonr(algae_flat,flat_counts))
        fig,ax=plt.subplots()

        ax.imshow(br.T,cmap='gray')
        ax.imshow(counts,alpha=0.5,vmax=100)
        #cv2.imwrite('./br.png',br)
        plt.savefig('./imgs3/' + basename + '.png',bbox_inches='tight',dpi=300)

    print('done!')
## At some point I need to move to R probably

# Binwise correlation of algae to space use

# Calculate repeatability of algae (sliding scale?) 
