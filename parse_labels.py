import sleap
import numpy as np
from matplotlib import pyplot as plt
#from matplotlib import colormaps
import scipy.stats as stats
import pandas as pd
import os

from statsmodels.formula.api import ols
import statsmodels.api as sm

#anns = sleap.load_file('/home/ammon/Downloads/pi11.Labeled/groupTracks.pi11.slp')
ann_dir = './labeledTracks/'

n_anns = len(os.listdir(ann_dir))

mega_df_list = []
fig,(ax1,ax2) = plt.subplots(2)
cmap = plt.get_cmap('viridis')

for a_ in range(n_anns):
    a = os.listdir(ann_dir)[a_]
    #print(anns)
    print(a_,a)
    ann_file = ann_dir + a
    anns = sleap.load_file(ann_dir + a)
    #import pdb;pdb.set_trace()
    pi = a.split('.')[1]

    n_vids = len(anns.videos)
    video_anns = [[] for v in range(n_vids)]

#vid_list = []
    vid_dict = {}
    vid_counts = {}
    vid_indices = {}
## Each video needs an array, with each fish having a little track
    all_tracks = np.full([n_vids,4,5],np.nan)
    count = 0
    for v in anns.videos:
        #vid_list.append(v.filename)
        vid_indices[v.filename] = count
        count += 1
        vid_dict[v.filename] = np.full([4,5],np.nan) 
        vid_counts[v.filename] = [0 for i in range(4)]

    for l in anns:
        if len(l.instances) == 0:
           continue
        for i in l.instances:

            vid = l.video.filename
            if i.track is None:
                continue
            fish_index = int(i.track.name) - 1
            count = vid_counts[vid][fish_index]
            if count > 4:
                #import pdb;pdb.set_trace()
                continue
            try:
                head_xy = np.array([i['Head'].x,i['Head'].y])
                body_xy = np.array([i['Body'].x,i['Body'].y])
                tail_xy = np.array([i['Tail'].x,i['Tail'].y])
            except:
                head_xy = np.array([i['head'].x,i['head'].y])
                body_xy = np.array([i['body'].x,i['body'].y])
                tail_xy = np.array([i['tail'].x,i['tail'].y])
            hb_dist = np.linalg.norm(head_xy - body_xy)
            bt_dist = np.linalg.norm(body_xy - tail_xy)
            length = hb_dist + bt_dist
            vid_dict[vid][fish_index,count] = length
            vid_index = vid_indices[vid]
            all_tracks[vid_index,fish_index,count] = length
            vid_counts[vid][fish_index] += 1

    mean_size = np.nanmean(all_tracks,axis=(1,2))
    size_variance_within = np.nanmean(np.nanstd(all_tracks,axis=2),axis=1)
    size_variance_among = np.nanmean(np.nanstd(all_tracks,axis=1),axis=1)
    rpt0 = size_variance_among / (size_variance_within + size_variance_among)

    means_quartiled = np.empty([len(all_tracks),4])
    rpt_list = []

    for f in range(len(all_tracks)):
        df_list = []
        sorted_sizes = np.sort(np.nanmean(all_tracks[f],axis=1))
        means_quartiled[f] = sorted_sizes
        for i in range(4):
            for j in range(5):
                df_list.append([i,j,all_tracks[f,i,j]])
        df = pd.DataFrame(df_list,columns=['fish','frame','size'])
        try:
            model = ols('size ~ C(fish)',data=df).fit()
            aov_table = sm.stats.anova_lm(model,typ=2)
            resid = aov_table['sum_sq'][1]
            var_a = aov_table['sum_sq'][0]
            rpt = var_a / (var_a + resid)
            rpt_list.append(rpt)
        except:
            rpt = np.nan
            var_a = np.nan
            resid = np.nan
            pass
            #import pdb;pdb.set_trace()

        vid_list = [pi,f,sorted_sizes[0],sorted_sizes[1],sorted_sizes[2],sorted_sizes[3],np.nanmean(sorted_sizes),np.nanmax(sorted_sizes),var_a,resid,rpt]
        vid_string = [str(v) for v in vid_list]
        mega_df_list.append(','.join(vid_string))

    if True:
        mean_size = mean_size * (23/800)
        means_quartiled = means_quartiled * (23/800)
    for i in range(4):
        #ax1.plot(np.nanmean(all_tracks[:,i],axis=1))
        ax1.plot(means_quartiled[:,i],color=cmap(a_/3),alpha=0.5)
    ax1.plot(mean_size,color=cmap(a_/3),linewidth=2,label=pi)
    ax1.set_xlabel('day')
    ax1.set_ylabel('Size (cm)')
#ax2.plot(size_variance_within)
#ax2.plot(size_variance_among)
#ax2.plot(rpt0)
    ax2.plot(rpt_list,color=cmap(a_/3))

ax2.set_xlabel('day')
ax2.set_ylabel('rpt')
ax1.legend()

fig.set_size_inches([4,8])
fig.tight_layout()
fig.savefig('pi_a11_track.png',dpi=300)
plt.show()
