import sleap
import numpy as np
from matplotlib import pyplot as plt
#from matplotlib import colormaps
import scipy.stats as stats
import pandas as pd
import os
import seaborn as sns

from datetime import datetime

from statsmodels.formula.api import ols
import statsmodels.api as sm

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

#anns = sleap.load_file('/home/ammon/Downloads/pi11.Labeled/groupTracks.pi11.slp')
ann_dir = './labeledTracks/'

n_anns = len(os.listdir(ann_dir))

mega_df_list = []
fig,(ax1,ax2) = plt.subplots(2)
cmap = plt.get_cmap('viridis')

day_dict = {
    'pi11':'18/04/05',
    'pi12':'18/03/04',
    'pi13':'18/03/20',
    'pi14':'18/03/01',
    'pi31':'18/03/30',
    'pi32':'18/02/16',
    'pi33':'18/03/27',
    'pi34':'18/03/27',
    'pi41':'18/03/21',
    'pi42':'18/03/21',
    'pi51':'18/03/20',
    'pi52':'18/03/21',
    'pi53':'18/03/20',
    'pi54':'18/03/21'}

file_to_day = {}

beh_df = pd.read_csv("JolleTracksAll_8.csv")
variance_list = []
for a_ in range(n_anns):
    a = sorted(os.listdir(ann_dir))[a_]
    #print(anns)
    print(a_,a)
    ann_file = ann_dir + a
    anns = sleap.load_file(ann_dir + a)
    pi = a.split('.')[1]

    n_vids = len(anns.videos)
    video_anns = [[] for v in range(n_vids)]

#vid_list = []
    vid_dict = {}
    vid_counts = {}
    vid_indices = {}
    vid_days = [] 
## Each video needs an array, with each fish having a little track
    all_tracks = np.full([n_vids,4,5],np.nan)
    count = 0

    #import pdb;pdb.set_trace()
    for v in anns.videos:
        #vid_list.append(v.filename)
        vid_indices[v.filename] = count
        count += 1
        basename = v.filename.split('/')[-1]
        pi,date_ = basename.split('_')
        date = date_.replace('.clip.mp4','')
        date = date.replace('.mp4','')
        start_day = datetime.strptime(day_dict[pi],'%y/%m/%d')
        vid_day = datetime.strptime(date,'%y%m%d')
        diff = (vid_day - start_day).days
        file_to_day[v.filename] = diff
        vid_dict[v.filename] = np.full([4,5],np.nan) 
        vid_counts[v.filename] = [0 for i in range(4)]
        vid_days.append(diff)

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
                continue

            node_names = sorted([n.name for n in i.nodes]) ## I don't want to talk about it.
            head_xy = np.array([i[node_names[1]].x,i[node_names[1]].y])
            body_xy = np.array([i[node_names[0]].x,i[node_names[0]].y])
            tail_xy = np.array([i[node_names[2]].x,i[node_names[2]].y])
            #try:
            #    head_xy = np.array([i['Head'].x,i['Head'].y])
            #    body_xy = np.array([i['Body'].x,i['Body'].y])
            #    tail_xy = np.array([i['Tail'].x,i['Tail'].y])
            #except:
            #    head_xy = np.array([i['head'].x,i['head'].y])
            #    body_xy = np.array([i['body'].x,i['body'].y])
            #    tail_xy = np.array([i['tail'].x,i['tail'].y])
            hb_dist = np.linalg.norm(head_xy - body_xy)
            bt_dist = np.linalg.norm(body_xy - tail_xy)
            length = hb_dist + bt_dist
            
            vid_dict[vid][fish_index,count] = length
            #import pdb;pdb.set_trace()
            vid_index = vid_indices[vid]
            all_tracks[vid_index,fish_index,count] = length
            vid_counts[vid][fish_index] += 1

    for vid in vid_dict.keys():
        expDay = file_to_day[vid]
        fish_means = np.nanmean(vid_dict[vid],axis=0)
        fish_err = np.nanmean(np.nanstd(vid_dict[vid],axis=0))
        fish_cv = fish_err / np.mean(fish_means)
        if sum(~np.isnan(fish_means)) == 0:
            continue
        day_std = np.nanstd(fish_means)
        mean_size = np.mean(fish_means)
        min_size = np.min(fish_means)
        max_size = np.max(fish_means)
        norm_std = day_std / mean_size
        pi = vid.split('/')[-1].split('_')[0]
        vid_df = beh_df[beh_df.Pi == pi]
        day_df = vid_df[vid_df.ExpDay == expDay]
        day_vel = day_df.vel_mean.values[0]
        day_dist = day_df.dist_mean.values[0]
        day_velC = day_df.velC_mean.values[0]
        day_pDist = day_df.pDist_mean.values[0]
        day_pDistC = day_df.pDistC_mean.values[0]
        day_angC = day_df.angleC_mean.values[0]
        variance_list.append([expDay,pi,mean_size,min_size,max_size,day_std,norm_std,fish_err,fish_cv,day_dist,day_vel,day_velC,day_pDist,day_pDistC,day_angC])

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
            #resid = aov_table['sum_sq'][1]
            #var_a = aov_table['sum_sq'][0]

            resid = aov_table['sum_sq'].iloc[1]
            var_a = aov_table['sum_sq'].iloc[0]
            rpt = var_a / (var_a + resid)
        except:
            rpt = np.nan
            var_a = np.nan
            resid = np.nan
            pass
            #import pdb;pdb.set_trace()

        rpt_list.append(rpt)
        day = vid_days[f]
        vid_list = [pi,f,day,sorted_sizes[0],sorted_sizes[1],sorted_sizes[2],sorted_sizes[3],np.nanmean(sorted_sizes),np.nanmax(sorted_sizes),var_a,resid,rpt]
        vid_string = [str(v) for v in vid_list]
        #mega_df_list.append(','.join(vid_string))
        mega_df_list.append(vid_list)

    vid_days = np.array(vid_days)
    rpt_list = np.array(rpt_list) 
    if True:
        mean_size = mean_size * (23/800)
        means_quartiled = means_quartiled * (23/800)
    if pi == 'pi11':
        a_scale = 1
    else:
        a_scale = 0.4
    if pi == 'pi11' and False:
        for i in range(4):
            #ax1.plot(np.nanmean(all_tracks[:,i],axis=1))
            ax1.plot(vid_days[~np.isnan(mean_size)],means_quartiled[~np.isnan(mean_size)][:,i],color='red',alpha=0.5*a_scale)
    if a_ == 0:
        cor = 'red'
    else:
        cor = 'black'
    cor = cmap(a_/9)
    if pi == 'pi11' or True:
        xs = vid_days[~np.isnan(mean_size)]
        ys = mean_size[~np.isnan(mean_size)]
        ys = ys[np.argsort(xs)]
        xs = xs[np.argsort(xs)]
        ax1.plot(xs,ys,color=cor,linewidth=2,label=pi,alpha=1*a_scale)
    else:
        pass
        
    ax2.plot(vid_days[~np.isnan(mean_size)],rpt_list[~np.isnan(mean_size)],color=cor)

variance_df = pd.DataFrame(variance_list,columns=['ExpDay','Pi','mean_size','min_size','max_size','std','norm_std','indv_err','indv_cv','dist_mean','vel_mean','velC_mean','pDist_mean','pDistC_mean','angC_mean'])

variance_df.to_csv('variance_df.csv',index=False)

#ax2.plot(size_variance_within)
#ax2.plot(size_variance_among)
#ax2.plot(rpt0)
ax1.set_xlabel('day')
ax1.set_ylabel('Size (cm)')
ax2.set_xlabel('day')
ax2.set_ylabel('rpt')
ax1.legend()

fig.set_size_inches([4,8])
fig.tight_layout()
fig.savefig('pi_a11_track.png',dpi=300)
#plt.show()
columns = ['Pi','index','ExpDay','size1','size2','size3','size4','mean_size','max_size','var_a','resid','rpt']

## Build correlation plots
size_df = pd.DataFrame(mega_df_list,columns=columns)

rs = []
ps = []
days = []

size_df.to_csv('size_df2.csv',index=False)

fig,axes = plt.subplots(3,1,sharex=True,sharey=True)
day_list = [0,24,54]
for d in range(len(pd.unique(size_df.ExpDay))):
#for d_ in range(len(day_list)):
    #d = pd.unique(size_df.day)[d_]
    #d = day_list[d_]
    day_df = size_df[size_df.ExpDay == d]
    day_beh = beh_df[beh_df.ExpDay == d]

    if len(day_df) <= 3:
        continue
    xs = []
    ys = [] 
    for _,row in day_df.iterrows():
        pi = row.Pi
        xs.append(row.mean_size)
        ys.append(day_beh[day_beh.Pi == row.Pi].dist_mean.values[0])

    if np.sum(np.isnan(xs)) > 0:
        continue
    r,p = stats.pearsonr(xs,ys)
    xs = np.array(xs) * (23/800)
    ys = np.array(ys) * (23/800)
    if d in day_list:
        if d == 0:
            d_ = 0
        elif d == 24:
            d_ = 1
        else:
            d_ = 2
        ax = axes[d_]
        if False:
            ax.scatter(xs,ys)
            ax.plot(np.unique(xs), np.poly1d(np.polyfit(xs, ys, 1))(np.unique(xs)))
        else:
            if p < 0.05:
                line_cor = 'red'
            else:
                line_cor = 'black'
            sns.regplot(x=xs,y=ys,ax=ax,
            scatter_kws={"color":"black"},
            line_kws={"color":line_cor})
        ax.set_title(' '.join([str(np.round(r,3)),str(np.round(p,3))]))
    rs.append(r)
    ps.append(p)
    days.append(d)
    fig.tight_layout()
    fig.set_size_inches([5.5,5.5])

print(ps)
if True:
    fig2,ax = plt.subplots()
    ax.plot(days,rs,color='black')
    ax.axhline(0,color='black',linestyle=':')
    ax.set_ylim([-1,1])
    ax.set_ylabel('Pearson coefficient')
    ax.set_xlabel('Day')
    fig2.set_size_inches([4,4])
    fig2.tight_layout()
#fig.savefig('pearson_time.png',dpi=300)
if True:
    plt.show()

