
import sys
from defineRoles import plot_video, get_tracks
from processTracks import deep_clean_track

csv_file = sys.argv[1]
vid_file = sys.argv[2]
track_array,track_polar, [n_frames,n_fish,fishIDs] = get_tracks(csv_file)

clean_array = deep_clean_track(track_array,drop_close = True)
plot_video(clean_array,vid_file,viz=True,trail=10)
