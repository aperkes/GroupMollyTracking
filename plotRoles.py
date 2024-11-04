
import sys
from defineRoles import plot_video, get_tracks


csv_file = sys.argv[1]
vid_file = sys.argv[2]
track_array,track_polar, [n_frames,n_fish,fishIDs] = get_tracks(csv_file)

plot_video(track_array,vid_file,viz=True)
