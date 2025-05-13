

for i in $(cat clips_to_check.txt); do
    echo "$i"
    vid_name=${i##*/}
    vid_name=${vid_name%%\"}
    echo $vid_name
    rclone copy aperkes:"$i" ./working_dir -P
### overlay tracks 
    i_file=$(basename $i)
    python overlayGroupTrack.py ./working_dir/$i_file
    done
