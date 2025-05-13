

for i in $(cat clips_to_check.txt); do
    echo "$i"
    vid_name=${i##*/}
    vid_name=${vid_name%%\"}
    echo $vid_name
    out_path="./check_clips/clips."$vid_name
    rclone copy aperkes:"$i" ./working_dir -P
### overlay tracks 
    loop_string="144 288 432 476 720"
    for s in $loop_string; do
        tmpfile="tmp$s.mp4"
        #echo /home/ammon/Documents/Scripts/GroupMollyTracking/$tmpfile >> myTimes.txt
        ffmpeg -i ./working_dir/overlay.$vid_name -ss $s -t 3 $tmpfile -y
        done
    ffmpeg -f concat -safe 0 -i myTimes.txt -c copy $out_path -y
    done < ./clips_to_check.txt
