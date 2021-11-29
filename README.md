# raytracer

## Creating video from frames
All the 1000 frames will be saved on /frames/frame(number of frame).jpg,
you can then run
`ffmpeg -r 60 frame%d.jpg -s 1920x1080 video.avi`
and it will encode the frames into a video.