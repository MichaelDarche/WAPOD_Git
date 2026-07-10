# https://pythonexamples.org/python-opencv-cv2-create-video-from-images/

import cv2
import os

image_folder = '.'
video_name = 'Bloch-Band.avi'

#images = [img for img in os.listdir(image_folder) if img.endswith(".jpg")]
images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
frame = cv2.imread(os.path.join(image_folder, images[0]))
height, width, layers = frame.shape
print(height, width)

video = cv2.VideoWriter(video_name, 0, 25, (width,height))

for image in images:
    print(image)
    video.write(cv2.imread(os.path.join(image_folder, image)))

cv2.destroyAllWindows()
video.release()
print('termine !')
