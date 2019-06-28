# [Generalized Video Deblurring for Dynamic Scenes (CVPR15)](http://cv.snu.ac.kr/publication/conf/2015/VD_CVPR2015.pdf "Dynamic Video Deblurring using a Locally Adaptive Linear Blur Model ") #

Contact Tae Hyun Kim (lliger9(at)snu.ac.kr).
 
[[Paper](http://cv.snu.ac.kr/publication/conf/2015/VD_CVPR2015.pdf "Paper")]
[[Project Page](http://cv.snu.ac.kr/research/~VD/ "Project Page")]
[[Code](http://cv.snu.ac.kr/research/~VD/VD-release_v1.0.zip "Code")]



***Any work that uses our code should cite our paper (see below).**


## *Instruction to run sample code ##

Pre-compilationa was done under windows7 64bit with matlab R2015a.
Our non-optimized code requires large memory. 
(We use 32GB RAM to deblur 50 frames with 1280X720 resolution.)


1. run video_deblurring.m to run code for "car" example.



## *Instruction to deblur your own examples ##

1. Name images from 0001.png. (Any extension that can be readable with MATLAB is ok.)
2. run 'TV\_L1\_OF\test\_motion.m' to obtain initial optical flows. 
(Chambolle, Antonin, and Thomas Pock. "A first-order primal-dual algorithm for convex problems with applications to imaging." Journal of Mathematical Imaging and Vision 40.1 (2011): 120-145.) You can use any other method if you want.

3. set proper params in video_deblurring.m

	lambda, tau; params in the paper
	
	n\_imlist; number of frames to deblur

	pyramid_levels: define pyramid levels used in coarse to fine approach. We usually use 17 levels for 1280X720 resolution and use less(~14) for smaller images.

4. run video_deblurring.m to deblur your own examples.


## *References ##

Tae Hyun Kim and Kyoung Mu Lee, "Generalized Video Deblurring for Dynamic Scenes," 
Computer Vision and Pattern Recognition (CVPR), 2015 IEEE Conference on. IEEE, 2015.

Tae Hyun Kim and Kyoung Mu Lee. "Segmentation-free dynamic scene deblurring." 
Computer Vision and Pattern Recognition (CVPR), 2014 IEEE Conference on. IEEE, 2014.