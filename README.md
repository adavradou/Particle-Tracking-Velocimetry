# About
Particle Tracking Velocimetry code applied in noisy environments, e.g. glass channel. 
It performs "noise mapping" of the recurrent noise noticed through the experiments and creates a background image. 
The background image is then deleted from all the images, resulting in data with less noise.

The functions bpass.m, pkfnd.m, cntrd.m and track.m were found here http://site.physics.georgetown.edu/matlab/code.html.

The PTV_create_backround_image_same_th.m creates the backgound image. It must be run first to create the background image and then run the main code.

The PTV_remove_noise.m is the main PTV code that utilizes all the above functions and the background image. 
