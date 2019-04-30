Each identified filament has its own txt file that is organized as (x y v) coordinates of the pixels along the corresponding filament. 

Two points I should mention: 
1. Because the c18o cube was too large, I divided it into two parts in order to be able to identify the filaments – hence the folders containing northern and southern filaments. The northern and southern c18o cubes are also attached. These will help you convert the identified filament pixel coordinates (using, for example, astropy wcs) to the pixel coordinates of whichever map you are using. 

2. Some txt files may be empty or have something like 2 or 5 pixels. You can eliminate these by setting up a threshold when you read the file, saying if it has less than 15 pixels (roughly 3 beam sizes), ignore the filament.