# Earthquake-to-Fault
This is a self-adaptive fault detection tool from high-precision earthquake catalog based on Hough transform

The E2F (Earthquake-to-Fault) is released and maintained at https://github.com/User0324-student/Earthquake-to-Fault

How to start? 
1): Please install hough-3d-lines (https://github.com/cdalitz/hough-3d-lines) at first. 
2): Matlab requires the Mapping toolbox and Statistics and Machine Learning Toolbox. 
3): Run the replace.py.

NOTE: 
1. We recommended that 'hough-3d-lines-master' and 'E2F_package' be placed in the same directory.
2. We have placed the ToC2ME on the path of '../E2F_package/SaveData/ToC2ME.txt', the earthquake catalog arrangement format like: 

# Event numbers | Year | Mon | Day | Hour | Minute | Second | Latitude | Longitude | Depth | Mag

