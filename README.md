# Earthquake-to-Fault
This is a self-adaptive fault detection tool from high-precision earthquake catalog based on Hough transform

The E2F (Earthquake-to-Fault) is released and maintained at https://github.com/User0324-student/Earthquake-to-Fault

How to start? 1): Please install hough-3d-lines (https://github.com/cdalitz/hough-3d-lines) at first. 2): Run the replace.py.

1.set the correct roads of "hough-3d-lines" and Earthquake-to-fault (E2F) for replacing the code of "hough3dlines.cpp", "pointcloud.cpp", "vector3d.cpp" and "vector3d.h"

2.Now, you can run the main program (E2F_1_0.m)

NOTE: please set the correct earthquake catalog file like:

# Event numbers | Latitude | Longitude | Depth | Year | Mon | Day | Hour | Minute | Second | Mag
  1	| 36.8081 | -97.6997 | 7.5201	| 2013 | 6 | 14 | 13 | 29	| 22.4900 | 0.960000000000000
  2	| 36.8085	| -97.7019 | 7.5300	| 2013 | 6 | 15	| 14 | 52	| 1.95000 | 0.910000000000000
  3	| 36.8088	| -97.7017 | 7.5700	| 2013 | 6 | 16	| 11 | 12	| 18.3800 | 0.610000000000000 
  ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ...

