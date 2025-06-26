E2F is based on the MATLAB framework, and it is tested on Linux Ubuntu 22.04.5 (suggested) and Widows  
The main script is E2F_V1.m  
# Earthquake-to-Fault  
This is a self-adaptive fault detection tool from high-precision earthquake catalog based on Hough transform  
The E2F (Earthquake-to-Fault) is released and maintained at https://github.com/Sustech-ZtHe/Earthquake-to-Fault  
The toolbox of MATLAB (R2023b) required:  
1.Mapping Toolbox  
2.Machine Learning Toolbox  
3.Optimization Toolbox  

How to start?  
*Install the hough-3d-lines (https://github.com/cdalitz/hough-3d-lines)  
*Install the Eigen (http://eigen.tuxfamily.org/)  
*Let the three zips (hough-3d-lines-master, eigen-3.4.0, Earthquake-to-fault-main) be placed in a directory such as 'E2F'.  
Then  
1): Setup the Eigen (C++ library) at first (For Windows systems, proceed directly to step 1.2.)  
&nbsp;&nbsp;1.1 mkdir build  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cd build  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cmake ..  
&nbsp;&nbsp;1.2 Enter the hough-3d-lines-master, find Makefile and change the road of 'LIBEIGEN='usr/include/eigen3' like 'LIBEIGEN='your_H3D_place/eigen-3.4.0'  
2): Enter the Earthquake-to-fault-main, run the 'replayce.py'.  
3): Setup hough-3d-lines  
&nbsp;&nbsp;make  
4) Run the main script is E2F_V1.m, and add the 'E2F_package' folder to the MATLAB path 

NOTE: 
1. We recommended that 'eigen-3.4.0','hough-3d-lines-master' and 'E2F_package' be placed in the same directory.
2. We have placed the test catalog (ToC2ME) on the path of '..Catalog_example/ToC2ME.txt', if you want use other catalog, the earthquake catalog arrangement format like:  
**ID | Year | Mon | Day | Hour | Minute | Second | Latitude | Longitude | Depth | Magnitude**

**Some promble maybe happen:**  
1.If you some problem like "GLIBCXX_3.4.29 not found" when open the E2F_v1 and run 'Catalog File'.  
  This issue occurs because your MATLAB environment lacks the required version of GLIBCXX_3.4.29, which is necessary to run hough-3d-lines.  
  You can check if your system has this version by running the following command in the terminal:  
&nbsp;&nbsp;*strings /usr/lib/x86_64-linux-gnu/libstdc++.so.6 | grep GLIBCXX_3.4.29*  
  If the required version is available, try copying the library to the MATLAB directory by running:  
&nbsp;&nbsp;*sudo cp /usr/lib/x86_64-linux-gnu/libstdc++.so.6 your_matlab_path/libstdc++.so.6*  

![Toc2_v1](https://github.com/user-attachments/assets/f86fbd3b-80e5-418f-acb2-c3cf0a2c2aed)

# Output File  
**Fault_Segment_Clustered.txt**  
#1 Event Easting (km)  
#2 Event Northing (km)  
#3 Depth (km)
#4 Magnitude  
#5-10 Year Month Day Hour Minute Second  
#11-13 R G B  
#14 Fault ID  
-117.239 54.3457 3.41 0.36 2016 11 5 19 5 54.8 0.814724 0.153814 0.436224 1  

**Fault_Segment_Modeling.txt**  
#1 Fault ID  
#2-4 Centroid X (km) Centroid Y (km) Centroid Z (km)  
#5-7 Major Axis Upper Endpoint X (km) Major Axis Upper Endpoint Y (km) Major Axis Upper Endpoint Z (km)  
#8-10 Major Axis Lower Endpoint X (km) Major Axis Lower Endpoint Y (km) Major Axis Lower Endpoint Z (km)  
#11-13 Intermediate Axis Upper Endpoint X (km) Intermediate Axis Upper Endpoint Y (km) Intermediate Axis Upper Endpoint Z (km)  
#14-16 Intermediate Axis Lower Endpoint X (km) Intermediate Axis Lower Endpoint Y (km) Intermediate Axis Lower Endpoint Z (km)  
#17-19 Short Axis Upper Endpoint X (km) Short Axis Upper Endpoint Y (km) Short Axis Upper Endpoint Z (km)  
#20-22 Short Axis Lower Endpoint X (km) Short Axis Lower Endpoint Y (km) Short Axis Lower Endpoint Z (km)  
#23 Strike (deg)  
#24 Dip Angle (deg)  
#25 Optimal C Value  
#26 Number of Events constituting the fault  
1 0.909659 1.57908 3.385 1.03479 1.81128 3.45246 0.753679 1.30613 3.31802 0.886865 1.54966 3.43464 0.901604 1.56775 3.33584 0.908689 1.55049 3.38589 0.879781 1.56692 3.38459 29.0955 87.7548 3.5 4489  


Feel free to contact us on WeChat (H-explorer) or via email (hezhengtao2001@163.com).  
The Python version will be launched soon.
