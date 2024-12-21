E2F is based on the MATLAB framework, the main program is E2F_V1.m  
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
1): Setup the Eigen (C++ library) at first  
  1.1 mkdir build  
      cd build  
      cmake ..  
  1.2 Enter the hough-3d-lines-master, find Makefile and change the road of 'LIBEIGEN='usr/include/eigen3' like 'LIBEIGEN='your_H3D_place/eigen-3.4.0'
2): Enter the Earthquake-to-fault-main, run the 'replayce.py'. 
3): Setup hough-3d-lines 
  make
4) Add the 'E2F_package' into the MATLAB road
5) Run the main project E2F_v1.m 

NOTE: 
1. We recommended that 'eigen-3.4.0','hough-3d-lines-master' and 'E2F_package' be placed in the same directory.
2. We have placed the test catalog (ToC2ME) on the path of '..Catalog_example/ToC2ME.txt', if you want use other catalog, the earthquake catalog arrangement format like: 
# Event numbers | Year | Mon | Day | Hour | Minute | Second | Latitude | Longitude | Depth | Mag

Some promble maybe happen:
1.If you some problem like "GLIBCXX_3.4.29 not found" when open the E2F_v1 and run 'Catalog File'.
  This is because your MATLAB lacks a patch of 'GLIBCXX_3.4.29' ro tun hough-3d-lines.
  You could run 'strings /usr/lib/x86_64-linux-gnu/libstdc++.so.6 | grep GLIBCXX' in terminal if it contains 'GLIBCXX_3.4.29'. 
  Try running 'sudo cp /usr/lib/x86_64-linux-gnu/libstdc++.so.6' your_matlab_path/libstdc++.so.6 '.
  

![Toc2_v1](https://github.com/user-attachments/assets/f86fbd3b-80e5-418f-acb2-c3cf0a2c2aed)


Chat with us on wechat (H-explorer) or send email (hezhengtao2001@163.com)
