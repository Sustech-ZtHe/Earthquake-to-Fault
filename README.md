# Earthquake-to-Fault
This is a self-adaptive fault detection tool from high-precision earthquake catalog based on Hough transform
The E2F (Earthquake-to-Fault) is released and maintained at https://github.com/Sustech-ZtHe/Earthquake-to-Fault
E2F is based on the matlab framework, the main program is E2F_1_0.m, the toolbox of matlab (R2023b) required:
1.Mapping Toolbox
2.Machine Learning Toolbox 
3.Optimization Toolbox
4.Text Analytics Toolbox

How to start? 
*Install the hough-3d-lines (https://github.com/cdalitz/hough-3d-lines)
*Install the Eigen (http://eigen.tuxfamily.org/) 
*Let the three zips (hough-3d-lines-master, eigen-3.4.0, Earthquake-to-fault-main) be placed in a directory such as 'H3D'.
Then
1): Setup the Eigen (C++ library) at first
  1.1 mkdir build
      cd build
      cmake ..
  1.2 Enter the hough-3d-lines-master, find Makefile and change the road of 'LIBEIGEN='usr/include/eigen3' like 'LIBEIGEN='your_H3D_place/eigen-3.4.0'
2): Enter the Earthquake-to-fault-main, run the 'replayce.py'. 
3): Setup hough-3d-lines 
  make
4) Enter the ../H3D/hough-3d-lines-master and Run main code E2F_1_0.m 

NOTE: 
1. We recommended that 'eigen-3.4.0','hough-3d-lines-master' and 'E2F_package' be placed in the same directory.
2. We have placed the test catalog (ToC2ME) on the path of '../E2F_package/SaveData/ToC2ME.txt', if you want use other catalog, the earthquake catalog arrangement format like: 
# Event numbers | Year | Mon | Day | Hour | Minute | Second |  Latitude | Longitude |  Depth | Mag |
3. Chat with us on wechat (H-explorer) or send email (hezhengtao2001@163.com)

Some promble maybe happen:
*`GLIBCXX_3.4.29‘ not found 
  The matlab (R2023b) need `GLIBCXX_3.4.29‘, check the libstdc++.so.6 ('../Matlab/sys/os/glnxa64') and cp `GLIBCXX_3.4.29‘ in this directory.

![image](https://github.com/user-attachments/assets/28cee59a-2ecd-4ad7-8986-cfe2d76a7d4c)
