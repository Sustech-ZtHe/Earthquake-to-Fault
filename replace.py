import shutil
import os
current_path=os.getcwd()
dir_path=os.path.dirname(current_path)
h3d_road=f"{dir_path}/hough-3d-lines-master"
fault_road=f"{h3d_road}/FaultSegment"
e2f_road=f"{dir_path}/Earthquake-to-Fault-main"
e2fpackage_road=f"{dir_path}/Earthquake-to-Fault-main/E2F_package"
os.makedirs(fault_road, exist_ok=True)
# 源文件路径
#origin_h3d_main= "/your_hough-3d-line_road/hough3dlines.cpp"
#origin_pointcloud= "/your_hough-3d-line_road/pointcloud.cpp"
#origin_vector3d= "/your_hough-3d-line_road/vector3d.cpp"
#origin_vector3d_h= "/your_hough-3d-line_road/vector3d.h"
origin_h3d_main= f"{h3d_road}/hough3dlines.cpp"
origin_pointcloud= f"{h3d_road}/pointcloud.cpp"
origin_vector3d= f"{h3d_road}/vector3d.cpp"
origin_vector3d_h= f"{h3d_road}/vector3d.h"
# 替换文件路径
replace_h3d_main = f"{e2f_road}/replace_code/hough3dlines_replaced.cpp"
replace_pointcloud= f"{e2f_road}/replace_code/pointcloud_replaced.cpp"
replace_vector3d= f"{e2f_road}/replace_code/vector3d_replaced.cpp"
replace_vector3d_h= f"{e2f_road}/replace_code/vector3d_replaced.h"


# 将 replace 的内容复制到 origin
shutil.copyfile(replace_h3d_main, origin_h3d_main)
shutil.copyfile(replace_pointcloud, origin_pointcloud)
shutil.copyfile(replace_vector3d, origin_vector3d)
shutil.copyfile(replace_vector3d_h, origin_vector3d_h)

shutil.move(e2fpackage_road,h3d_road)
