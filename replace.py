import shutil
import os
current_path=os.getcwd()
dir_path=os.path.dirname(current_path)
h3d_road=f"{dir_path}/hough-3d-lines-master"
e2f_road=f"{dir_path}/Earthquake-to-Fault-main"
fault_road=f"{dir_path}/OutputFile"
Catalog_example_road=f"{e2f_road}/Catalog_example"
e2fpackage_road=f"{dir_path}/Earthquake-to-Fault-main/E2F_package"
e2fgui_road=f"{dir_path}/Earthquake-to-Fault-main/E2F_v1.m"
e2fguifig_road=f"{dir_path}/Earthquake-to-Fault-main/E2F_v1.fig"
os.makedirs(fault_road, exist_ok=True)

# origin file road
#origin_h3d_main= "/your_hough-3d-line_road/hough3dlines.cpp"
#origin_pointcloud= "/your_hough-3d-line_road/pointcloud.cpp"
#origin_vector3d= "/your_hough-3d-line_road/vector3d.cpp"
#origin_vector3d_h= "/your_hough-3d-line_road/vector3d.h"
origin_h3d_main= f"{h3d_road}/hough3dlines.cpp"
origin_pointcloud= f"{h3d_road}/pointcloud.cpp"
origin_vector3d= f"{h3d_road}/vector3d.cpp"
origin_vector3d_h= f"{h3d_road}/vector3d.h"
# replace file road
replace_h3d_main = f"{e2f_road}/replace_code/hough3dlines_replaced.cpp"
replace_pointcloud= f"{e2f_road}/replace_code/pointcloud_replaced.cpp"
replace_vector3d= f"{e2f_road}/replace_code/vector3d_replaced.cpp"
replace_vector3d_h= f"{e2f_road}/replace_code/vector3d_replaced.h"

# replace file
shutil.copyfile(replace_h3d_main, origin_h3d_main)
shutil.copyfile(replace_pointcloud, origin_pointcloud)
shutil.copyfile(replace_vector3d, origin_vector3d)
shutil.copyfile(replace_vector3d_h, origin_vector3d_h)

#move GUI and .fig
shutil.move(e2fpackage_road,dir_path)
shutil.move(e2fgui_road,dir_path)
shutil.move(e2fguifig_road,dir_path)
shutil.move(Catalog_example_road,dir_path)
