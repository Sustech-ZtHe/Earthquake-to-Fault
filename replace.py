import shutil

# 源文件路径
origin_h3d_main= "/home/me/Documents/Ridgecrest/hough3dlines.cpp"
origin_pointcloud= "/home/me/Documents/Ridgecrest/pointcloud.cpp"
origin_vector3d= "/home/me/Documents/Ridgecrest/vector3d.cpp"
origin_vector3d_h= "/home/me/Documents/Ridgecrest/vector3d.h"
# 替换文件路径
replace_h3d_main = "/home/me/Documents/replace_code/hough3dlines_replaced.cpp"
replace_pointcloud= "/home/me/Documents/replace_code/pointcloud_replaced.cpp"
replace_vector3d= "/home/me/Documents/replace_code/vector3d_replaced.cpp"
replace_vector3d_h= "/home/me/Documents/replace_code/vector3d_replaced.h"

# 将 replace 的内容复制到 origin
shutil.copyfile(replace_h3d_main, origin_h3d_main)
shutil.copyfile(replace_pointcloud, origin_pointcloud)
shutil.copyfile(replace_vector3d, origin_vector3d)
shutil.copyfile(replace_vector3d_h, origin_vector3d_h)
