import shutil

hough3dlines_road="..."
E2F_road="..."

# 源文件路径
origin_h3d_main= f"{hough3dlines_road}/hough3dlines.cpp"
origin_pointcloud= f"{hough3dlines_road}/pointcloud.cpp"
origin_vector3d= f"{hough3dlines_road}/vector3d.cpp"
origin_vector3d_h= f"{hough3dlines_road}/vector3d.h"
# 替换文件路径
replace_h3d_main = f"{E2F_road}/replace_code/hough3dlines_replaced.cpp"
replace_pointcloud= f"{E2F_road}/replace_code/pointcloud_replaced.cpp"
replace_vector3d= f"{E2F_road}/replace_code/vector3d_replaced.cpp"
replace_vector3d_h= f"{E2F_road}/replace_code/vector3d_replaced.h"

# 将 replace 的内容复制到 origin
shutil.copyfile(replace_h3d_main, origin_h3d_main)
shutil.copyfile(replace_pointcloud, origin_pointcloud)
shutil.copyfile(replace_vector3d, origin_vector3d)
shutil.copyfile(replace_vector3d_h, origin_vector3d_h)
