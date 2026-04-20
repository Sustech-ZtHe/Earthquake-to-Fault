import os
import shutil

# 当前目录（就是 replace.py 所在目录）
base_path = os.getcwd()

# 源（替换文件）
replace_dir = os.path.join(base_path, "replace_code")

# 目标（被替换文件）
target_dir = os.path.join(base_path, "hough-3d-lines-master")

# 文件映射关系
file_map = {
    "hough3dlines_replaced.cpp": "hough3dlines.cpp",
    "pointcloud_replaced.cpp": "pointcloud.cpp",
    "vector3d_replaced.cpp": "vector3d.cpp",
    "vector3d_replaced.h": "vector3d.h",
}

# 执行替换
for src_name, dst_name in file_map.items():
    src = os.path.join(replace_dir, src_name)
    dst = os.path.join(target_dir, dst_name)

    print(f"复制: {src} -> {dst}")

    if not os.path.exists(src):
        print(f"源文件不存在: {src}")
        continue

    shutil.copyfile(src, dst)

print("替换完成")
