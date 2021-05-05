# mpiExp
## mpi程序编译
1. mpic++ PSRSSort.cpp -o PRRSSortFloat
2. mpic++ [cpp_file] -o [target name]
## 运行
1. 编辑host文件，这里是mpi_config,指定主机名和核数
如：
![image](https://user-images.githubusercontent.com/33123364/117093025-5a7e7700-ad92-11eb-8141-dd31fb4f8e97.png)
2.PRRS排序：
``mpirun -np [process num] -f [host_file] ./PRRSSortFloat [pow] [file_name]``
3.奇偶排序
``mpirun -np 4 -f [host_file] ./OddEvenSort [file_name]``
