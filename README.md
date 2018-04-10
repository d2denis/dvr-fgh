# discrete  variable representation Fourier grid Hamiltonian (DVR-FGH)
# 离散变量表示傅里叶格点哈密顿量

对于多模式(nmode>2)的情况下，location数组过大，厚道的python版无法进行。
此程序为C++版本，丰富了一些功能，能一次性的产生DEOM所需输入文件。
另外改进的主要是每个模式挑选各自的最低能级：
1) 多模式格点的编码(encode)和解码(decode)分别有一套构建哈密顿量的程序；
2) 沿用批处理法，构建稀疏矩阵哈密顿；
3) 厚道师兄的算法解出来的本征值从小到大排列，c++版的arpack则只能直接调节tol来控制收敛。
目前解码法构建哈密顿量目前最优。