# Pattern Recognition-2023-HMRMSC  

This repo contains the MATLAB code of our Pattern Recognition'2023 article: [High-order manifold regularized multi-view subspace clustering with robust affinity matrices and weighted TNN](https://www.sciencedirect.com/science/article/abs/pii/S0031320322005477).

## Requirements

- MATLAB R2018b

## Run

```
Demo_main.m
```

## Note

If you have any question, please contact this e-mail: <bingly@foxmail.com>.

```latex
@article{CAI2023109067,
title = {High-order manifold regularized multi-view subspace clustering with robust affinity matrices and weighted TNN},
journal = {Pattern Recognition},
volume = {134},
pages = {109067},
year = {2023},
issn = {0031-3203},
doi = {https://doi.org/10.1016/j.patcog.2022.109067},
url = {https://www.sciencedirect.com/science/article/pii/S0031320322005477},
author = {Bing Cai and Gui-Fu Lu and Liang Yao and Hua Li},
keywords = {High-order manifold regularization, Robust affinity matrices, Multi-view subspace clustering, Weighted TNN},
abstract = {Multi-view subspace clustering achieves impressive performance for high-dimensional data. However, many of these models do not sufficiently mine the intrinsic information among samples and consider the robustness problem of the affinity matrices, resulting in the degradation of clustering performance. To address these problems, we propose a novel high-order manifold regularized multi-view subspace clustering with robust affinity matrices and a weighted tensor nuclear norm (TNN) model (termed HMRMSC) to characterize real-world data. Specifically, all the similarity matrices of different views are first stacked into a third-order tensor. However, the constructed tensor may contain an additional inter-class representation since the data are usually noisy. Then, we use a technique similar to tensor principal component analysis (TPCA) to obtain a more robust similarity tensor, which is constrained by the so-called weighted TNN since the original TNN treats each singular value equally and usually considers no prior information of singular values. In addition, a high-order manifold regularized term is also added to utilize the manifold information of data. Finally, all the steps are unified into a framework, which is resolved by the augmented Lagrange multiplier (ALM) method. Experimental results on six representative datasets show that our model outperforms several state-of-the-art counterparts.}
}
```

