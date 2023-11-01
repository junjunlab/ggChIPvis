# ggChIPvis <img src="man/ggChIPvis-logo.png" align="right" height="200" />

<!-- badges: start -->

**ggChIPvis** is designed to re-produce the ***Profile*** and ***Heatmap*** plots uisng
**ggplot2** package with more parameter controls and settings. **ggChIPvis** can accept
the data from **EnrichedHeatmap::normalizeToMatrix**, **ChIPseeker::getTagMatrix** and deeptools'
**computeMatrix** output data which allow you re-draw the plot in R session with other
graphic modifications.

<!-- badges: end -->

## Installation

You can install the development version of ggChIPvis like so:

``` r
# install.packages("devtools")
devtools::install_github("junjunlab/ggChIPvis")

# or
remotes::install_github("junjunlab/ggChIPvis")
```

---

## Citation

> Jun Zhang (2023). *ggChIPvis: Profile and Heatmap Visualization by Using ggplot2.*  https://github.com/junjunlab/ggChIPvis

---

## News

- 2023/10/26  (Changing **grid.xaxis2/grid.yaxis2/grid.colorkey** into grobs.)
- 2023/10/24  (Adding **ChipHeatmap** and **multiHeatmap** for visualization based on ***grid*** systerm.)
- 2023/10/20  (Setting **independent** arg in **if** loop.)
- 2023/10/18  (Add **group.split** in retriveData and **rowgroup.order** in ChipVis.)
---

## Documentation

> ***[https://junjunlab.github.io/ggChIPvis-manual/](https://junjunlab.github.io/ggChIPvis-manual/)***

## Related blogs

> [ChIP-seq 可视化探索](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247510110&idx=1&sn=5b8a12b864509d65c243b61d080350a9&chksm=c184922ff6f31b3943a3681def268998388a3a930b2c384f3d8ea33ff17a437cc4513fc1b286&token=353264504&lang=zh_CN#rd)
> [ggChIPvis 可视化基因组富集信号](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247510122&idx=1&sn=f06cfbecb15ee4f133c9b3d49791743a&chksm=c184921bf6f31b0d471e37ef63fdc4c7caf63df2a6598bd1f3c04d9b1505e2e103ef4f60eb66&token=353264504&lang=zh_CN#rd)
> [ChIP-seq 绘图设计](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247510151&idx=1&sn=de70b667d6eb7f581a68fd3fc231b1fc&chksm=c18492f6f6f31be05a57cb681d802fb813e36ac67f3d70fb05d26a6ea1aeee14bd17a91844e7&token=353264504&lang=zh_CN#rd)
> [ChIP-seq 绘图设计<<续>>](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247510209&idx=1&sn=ea595989806aeaf116fab44f808d295e&chksm=c18492b0f6f31ba6832f666ee1734d258dcb0ad852a2982443668f735e46e86f2a9456a1059d&token=353264504&lang=zh_CN#rd)

