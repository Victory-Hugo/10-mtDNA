cat /home/luolintao/schmutzi/script/0_提取MT_error.log|parallel -j 8 \
    /home/luolintao/schmutzi/script/1.5_BAM添加样本信息_核心代码.sh {} 