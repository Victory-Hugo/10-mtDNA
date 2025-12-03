cat /home/luolintao/schmutzi/script/input.list|parallel -j 16 \
    /home/luolintao/schmutzi/script/0_提取线粒体DNAreads_核心代码.sh --inputfile {}\
    --output_dir /home/luolintao/schmutzi/INPUT/处理汇总/ \
    --errlorlog /home/luolintao/schmutzi/script/0_提取MT_error.log
