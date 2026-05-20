# 研究背景
AMOVA（Analysis of Molecular Variance）是一种用于分析分子数据中遗传变异的统计方法。
目前最正宗的计算软件是Arlequin，但其界面较为复杂，且不支持命令行操作。
因此，我需要一个基于Python的AMOVA计算工具,CLI工具。
# 研究目标
原始的Arlequin软件`/mnt/e/Scientifc_software/arlecore_linux`，其中有各种核心的算法。你需要将其进行重构，并形成1套pipline。
# 代码风格示例
`/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/5-经典群体遗传学计算`。
# 示例数据
## 输入
- `/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/15-AMOVA/input/Altaic_本次研究_chrM.vcf.gz`。
- 分组文件：`/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/15-AMOVA/input/1-group.tsv`;`/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/15-AMOVA/input/2-group.tsv`。

## 预期输出
`/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/15-AMOVA/tmp/预期产出.md`：用arlequin软件计算的AMOVA结果。
