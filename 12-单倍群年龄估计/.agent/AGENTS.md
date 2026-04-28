# 环境与软件安装

- 请使用名为`BigLin`的conda环境进行分析。
- 若该`conda`环境中没有安装相关工具，请先安装相关工具。
- `conda`路径是`/home/luolintao/miniconda3/condabin/conda/`。
- 软件安装方式：优先使用`conda`安装。
- 若需计算群体遗传相关指标，**优先使用成熟的生物信息学工具包或者软件（例如pysam, scikit-allel, scikit-learn）等**
- 若不存在相关工具，可以自己编写代码实现相关功能；但是**自己编写代码不作为首选方案**。
- 若拟安装软件与`BigLin`环境中已安装的软件存在版本冲突，则新建一个`conda`环境。



# 编写代码的要求

1. 你需要编写的是清晰的生物信息学**pipeline代码**(使用生信pipeline.skill.md)，注释要清晰，变量命名要有意义（注释要中文）。
2. 所有的python代码均需使用**双模块**(使用生信python双模块.skill.md)。