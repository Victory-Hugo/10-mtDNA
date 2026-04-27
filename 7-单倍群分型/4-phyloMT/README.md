# phyloMT 使用教程

phyloMT 来自Haplogrep3的单倍群分型模块，并采用`python3`实现，当前版本已**完全脱离Java运行环境的依赖**。

## 1. 当前正确使用方式


请直接在项目根目录下使用：

```bash
python3 phyloMT.py --tree ... --input ... --output ...
```


## 2. 查看帮助

查看主帮助：

```bash
python3 phyloMT.py --help
```

当前只有一个命令接口，`--help` 会显示全部分类参数。

## 3. 树资源说明

当前版本直接使用 `data/trees/` 下已有树资源，运行时不再依赖 `Config.yaml`，也不再包含树安装功能。

先查看当前可用树：

```bash
python3 phyloMT.py --list-trees
```


## 4. 核心分类命令

通用格式：

```bash
python3 phyloMT.py \
  --tree <树ID> \
  --input <输入文件> \
  --output <输出文件>
```

当前支持输入格式：

- `HSD`
- `VCF`
- `FASTA`


