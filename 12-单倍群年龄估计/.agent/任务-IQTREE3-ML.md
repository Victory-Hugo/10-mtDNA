# IQ-TREE 3 ML 单倍群年龄估计规划（修订版 v2）

## 背景与目标

### 为什么需要 ML 方法

ρ 方法已对全部 5,439 个单倍群完成年龄估计，但存在已知统计弱点：

| 样本量区间 | 单倍群数 | ρ 方法可信度 |
|-----------|---------|------------|
| n < 2     | 183     | 无 CI，完全不可信 |
| 2 ≤ n < 10 | 3,082  | 方差极大，低可信 |
| 10 ≤ n < 50 | 1,228 | 中等可信 |
| n ≥ 50    | 948     | 统计稳健 |

**3,265 个小样本单倍群（60% 总量）的 ρ 估计不可靠。**

IQ-TREE 3 方法使用的信息来源不同：3,019 条叶节点祖先重建序列，每个节点的年龄估计**完全不依赖其现代后代样本数量**，天然绕开小样本问题。

### 目标

1. 对全部 5,439 个单倍群提供独立于 ρ 的 ML 年龄估计
2. 重点为 n<10 的 3,265 个小样本单倍群提供可信年龄
3. 与 ρ 结果交叉验证，大样本节点两者应高度吻合

---

## 数据现状

| 文件 | 路径 | 说明 |
|------|------|------|
| 祖先序列 | `data/ancestral_sequences.fasta` | 5,439 条，含全部节点（叶+内部）|
| PhyloTree JSON | `data/phylotree_index_withacc.json` | 单倍群图（含 parent/lineage/mutations）|
| rCRS | `data/rcrs.fasta` | 16,569 bp 参考序列 |
| 祖先样本列表 | `output/stage_1/intermediate/ancestor_samples.txt` | 5,439 行 |
| ρ 结果 | `output/stage_1/results/rho_dating.tsv` | 全量结果，用于校准锚点 |
| 单倍群图 | `output/stage_1/intermediate/haplogroup_graph.json` | parent/children 结构 |

### 工具环境

- **IQ-TREE 3.0.1**：`/home/luolintao/miniconda3/envs/BigLin/bin/iqtree3`
- **halign4**：`bin/halign4`
- **conda 环境**：`/home/luolintao/miniconda3/envs/BigLin`（含 ete3、biopython、numpy、pandas）

---

## 方法概述与已知问题（v2 修订）

```
Step A: s08_prepare_ml_inputs.py
    ancestral_sequences.fasta
    + haplogroup_graph.json（识别叶节点）
    → ancestors_sanitized_leaves.fasta（3,019 条叶节点序列，净化名称）
    → name_sanitize_map.tsv（5,439 行全量映射，供 s10 反映射）

Step B: halign4 比对
    ancestors_sanitized_leaves.fasta
    → ancestors_aligned_leaves.fasta（统一长度 16,634 bp）

Step C: s09_phylotree_to_newick.py
    haplogroup_graph.json → phylotree_fixed.nwk
    （5,439 节点，3,019 叶，2,420 内部；叶与 FASTA 完全匹配）

Step D: IQ-TREE 3
    ancestors_aligned_leaves.fasta + phylotree_fixed.nwk
    -te 固定拓扑，-m GTR+G，-T AUTO
    → anc_tree.treefile（含 ML 枝长）
    ⚠️ 已知行为：IQ-TREE 合并 1,059 个单子节点、删去根节点标签

Step E: s10_ml_dating.py（局部路径插值校准）
    anc_tree.treefile + name_sanitize_map.tsv
    + rho_dating.tsv（校准锚点） + haplogroup_graph.json（单子节点处理）
    → ml_dating.tsv（5,439 行 ML 年龄 + CI）
    → ml_vs_rho_comparison.tsv（比较）
```

### 为什么只用 3,019 条叶序列（不用全部 5,439）？

IQ-TREE `-te` 要求 FASTA 序列 = 树叶节点。内部节点序列不能直接输入。
叶序列携带了其完整祖先链所有定义突变的累积信息，ML 可以从叶序列正确推断所有内部枝长。
全部 5,439 个节点（包括 1,060 个 IQ-TREE 折叠节点）的年龄均可从校准后的树中恢复。

---

## 关键技术问题与解决方案（v2 新增）

### 问题 1：IQ-TREE 删去根节点标签

**原因**：IQ-TREE 内部将有根树当作无根树优化，输出时根标签 `mt-MRCA_RSRS` 丢失。

**解决**：s10 中用 ete3 的 `set_outgroup("L0")` 重新定根（`mt-MRCA_RSRS` 的两个子节点是 `L0` 和 `L1'2'3'4'5'6'`），然后将新根命名为 `mt-MRCA_RSRS`。

### 问题 2：1,059 个单子节点被 IQ-TREE 合并

**原因**：单子节点（如 A1 → A1a 的直链）上下两段枝长不可独立估计，IQ-TREE 将其合并为一条枝。

**解决**：s10 中对这些节点做**父子中点插值**：
```
age_V = (age_parent + age_child) / 2
```
父子年龄优先取 ML 树的结果，无 ML 结果则用 ρ 兜底。

### 问题 3：原校准公式错误（最关键）

**原公式**（已废弃）：`age_V = root_age × (1 - D_V / D_mean)`

**错误原因**：此公式假设所有叶节点在同一时间点（如现代样本）。但我们的叶序列代表不同历史时间的祖先（从 ~0 kya 到 ~130 kya），叶节点到根的距离天然不同（CV = 0.4977），与速率无关，是时间跨度异质性的必然结果。强行套用此公式导致 H 节点偏差高达 65 kya。

### 解决：局部路径插值（Local Clock Interpolation）

对 ML 树中每个节点 V：

1. 向上找最近的**高可信 ρ 锚点祖先 A**（n ≥ 10）
2. 找 V 子树中所有**高可信 ρ 锚点后代 D_i**（n ≥ 10）
3. 用 ML 枝长做比例插值：

```
f_i = D_{A→V}^{ML} / D_{A→D_i}^{ML}   （V在A到D_i路径上的比例位置）

age_V (from D_i) = age_A − f_i × (age_A − age_D_i)

age_V_final = median(age_V from all D_i)
```

4. CI 同比例传播：
```
ci_lower_V = ci_lower_A − f × (ci_lower_A − ci_lower_D)
ci_upper_V = ci_upper_A − f × (ci_upper_A − ci_upper_D)
```

**优势**：
- ML 枝长只用于确定**相对位置**，不假设绝对速率
- 绝对时间由 ρ 锚点（高可信节点）定标
- 规避了全局速率假设和叶时间异质性问题
- CI 从 ρ 自然传播，无需 bootstrap

### 关于 bootstrap（`-b 100` / `-bb 1000`）

**结论：不添加 bootstrap。**

原因：
- `-te` 固定拓扑，bootstrap 不提供有用的拓扑支持值
- 枝长 CI 需要用 `-wbt`（写出所有 bootstrap 树）再后处理，流程复杂
- 我们的校准方法（局部路径插值）中，枝长仅作为插值权重使用
- **年龄 CI 的主要来源是 ρ 的 CI**，不是 ML 枝长的不确定性
- 标准 bootstrap 需额外 ~1.6 小时，ultrafast bootstrap 需 ~3-6 分钟，均不划算

---

## 脚本规划（v2）

### 已完成并验证通过

| 脚本 | 状态 | 说明 |
|------|------|------|
| `s08_prepare_ml_inputs.py` | ✅ 完成 | 净化名称 + 提取叶节点（3,019条） |
| `s09_phylotree_to_newick.py` | ✅ 完成 | 生成固定拓扑 Newick |
| IQ-TREE 3 运行 | ✅ 完成（58 秒）| anc_tree.treefile 已就绪 |

### s10_ml_dating.py（v2，局部插值校准）

**输入**：
- `output/stage_2/intermediate/ml/anc_tree.treefile`
- `output/stage_2/intermediate/ml/name_sanitize_map.tsv`
- `output/stage_1/results/rho_dating.tsv`（校准锚点 + 兜底）
- `output/stage_1/intermediate/haplogroup_graph.json`（单子节点处理）

**算法**：

```python
# 1. 重新定根
t.set_outgroup("L0")
root_node.name = "mt-MRCA_RSRS"

# 2. 预计算根到所有节点的 ML 距离（O(n)）
dist_from_root = {id(n): n.get_distance(root) for n in t.traverse()}

# 3. 构建后序子树校准叶集合（O(n)，避免 O(n²) 遍历）
calib_leaves_below = {id(leaf): [entry] for leaf in calibrated_leaves}
# 向上合并

# 4. 预序构建最近校准祖先索引（O(n)）
nearest_calib_anc = {}  # node_id → nearest calibrated ancestor node

# 5. 对每个节点做局部路径插值

# 6. 处理 1,059 个单子节点（父子中点）

# 7. 输出
```

**输出**：

`output/stage_2/results/ml_dating.tsv`（5,439 行）：

| 列名 | 说明 |
|------|------|
| `haplogroup` | 单倍群名 |
| `ml_age_years` | ML 校准年龄（年） |
| `ml_ci95_lower_years` | 95% CI 下界（年） |
| `ml_ci95_upper_years` | 95% CI 上界（年） |
| `ml_age_kya` | ML 年龄（kya） |
| `ml_ci95_lower_kya` | CI 下界（kya） |
| `ml_ci95_upper_kya` | CI 上界（kya） |
| `n_samples` | ρ 样本数（0=仅内部节点）|
| `calibration_method` | 校准方法标记 |

`output/stage_2/results/ml_vs_rho_comparison.tsv`（5,439 行）：
- 包含 ρ 和 ML 年龄、CI、差值、方法标记

---

## Pipeline 整合

`pipe/run_ml.sh` → `conf/stage_2_ml.yaml`，步骤：
1. s08（净化+叶节点提取）
2. halign4（比对）
3. s09（Newick 生成）
4. IQ-TREE 3（固定拓扑 ML，GTR+G）
5. s10（局部路径插值校准）

---

## 预期结果与验证（修订）

### 关键节点验证（预期 ML vs ρ 偏差）

| 节点 | ρ 年龄（kya）| 可接受 ML 偏差 | 预期校准方法 |
|------|-------------|--------------|------------|
| mt-MRCA(RSRS) | 163.9 | 0（根定义）| root_calibration |
| L3 | 64.4 | ±3 kya | rho_direct（n=71953）|
| M | 50.2 | ±3 kya | rho_direct（n=32778）|
| H | 14.7 | ±5 kya | rho_direct（n=4242）|
| D4 | 28.2 | ±5 kya | rho_direct（n=8109）|

大样本节点（n≥50）应直接返回 ρ 结果（`rho_direct` 方法）。
小样本节点（n<10）：若子树内有校准后代则 `local_interpolation`，否则 `rho_fallback_no_desc`（ρ 兜底）。

### 校准方法分布（实际结果 v2）

| 方法 | 实际节点数 | 说明 |
|------|----------|------|
| `rho_direct` | 1,677 | n≥10 且非单子节点的 ML 树节点 |
| `single_child_interpolation` | 1,059 | 单子节点父子中点（含 480 个 n≥10 单子节点）|
| `rho_fallback_no_desc` | 2,702 | n<10 无校准后代，ρ 兜底 |
| `root_calibration` | 1 | mt-MRCA(RSRS) 固定锚点 |
| `local_interpolation` | 0 | 见下方说明 |

**`local_interpolation = 0` 说明**：

ML 树中 3,019 个叶节点（终端单倍群）天然无后代，无法路径插值。195 个非校准内部节点的子树也全为 n<10 稀有节点，无校准后代可用。这是方法的结构性限制：路径插值只能用于有校准后代的节点，终端小样本单倍群无此条件，ρ 兜底是正确行为。

**额外已知 bug（已修复）**：IQ-TREE 将原始 Newick 的根节点标签（`mt-MRCA_RSRS`）改为 `A`（haplogroup A 的节点名）作为 treefile 根，ete3 `set_outgroup` 重定根后创建匿名节点代替 haplogroup A 的位置。s10 通过检测匿名内部节点（子节点与原 treefile 根节点子节点有重叠）自动恢复正确标签。

---

## 注意事项（v2 更新）

1. **名称净化一致性**：s08/s09/s10 使用相同净化规则（已验证无碰撞）
2. **重新定根**：s10 必须在读取 treefile 后立即重新定根，否则距离计算错误
3. **CI 下界夹零**：ρ CI 下界有时为负（小样本），s10 应 clip 到 0
4. **局部速率有效性**：路径插值假设同一线系内速率平稳；跨主要分支（如 L→M→N）时应找到中间校准锚点
5. **单子节点 ≠ 零信息**：这些节点有 ρ 估计，兜底用 ρ 不会产生空洞
6. **IQ-TREE 根标签替换 bug**：IQ-TREE 有时将 treefile 根标签替换为树内某个节点名（本例为 `A`），需在 s10 重定根时检测并恢复匿名节点标签
