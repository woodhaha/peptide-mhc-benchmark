# 姐妹篇评估与大纲 (2026-08-03)

## 判断：可独立成篇，但需与主稿内容边界清晰

主稿（Bioinformatics，4,409 词）已包含多等位基因**核心结果**（§3.6 泛化、§3.7 平台、Figure 9、表 4/5）。姐妹篇必须**不重复**这些，只深挖被压缩掉的**完整数据与深入分析**。

## 目标期刊
**PLOS Computational Biology** (IF~4.5) — 记忆里推荐的多等位基因扩展文目标。收深度方法学研究，字数宽松（~6,000-8,000），适合完整呈现多等位基因基准。

## 核心叙事（一句话）
> 一个重平衡 pan 模型能在多大范围内取代 per-allele 面板？结构预测能否突破非规范结合子平台？——通过 5 等位基因 + 15 等位基因 + 时间 holdout + 结构预测器的完整多等位基因基准回答。

## 与主稿的边界
| 内容 | 主稿（已有） | 姐妹篇（新增） |
|------|------------|--------------|
| 重平衡 pan 5 等位基因 5 种子 | ✅ §3.6 | 逐等位基因完整曲线、误差棒细节 |
| pan16 0.819 pooled | ✅ §3.7 | 15 等位基因**逐个**完整表、数据稀缺 tail 分析 |
| Type A/B 平台 | ✅ 表 4 | Type A/B 定义协议、broad-motif retrain 全结果 |
| MUNIS 结构预测 | ✅ 提及 | MUNIS 全人群测试（n=3,112 Type B）完整方法+结果 |
| 时间 holdout | ✅ 结果 | 完整协议、pre-2021 vs post-2022 全表、打破标签循环论证 |
| 盐桥 F10 | ✅ | 结构基础完整呈现 |

## 大纲

### Title
"A Rebalanced Pan-Specific Model and the Non-Canonical Frontier: A Multi-Allele Benchmark of Peptide-MHC-I Binding Prediction"

### Abstract
多等位基因基准：5 等位基因 → 15 等位基因面板、结构预测器对比、时间 holdout 打破标签循环。

### 1 Introduction
- pan-specific 范式（NetMHCpan、MHCflurry、MUNIS）
- 逐等位基因 vs pan 的争议（per-allele 优势是数据量伪影）
- 非规范结合子的真实难度（结构预测能否解决）

### 2 Methods
- 2.1 多等位基因数据（5 等位基因 n=14,488 + 15 等位基因 pan16）
- 2.2 Type A/B 定义与分类协议
- 2.3 重平衡 pan 模型（5 种子）+ per-allele 对照
- 2.4 比较器：MHCflurry、MUNIS（结构）
- 2.5 时间 holdout 协议
- 2.6 评估：蛋白背景负样本、固定阈值、Wilson CI

### 3 Results
- 3.1 重平衡 pan = per-allele 持平（5 等位基因，5 种子）——修正 per-allele>pan 伪影
- 3.2 规范结合子 ~0.98（所有预测器）
- 3.3 非规范（Type B）平台 ~0.65（DFFN/MHCflurry/MUNIS）
- 3.4 Type A 固有难（n=151，broad-motif retrain 不恢复）
- 3.5 硬尾跨预测器抵抗 + 结构单例（QYDPVAALF MUNIS 0.972）
- 3.6 pan16 15 等位基因完整表 + 数据稀缺 tail
- 3.7 时间 holdout 打破标签循环

### 4 Discussion
- 修正伪影的方法学意义（per-allele>pan = benchmark artifact）
- 结构预测平台：挑战"结构解决一切"假设
- pan16 ≈ NetMHCpan 实用价值（简单 DFFN）
- 三阶梯 pipeline 作为假说
- 局限：IEDB 派生、非-IEDB 实验验证待做

### 5 Conclusion
一个重平衡 pan 模型取代 per-allele 面板；规范解决，非规范尾定义前沿。

## 数据与图
- 图 1：重平衡 pan 5 种子 per-allele 曲线
- 图 2：pan16 15 等位基因完整条形图
- 图 3：预测器平台（DFFN/MHCflurry/MUNIS）
- 图 4：时间 holdout（pre-2021 vs post-2022）
- 图 5：Type A/B 分类 + broad-motif retrain
- 表：15 等位基因全表、Type A/B 计数、时间 holdout 全表

## 现有代码/数据
- `src/eval_pan16.py`、`eval_pan_multiseed.py`、`collect_type_ab.py`、`temporal_holdout.py`、`score_deepffn_type_ab.py`
- `06_Logs/type_ab_*.json`、`multiseed_results.json`、`05_Submission/munis_results/`
- 图 9 已有（可复用为姐妹篇图 1）

## 复用
主稿 Figure 9（重平衡 pan 5 种子）可直接作为姐妹篇核心图；findings_panel 的 4 面板（temporal/pan16/P2/plateau）已生成，可拆分使用。
