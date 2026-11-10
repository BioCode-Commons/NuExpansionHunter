# NuExpansionHunter 索引压缩候选验证报告

生成时间（UTC）：2026-11-10T09:23:11.195604+00:00
分支：Jiapeng_optimization；原始基线：cdc6dd9fc6ef0a3bc2fc24f77f9fe9f6597a4f8c

## 当前状态

**严格不退化准入尚未证明：端到端观测存在波动**

部分端到端中位数出现正向变化，但逐对方向不一致或重复证据不足；不能仅凭任意小的正数宣称已确认稳定退化，也不能据此证明严格不退化。当前结果应如实保留为有限本地观察，不能承诺目标 WGS 不变慢或声称已完成该目标场景的性能验收。

用户目标是同一台机器并行处理更多样本，并严格保留输出。32 核、118 GB RAM、5,577,376 位点真实 WGS 尚未验收；没有实测 samples/hour、目标峰值内存或最大样本并发数。

## 候选改动范围

候选版本：SoA + narrow buckets + forceinline；PathView 延迟读取 count

源码实现及已测输出检查的完成情况，与严格性能保证分开报告。性能准入未证明不妨碍提交代码供审查，但不得把它写成目标 WGS 已获保证。

第一阶段包含以下改动：

- `--decode-threads auto|N` 与 `--encode-threads auto|N` 分别设置读解码和 BAM 压缩线程。auto 保留原默认规模；`--threads` 仍表示计算 worker，并非进程全部线程的硬上限。
- `--bamlet-queue-mb N` 可选限制待写 BAM 记录的字节预算；默认 0 保留无界队列。正值使用背压，不丢弃记录，单条超大记录仍可推进。它不是整个进程 RAM 上限，启用后需单独测试吞吐量。
- JSON 按位点序列化，保留字段与位点顺序、数值、转义、缩进和重复 ID 最后覆盖行为，避免构建完整输出 DOM；SampleFindings 和轻量排序表仍全量驻留。

第一阶段独立 JSON 实验（100,000 位点、每位点一个变异）中，writer 时间中位数 0.762234→0.565990 秒；包含输入 fixture 的进程峰值 604.17→233.08 MiB。该实验不是整个 WGS 的 RAM 或总耗时预测，原始记录随验证包附上。

本轮候选将索引拆为连续 key、payload 和完整路径 count 数组，省略 MiniPath 的隐式 count；唯一 key 数不超过 65,535 时使用 16 位桶，否则使用 32 位桶。完整 key、路径、节点、顺序和重复次数保留。另包含批量序列化和查询内联。这些是候选实现说明，不表示性能验收完成。

## 端到端性能：准入依据

对照是已含第一阶段改动、尚未压缩索引的快照。以下主实验使用 4 个计算线程，每配置 1 次预热 pair + 5 次测量 pair，交替执行顺序。计时包括启动、catalog 加载、分析和输出，不包括缓存准备及结果比较。微基准结果不能替代本节。

| 配置 | 对照中位秒数 | 候选中位秒数 | 中位数之比变化 | 逐对变化中位数 | 配对更慢/更快 |
|---|---:|---:|---:|---:|---:|
| seeking-owned | 1.5371 | 1.5314 | -0.37% | +0.41% | 3/2 |
| seeking-cached | 1.6210 | 1.6272 | +0.38% | -0.67% | 2/3 |
| streaming-owned | 1.6400 | 1.6604 | +1.24% | +1.30% | 4/1 |
| streaming-cached | 1.7369 | 1.6853 | -2.97% | -2.97% | 0/5 |

补充端到端实验：1 计算线程；32 次运行；状态 complete。
证据目录：`paired-end-to-end-lazy-t1`

| 配置 | 对照中位秒数 | 候选中位秒数 | 中位数之比变化 | 逐对变化中位数 | 配对更慢/更快 |
|---|---:|---:|---:|---:|---:|
| seeking-owned | 5.4051 | 5.3983 | -0.13% | -0.13% | 1/2 |
| seeking-cached | 5.3980 | 5.3859 | -0.22% | -0.22% | 1/2 |
| streaming-owned | 5.2687 | 5.2558 | -0.25% | -0.29% | 1/2 |
| streaming-cached | 5.2509 | 5.2482 | -0.05% | -0.08% | 1/2 |

正值表示该统计量下候选更慢。小幅正向中位数变化可能来自测量波动，不能直接写成稳定退化；需同时看配对方向及重复次数。报告保守区分“尚未证明严格不退化”和“所测配对持续更慢”，没有调整数据或放宽阈值来过关。

所选记录：48 次运行、20 个测量 pair。No OS page-cache clearing; one warmup pair before each configuration's measured pairs.
实验范围：Small synthetic dense32 example; this is not a whole-genome, multi-sample, or 32-fold genome-wide benchmark.

## 索引微基准：仅用于局部解释

记录说明：50 deterministic synthetic indexes; 72 hot-512-query cells; 5 alternating process pairs

单位操作时间的候选/对照中位数比范围 0.604–0.987；所列 72 个组合中 0 个比值大于 1。

6 项构建计时的中位比值范围 0.963–1.015；wide14 构建约慢 1.54%，该局部差异保留在原始数据中。

该测试反复查询固定 key 集合，分别统计命中、未命中、计数和取回首条路径。它不包含实际 read 比对中全部路径扩展、其它数据结构与多线程之间的缓存影响，也不是全 catalog、冷缓存或多样本测试。不得把局部加速写成 WGS 总加速。

## 归档容量

| 数据 | 对照字节 | 候选字节 | 减少 |
|---|---:|---:|---:|
| example | 79,645 | 55,369 | 30.48% |
| dense32 | 2,536,341 | 1,759,633 | 30.62% |
| 合成索引归档 | 5,360,343 | 4,658,347 | 13.10% |

这里统计归档容量，不是整个样本的峰值 RSS 降幅。相同文件的 mmap 页可跨进程共享，应结合 PSS/USS 和整机内存；不能将 RSS 简单相加，也不能把小数据压缩比例直接乘到 558 万位点 WGS。

## 输出一致性与边界检查

口径：VCF/JSON 逐字节一致；BAM 完整 header、记录、标签及重复次数一致，仅允许记录顺序不同。

- 所选完整回归记录：40 次，其中 40 次标记输出一致；要求至少 40 次。
- 配对计时的 48 次四线程及 32 次单线程运行，也均在计时后完成完整输出比较并通过。
- 最终 EH 单元测试日志：119/120 通过；失败项：PreparedCatalog.WholeCatalogRoundTripAndIndependentHandleLifetime。
- GraphTools 的 30 个测试程序全部通过；EH 的既有失败单独记录。
- 新增 9 项索引测试覆盖有序路径/重复、空索引、16-mer 完整 key、位置和节点 ID 越界回退、65,535/65,536 桶宽切换、宽索引、mapped 生命周期、旧格式及 checksum 拒绝。
- ASan/UBSan 日志记录验证 50 个索引。仅变更的 KmerIndex.cpp 与独立 probe 启用 sanitizer，其它对象仍按普通 Release 链接；未做 LeakSanitizer 检查。

上述 PreparedCatalog 测试是该 macOS 环境的既有 graph edge 标签顺序差异；没有把它写成全部通过。原始基线对照证据应随验证包保存。

索引语义差分：50 个索引、25 种图；规范文本 19902304 字节，SHA256 `c6b665dec8314fbcab93f349b69193634c3a55b61b12387d578a1efc1d87ad37`。范围：all keys, ordered encode, ordered paths and sequences, counts, node/edge unique counts, copies, moves, mmap and mapped lifetime, membership/absent PathView exceptions。

有限测试支持已测范围内的一致性，不能无条件保证所有未见输入或平台。

## 缓存兼容与交付状态

候选使用 v3 归档格式，旧缓存不能直接复用。试验该候选时使用新文件名重新准备缓存；各样本共享同一匹配文件。缓存校验规则没有放宽。代码可提交供审查和继续试验；本地性能验收状态不得改写为真实 WGS 严格不退化的保证。

## 证据来源与二进制对应关系

本地环境：macOS arm64，Apple clang 21.0.0，Release；Boost 1.73.0，HTSlib 1.10.2（仓库已有补丁）。双方使用相同依赖与优化参数；未在 Linux、32 核或多样本生产环境验证。

- paired_results：`work/index-validation/paired-end-to-end-lazy/paired-results.json`
- correctness_results：`work/index-validation/end-to-end-lazy/results.json`
- micro_results：`work/index-validation/paired-lazy-five/final-summary.json`
- ctest_log：`work/index-validation/ctest-lazy-final.log`
- baseline_ctest_log：`test-baseline.log`
- sanitizer_log：`work/index-validation/sanitized-lazy-probe.log`
- sanitizer_build_log：`work/index-validation/build-sanitizer-lazy.log`
- baseline_binary：当前 SHA256 `079ad0a47319f8caf169cfde52886b0341c322a9c63df9613c462c9560f07f9c`；计时 SHA256 `079ad0a47319f8caf169cfde52886b0341c322a9c63df9613c462c9560f07f9c`；匹配=True。
- candidate_binary：当前 SHA256 `419d88a0837ff56555d6e3d5e9a9ddfc81faf42dfe8427538ecea0fdd658833b`；计时 SHA256 `419d88a0837ff56555d6e3d5e9a9ddfc81faf42dfe8427538ecea0fdd658833b`；匹配=True。

报告读取显式选择的结果文件；后续重编译或更换实验目录时，必须更新配置并重新生成，不能沿用旧报告。
