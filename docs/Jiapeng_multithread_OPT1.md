# Jiapeng_multithread_OPT1

本分支基于 `Jiapeng_optimization` 的 commit `a97a36451e97666d471bdcceef1790a6dc5c4a90`，包含 P1 图遍历优化候选，并将 I/O 默认配置固定为 `--decode-threads 0 --encode-threads 0`。

## 默认配置

`--threads` 仍控制计算 worker。无论设为 1、4 或 32，省略 I/O 参数时都采用 `decode=0、encode=0`，不额外创建 HTS codec worker 池。输入解码和 BAM 压缩仍正常执行，协调线程和 BAM writer 仍然存在；这不是进程线程总数的硬上限。

例如以下命令自动使用固定 I/O 默认值：

```sh
ExpansionHunter --reads sample.bam --reference reference.fa \
  --variant-catalog catalog.json --output-prefix sample \
  --analysis-mode streaming --threads 32
```

显式数字仍可覆盖默认值。使用 `--decode-threads auto --encode-threads auto` 可恢复原来的自动规模：解码 `min(--threads, 12)`，编码 `--threads`。decoder 参数只用于 streaming；seeking 省略该参数仍可正常运行，保留原有显式数字参数校验。

## P1 改动

- mapped `Graph::hasEdge` 直接查找已校验边界的邻接数组。
- `PathOperations` 与 DAG 遍历使用局部借用的邻接视图，减少共享 archive owner 的引用计数操作；调用期间 Graph 和映射继续由 analyzer 持有。
- 保留原有 owning 邻接 API 的生命周期、邻接顺序、节点校验及自定义 Graph 类型行为。没有增加 Graph 成员、索引副本或线程缓冲，也没有调整队列容量。

缓存格式没有改变，但源码指纹会改变。请为此版本重新准备兼容缓存，不要将此前二进制的缓存直接复用。

## 验证记录与限制

以下性能记录来自修改默认值之前的 P1 二进制，均显式指定 I/O 参数；不能将它们标记为此次重新编译的二进制重新测速结果。环境为 macOS arm64、10 个逻辑 CPU、AppleClang 21、Release。使用仓库 example，以及重复同一坐标、使用不同位点和变异 ID 的 dense32/dense128 合成 catalog。32 worker 是超额订阅测试，不是真实 32 核 WGS、NUMA 或常驻多 sample 性能验证。

P1 完整矩阵及重复实验共 104 次 EH 分析，配对输出全部通过；另有 48 项跨线程比较通过。VCF/JSON 逐字节相同；BAM 完整解压 header 和二进制 record multiset 相同，保留所有字段、tag 类型、浮点位及重复次数，仅忽略记录顺序与 BGZF 压缩表示。

32 worker、dense128、streaming、prepared cache 的结果如下。每种配置预热一对后交替测量三对；缓存准备和结果比较不计入时间。

| 显式 I/O 配置 | Read processing 中位秒数，基线→P1 | 整进程峰值 RSS 中位 MiB，基线→P1 |
|---|---:|---:|
| decode=0, encode=0 | 6.607 → 6.299 | 51.094 → 51.828 |
| decode=auto, encode=auto | 6.707 → 6.700 | 52.719 → 54.422 |

固定 I/O 的阶段耗时中位数下降约 4.65%，但逐对变化为 +8.01%、−4.65%、−0.98%，不能声称稳定获得此幅度加速。自动 I/O 的阶段耗时中位数基本不变。两组配置均三对 RSS 全部上升；**严格 RAM 不增加条件未通过**。默认配置调整并不改变这项结论。

共享 mmap 图的独立比对压力测试中，32 worker 阶段时间中位数为 14.230→1.829 秒（约 7.78 倍吞吐），但不能外推为完整程序/WGS 收益。owned 图对照在 32 worker 下反而慢约 2.69%；mapped32 的 RSS 也在三对运行中全部上升。没有队列峰值遥测，不能确定将 RSS 差值归因于输出积压或工作区容量。

修改默认值后的单独检查：

- 重新编译成功；CLI `--help` 显示两个 I/O 默认值为 0。
- 7/7 参数测试通过，覆盖省略参数、显式数字、`auto`、1/4/32 计算线程、独立覆盖及 seeking 兼容行为。
- dense32、32 worker、seeking/streaming、owned/cache 四种配置下，分别运行原 P1 显式 0/0、新版本省略 I/O 参数、新版本显式 0/0，共 12 次分析、8 项完整输出比较，全部一致。streaming 日志均为 `compute=32, decode=0, encode=0`。
- GraphTest 25/25 通过。完整 CTest 为 30/31 个测试程序通过；UnitTests 为 120/121 用例通过。唯一失败仍是基线已有的 `PreparedCatalog.WholeCatalogRoundTripAndIndependentHandleLifetime`：JSON 标签数组的顺序不同。
- 先前 P1 的 ASan/UBSan focused 图测试 24/24 通过；完整图比对测试遇到既有 DAG vector container-overflow，已在基线与 P1 中用相同 owned 图复现。没有屏蔽该问题，也未将完整 sanitizer 套件记为通过。

本分支供继续验证与优化；未叠加 P2/P3，没有声称真实 32 核性能或 RAM 验收完成。
