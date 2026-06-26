# Mini Hybrid Switched System（混合切换系统）

一套**从零构建、零外部依赖的纯 C 语言实现**，涵盖混合与切换动力系统理论。混合系统将连续动力学（微分方程）与离散模态切换（自动机、守卫条件、重置映射）相结合，建模对象从弹跳球、恒温器到先进电力电子、网络化控制和自主系统。

## 子模块

| 子模块 | 主题 | 对应课程 |
|--------|------|----------|
| [mini-dwell-time-analysis](mini-dwell-time-analysis/) | 驻留时间稳定性、平均驻留时间、基于 LMI 的计算、公共/多重李雅普诺夫函数、慢切换 | MIT 6.241J, Stanford AA203 |
| [mini-event-triggered-control](mini-event-triggered-control/) | 事件触发采样、自触发控制、Zeno 行为检测、事件间隔分析、相对/绝对阈值 | MIT 6.241J, UC Santa Barbara (Tabuada) |
| [mini-hybrid-automata](mini-hybrid-automata/) | 混合自动机建模（Henzinger）、前向/后向可达性、障碍证书、混合仿真、Zeno 分析 | MIT 6.841, Stanford CS359, Berkeley EECS 291E |
| [mini-impulsive-system](mini-impulsive-system/) | 脉冲动力学、跳变映射、B-等价方法、脉冲镇定、基于李雅普诺夫的脉冲控制 | MIT 6.241J, Caltech CDS140 |
| [mini-piecewise-affine-system](mini-piecewise-affine-system/) | 分段仿射系统辨识、多面体几何、PWA 仿真/可达性、显式模型预测控制 | ETH 227-0216, MIT 6.832 |
| [mini-reset-control-system](mini-reset-control-system/) | Clegg 积分器、重置带、一阶重置元件（FORE）、突破 Bode 积分约束、混合回路分析 | Cambridge 4F2, MIT 6.241J |
| [mini-supervisory-control](mini-supervisory-control/) | 离散事件系统自动机、Ramadge-Wonham 监督理论、最大可控子语言、非阻塞性、模块化监督 | MIT 6.241J, Cambridge 4F3 |
| [mini-switched-stability](mini-switched-stability/) | 公共李雅普诺夫函数、多重李雅普诺夫函数（Branicky）、李代数条件、慢切换、李雅普诺夫-Metzler 不等式 | Berkeley EECS 291E, MIT 6.241J |

## 设计哲学

- **零外部依赖** — 纯 C 语言（C99/C11），仅使用 `libc` 和 `libm`
- **模块自包含** — 每个目录独立包含 `Makefile`、`include/`、`src/`、`examples/`、`demos/`、`tests/`
- **理论到代码的映射** — 每个模块均包含 `docs/` 目录，提供课程对齐说明与知识层级文档
- **实用演示** — 包含弹跳球、恒温器、DC-DC 转换器、车辆编队、制造单元等完整演示

## 构建

每个子模块独立构建。进入子模块目录后执行：

```bash
cd mini-dwell-time-analysis
make all    # 构建全部目标
make test   # 运行测试
```

或构建顶层统一框架：

```bash
make          # 构建库、测试和示例
make test     # 运行所有顶层测试
make test-all # 运行顶层测试 + 所有子模块测试
```

需要 **GCC** 和 **GNU Make**。

## 项目结构

```
mini-hybrid-switched-system/
├── include/                          # 顶层统一 HSS 框架头文件
├── src/                              # 顶层统一 HSS 框架源文件
├── demos/                            # 演示程序
├── examples/                         # 使用示例
├── tests/                            # 顶层测试套件（48 个测试）
├── benches/                          # 性能基准测试
├── docs/                             # 课程对齐与理论文档
├── mini-dwell-time-analysis/         # 驻留时间稳定性分析与 LMI 计算
├── mini-event-triggered-control/     # 事件触发与自触发控制
├── mini-hybrid-automata/             # 混合自动机建模与可达性分析
├── mini-impulsive-system/            # 脉冲动力学与镇定
├── mini-piecewise-affine-system/     # 分段仿射系统辨识
├── mini-reset-control-system/        # 重置控制（Clegg 积分器、FORE）
├── mini-supervisory-control/         # 监督控制（Ramadge-Wonham）
├── mini-switched-stability/          # 切换稳定性（CLF、MLF、李代数）
└── Makefile                          # 顶层构建系统
```

## 许可证

MIT
