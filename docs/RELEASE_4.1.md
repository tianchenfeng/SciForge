# SciForge 4.1 · Desktop Edition

面向 Windows 10/11 x64 的科研数据分析工作台。

## 下载与使用

1. 下载 `SciForge-v4.1-desktop-windows-x64.zip`。
2. 将压缩包完整解压到本地文件夹。
3. 双击 `SciForge/SciForge.exe`，无需安装 Python。

请保留 EXE 同目录的 `_internal` 文件夹，发送给同学时直接发送整个 ZIP。

## 本次改动

- 统一浅灰白界面、蓝色操作按钮、侧栏导航和模块卡片。
- 八个页面增加滚动保护；首页卡片、配方约束和代谢组面板按窗口宽度调整排列。
- 修复配方目标值加减按钮重复创建、按钮状态与实际控件不一致的问题。
- 为配方表类型列保留足够宽度，结果标签支持滚动。
- 统计、代谢组和配方输入提前校验无穷大、非有限值及不安全的极端数值；统计分组不能为空。
- 修正统计失败时仍显示“完成”日志的问题。
- 窗口版 EXE 的分析子进程恢复标准输出管道，保留进度和结果回传。
- 发布包仅包含应用资源，不包含本地实验数据、测试输出或 Python 环境。

ANOVA、Tukey、多重比较字母、代谢组建模和其他已有统计计算算法保持原样。

## 修改文件

- 界面：`ui/main_window.py`、`ui/home_page.py`、`ui/theme.py`、`ui/page_feed_formula.py`、`ui/page_metabolomics.py`、`ui/page_stats.py`。
- 输入检查：`modules/module_stats/service.py`、`modules/module_feed_formula/validators.py`、`modules/module_metabolomics/validation.py`；版本记录：`modules/module_metabolomics/service.py`。
- 运行与打包：`app.py`、`common/runtime_check.py`、`build_release.py`。
- 验收：`tests/conftest.py`、`tests/test_desktop_layout.py`、`tests/test_stats_boundaries.py`、`tests/test_module_acceptance.py`、`tests/test_16s_pipeline.py`、`tests/test_local_study_mutations.py`。
- 说明：`README.md`、`.gitignore`、`docs/RELEASE_4.1.md`、`docs/VALIDATION_4.1.md`、`docs/desktop-preview.png`。

## 风险与兼容性

有效的常规输入仍调用原有计算函数。新校验会拒绝以往可能继续计算的空分组、非有限数值或超出安全范围的输入，需要补全分组或换算单位。小窗口中的长内容通过滚动访问；配方表在列宽不足时使用表格本身的水平滚动条。

## 检验范围

完整测试结果、复现命令和边界说明见发行包中的 `VALIDATION_4.1.md`，或仓库 `docs/VALIDATION_4.1.md`。

需要诊断本机运行环境时，可在命令行执行：

```powershell
.\SciForge.exe --self-test --output .\diagnostics
```

诊断会写出 JSON 结果及界面截图。原始 FASTQ 到 QIIME2/PICRUSt2 的外部环境安装不属于此 Windows 包；16S 功能读取已导出的对应输入文件。
