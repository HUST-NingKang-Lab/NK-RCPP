#!/usr/bin/env python3
"""
修正负向对接的标签错位 bug。

原 bug: 生成负向输入时用 `comp{i:02d}` (i 从 0 开始) => comp00..comp09
        而拆分批次时用 `comp{i:02d}` (i 从 1 开始) => comp01..comp10
        结果: 批次 comp01 装的是 2 号肽, comp02 装的是 3 号肽, ..., comp09 装的是 10 号肽,
             1 号肽漏掉, 批次 comp10 为空。

修正策略 (保持对接结果不变, 仅重命名——计算内容本身是正确的):
  1. 结果目录 `*_compNN` (NN=09..01) 逆序重命名为 `*_comp(NN+1)`
  2. 目录内文件名中的 compNN 同步重命名
  3. 1 号肽 (comp01) 需重新对接 (本脚本只负责重命名)
"""
import shutil
import os
from pathlib import Path

# --- paths ---------------------------------------------------------------
# Repository-relative: the scripts run from anywhere inside a checkout.
# Raw RAPiDock run directories and the competition workbook are not shipped;
# set NK_RCPP_RAW / NK_RCPP_COMP_XLSX (or NK_RCPP_BASE) to point at them.
# See scripts/README.md.
BASE = Path(os.environ.get("NK_RCPP_BASE", Path(__file__).resolve().parent.parent))
RAW = Path(os.environ.get("NK_RCPP_RAW", BASE / "raw"))
DOCKING = RAW / "peptide_docking"
RAPIDOCK = DOCKING / "RAPiDock"
RESULTS_DIR = RAPIDOCK / "results"
# ---------------------------------------------------------------------------

NEG = RESULTS_DIR / "COMP_TOP10_neg"


def main():
    # reverse order to avoid collisions
    for nn in range(9, 0, -1):
        old = f"comp{nn:02d}"
        new = f"comp{nn+1:02d}"
        renamed = 0
        for d in sorted(NEG.glob(f"*_{old}")):
            if not d.is_dir():
                continue
            target = d.parent / (d.name[:-len(old)] + new)
            if target.exists():
                print(f"  [SKIP exists] {d.name} -> {target.name}")
                continue
            d.rename(target)
            # rename inner files
            for f in list(target.iterdir()):
                if old in f.name:
                    f.rename(target / f.name.replace(old, new))
            renamed += 1
        print(f"[{old} -> {new}] renamed {renamed} dirs")


if __name__ == "__main__":
    main()
