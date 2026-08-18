from pathlib import Path
import shutil
import subprocess
import sys


PROJECT_ROOT = Path(__file__).resolve().parent
APP_NAME = "SciForge"

DIST_DIR = PROJECT_ROOT / "dist"
BUILD_DIR = PROJECT_ROOT / "build"
RELEASE_DIR = PROJECT_ROOT / "release" / APP_NAME


def remove_path(path: Path):
    if path.exists():
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def copy_dir(src: Path, dst: Path):
    if src.exists():
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(src, dst)


def main():
    print("===== 清理旧打包文件 =====")
    remove_path(DIST_DIR)
    remove_path(BUILD_DIR)
    remove_path(RELEASE_DIR)

    print("===== 开始 PyInstaller 打包 =====")

    cmd = [
        sys.executable,
        "-m",
        "PyInstaller",
        "--noconfirm",
        "--onedir",
        "--windowed",
        "--name",
        APP_NAME,
        "app.py",
        "--add-data",
        "assets;assets",
        "--add-data",
        "templates;templates",
        "--add-data",
        "example_data;example_data",
        "--add-data",
        "resources;resources",
        "--add-data",
        "modules;modules",
        "--add-data",
        "ui;ui",
        "--add-data",
        "common;common",
        "--add-data",
        "config;config",
        "--hidden-import",
        "openpyxl",
        "--hidden-import",
        "xlrd",
        "--hidden-import",
        "docx",
        "--hidden-import",
        "scipy",
        "--hidden-import",
        "scipy.optimize",
        "--hidden-import",
        "scipy.spatial",
        "--hidden-import",
        "scipy.stats",
        "--hidden-import",
        "statsmodels",
        "--hidden-import",
        "statsmodels.api",
        "--hidden-import",
        "statsmodels.formula.api",
        "--hidden-import",
        "statsmodels.stats.multicomp",
        "--hidden-import",
        "statsmodels.stats.multitest",
        "--hidden-import",
        "matplotlib",
        "--hidden-import",
        "matplotlib.backends.backend_agg",
        "--hidden-import",
        "matplotlib.backends.backend_pdf",
        "--hidden-import",
        "matplotlib.backends.backend_svg",
        "--hidden-import",
        "matplotlib.backends.backend_qtagg",
        "--hidden-import",
        "networkx",
        "--hidden-import",
        "seaborn",
        "--hidden-import",
        "sklearn",
        "--hidden-import",
        "sklearn.decomposition",
        "--hidden-import",
        "sklearn.cross_decomposition",
        "--hidden-import",
        "sklearn.model_selection",
        "--hidden-import",
        "sklearn.metrics",
        "--hidden-import",
        "sklearn.impute",
        "--hidden-import",
        "scipy.cluster.hierarchy",
        "--hidden-import",
        "statsmodels.stats.anova",
        "--hidden-import",
        "patsy",
    ]

    subprocess.run(cmd, cwd=PROJECT_ROOT, check=True)

    built_dir = DIST_DIR / APP_NAME
    if not built_dir.exists():
        raise FileNotFoundError(f"打包结果不存在：{built_dir}")

    print("===== 整理 release 文件夹 =====")
    RELEASE_DIR.parent.mkdir(exist_ok=True)
    shutil.copytree(built_dir, RELEASE_DIR)

    print("===== 复制额外说明文件 =====")
    for name in ["README.md", "requirements.txt"]:
        src = PROJECT_ROOT / name
        if src.exists():
            shutil.copy2(src, RELEASE_DIR / name)

    print("===== 创建输出目录 =====")
    for name in ["output", "logs"]:
        (RELEASE_DIR / name).mkdir(exist_ok=True)

    print("")
    print("[OK] 发布包已生成：")
    print(RELEASE_DIR)
    print("")
    print("双击运行：")
    print(RELEASE_DIR / f"{APP_NAME}.exe")


if __name__ == "__main__":
    main()
