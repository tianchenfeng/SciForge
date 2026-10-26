from pathlib import Path
import json
import os
import shutil
import subprocess
import sys


PROJECT_ROOT = Path(__file__).resolve().parent
APP_NAME = "SciForge"
VERSION = "4.1-desktop"
PACKAGE_NAME = f"{APP_NAME}-v{VERSION}-windows-x64"

DIST_DIR = PROJECT_ROOT / "dist" / VERSION
BUILD_DIR = PROJECT_ROOT / "build" / VERSION
RELEASE_DIR = PROJECT_ROOT / "release" / PACKAGE_NAME / APP_NAME


def remove_path(path: Path):
    resolved = path.resolve()
    if not resolved.is_relative_to(PROJECT_ROOT) or resolved == PROJECT_ROOT:
        raise ValueError(f"Refusing to remove a path outside the project: {resolved}")
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

    # Ship tracked application resources only, excluding local research and tests.
    bundle_dir = BUILD_DIR / "bundle_data"
    resource_names = ["assets", "templates", "example_data", "resources", "modules", "ui", "common", "config"]
    tracked = subprocess.check_output(["git", "ls-files", "-z", "--", *resource_names], cwd=PROJECT_ROOT)
    paths = {Path(name.decode("utf-8")) for name in tracked.split(b"\0") if name}
    paths.update({Path("ui/theme.py"), Path("common/runtime_check.py")})
    for relative in sorted(paths):
        if any(part.startswith("test") or part == "__pycache__" for part in relative.parts):
            continue
        source = PROJECT_ROOT / relative
        if source.is_file():
            target = bundle_dir / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
    for name in resource_names:
        (bundle_dir / name).mkdir(parents=True, exist_ok=True)

    print("===== 开始 PyInstaller 打包 =====")

    cmd = [
        sys.executable,
        "-m",
        "PyInstaller",
        "--noconfirm",
        "--onedir",
        "--windowed",
        "--distpath", str(DIST_DIR),
        "--workpath", str(BUILD_DIR / "pyinstaller"),
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

    for index, arg in enumerate(cmd):
        if index and cmd[index - 1] == "--add-data":
            source, destination = arg.split(";", 1)
            cmd[index] = f"{bundle_dir / source};{destination}"
    # Conda DLLs live outside site-packages. Make them visible to dependency
    # discovery so the release works without Conda on the recipient's PATH.
    build_env = os.environ.copy()
    runtime_dirs = [Path(sys.prefix) / "Library" / "bin", Path(sys.prefix), Path(sys.prefix) / "Scripts"]
    build_env["PATH"] = os.pathsep.join(str(path) for path in runtime_dirs if path.exists()) + os.pathsep + build_env.get("PATH", "")
    build_env["PYTHONIOENCODING"] = "utf-8"
    subprocess.run(cmd, cwd=PROJECT_ROOT, env=build_env, check=True)

    built_dir = DIST_DIR / APP_NAME
    if not built_dir.exists():
        raise FileNotFoundError(f"打包结果不存在：{built_dir}")

    print("===== 整理 release 文件夹 =====")
    RELEASE_DIR.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(built_dir, RELEASE_DIR)

    print("===== 复制额外说明文件 =====")
    for name in ["README.md", "requirements.txt", "docs/RELEASE_4.1.md", "docs/VALIDATION_4.1.md"]:
        src = PROJECT_ROOT / name
        if src.exists():
            shutil.copy2(src, RELEASE_DIR / Path(name).name)

    print("===== 创建输出目录 =====")
    for name in ["output", "logs"]:
        (RELEASE_DIR / name).mkdir(exist_ok=True)

    files = [str(path.relative_to(RELEASE_DIR)) for path in RELEASE_DIR.rglob("*") if path.is_file()]
    (RELEASE_DIR / "package_manifest.json").write_text(json.dumps({"version":VERSION,"files":files},ensure_ascii=False,indent=2),encoding="utf-8")

    print("")
    print("[OK] 发布包已生成：")
    print(RELEASE_DIR)
    print("")
    print("双击运行：")
    print(RELEASE_DIR / f"{APP_NAME}.exe")


if __name__ == "__main__":
    main()
