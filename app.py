import os
import sys
import runpy
import argparse

from PySide6.QtWidgets import QApplication
from ui.main_window import MainWindow


def get_base_dir() -> str:
    """
    源码运行时：返回项目根目录
    exe 运行时：返回 PyInstaller 解包目录
    """
    if getattr(sys, "frozen", False):
        return sys._MEIPASS
    return os.path.dirname(os.path.abspath(__file__))


def run_embedded_script() -> bool:
    """
    打包后子进程运行 16S 脚本时使用：
    app.exe --run-script modules/module_16s/scripts/01_qc.py --input_dir ... --output_dir ... --mode ...
    不启动 GUI，直接执行对应脚本。
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--run-script", dest="run_script", default=None)
    parser.add_argument("--input_dir", default=None)
    parser.add_argument("--output_dir", default=None)
    parser.add_argument("--mode", default=None)
    parser.add_argument("--app_base_dir", default=None)
    args, _ = parser.parse_known_args()

    if not args.run_script:
        return False

    base_dir = args.app_base_dir if args.app_base_dir else get_base_dir()
    script_path = os.path.join(base_dir, args.run_script)

    if not os.path.exists(script_path):
        raise FileNotFoundError(f"找不到脚本文件：{script_path}")

    sys.argv = [
        args.run_script,
        "--input_dir", args.input_dir or "",
        "--output_dir", args.output_dir or "",
        "--mode", args.mode or "two_factor",
        "--app_base_dir", base_dir,
    ]

    runpy.run_path(script_path, run_name="__main__")
    return True


def main() -> int:
    if run_embedded_script():
        return 0

    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())