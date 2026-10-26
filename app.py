import os
import sys
import runpy
import argparse


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
    ]
    if os.path.basename(args.run_script) == "05_picrust2.py":
        sys.argv.extend(["--app_base_dir", base_dir])

    runpy.run_path(script_path, run_name="__main__")
    return True


def run_metabolomics() -> int | None:
    """Run the process-isolated metabolomics pipeline without starting Qt."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--run-metabolomics", action="store_true")
    parser.add_argument("--config", default=None)
    args, _ = parser.parse_known_args()
    if not args.run_metabolomics:
        return None
    if not args.config:
        raise ValueError("--config is required with --run-metabolomics")
    from modules.module_metabolomics.runner import main as metabolomics_main
    return metabolomics_main(["--config", args.config])


def main() -> int:
    if getattr(sys, "frozen", False) and any(flag in sys.argv for flag in ("--run-script", "--run-metabolomics", "--self-test")):
        _restore_child_streams()
    if "--self-test" in sys.argv:
        from common.runtime_check import main as runtime_check
        return runtime_check(sys.argv[1:])
    metabolomics_exit = run_metabolomics()
    if metabolomics_exit is not None:
        return metabolomics_exit
    if run_embedded_script():
        return 0

    from PySide6.QtWidgets import QApplication
    from ui.main_window import MainWindow
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    return app.exec()


def _restore_child_streams():
    """Windowed PyInstaller children still need their inherited QProcess pipes."""
    if os.name != "nt":
        return
    import ctypes
    from ctypes import wintypes
    import msvcrt
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.GetStdHandle.argtypes = [wintypes.DWORD]
    kernel.GetStdHandle.restype = wintypes.HANDLE
    for name, handle_id in [("stdout", -11), ("stderr", -12)]:
        if getattr(sys, name) is not None:
            continue
        handle = kernel.GetStdHandle(handle_id & 0xffffffff)
        if handle and handle != wintypes.HANDLE(-1).value:
            try:
                fd = msvcrt.open_osfhandle(handle, os.O_WRONLY | os.O_BINARY)
                setattr(sys, name, os.fdopen(fd, "w", encoding="utf-8", buffering=1))
            except OSError:
                pass


if __name__ == "__main__":
    sys.exit(main())
