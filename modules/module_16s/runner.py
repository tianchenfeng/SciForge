import os
import sys
import tempfile

from PySide6.QtCore import QProcess, QProcessEnvironment

from modules.module_16s.service import SixteenSService


class SixteenSRunner:
    def __init__(self, parent=None):
        self.parent = parent

    def _build_process_environment(self) -> QProcessEnvironment:
        """
        为子进程构建稳定的环境变量，重点修复：
        1. Path.home() 无法识别 HOME/USERPROFILE
        2. matplotlib 无法创建配置目录
        """
        env = QProcessEnvironment.systemEnvironment()

        # 先尽量拿到一个可靠的用户主目录
        home_dir = (
            os.environ.get("USERPROFILE")
            or os.environ.get("HOME")
            or (
                (os.environ.get("HOMEDRIVE", "") + os.environ.get("HOMEPATH", ""))
                if os.environ.get("HOMEDRIVE") and os.environ.get("HOMEPATH")
                else ""
            )
        )

        if not home_dir:
            home_dir = os.path.expanduser("~")

        if not home_dir or home_dir == "~":
            home_dir = tempfile.gettempdir()

        # 保证这些变量在子进程里存在
        env.insert("USERPROFILE", home_dir)
        env.insert("HOME", home_dir)

        if os.path.splitdrive(home_dir)[0]:
            env.insert("HOMEDRIVE", os.path.splitdrive(home_dir)[0])

        home_tail = home_dir[len(os.path.splitdrive(home_dir)[0]):]
        if home_tail:
            env.insert("HOMEPATH", home_tail)

        # 给 matplotlib 单独指定可写目录，避免它去猜 home 出问题
        mpl_config_dir = os.path.join(home_dir, ".matplotlib")
        os.makedirs(mpl_config_dir, exist_ok=True)
        env.insert("MPLCONFIGDIR", mpl_config_dir)

        # 编码稳一点
        env.insert("PYTHONIOENCODING", "utf-8")
        env.insert("PYTHONUTF8", "1")

        return env

    def create_process(
        self,
        script_name: str,
        input_dir: str,
        output_dir: str,
        mode: str,
        app_base_dir: str,
    ) -> QProcess:
        SixteenSService.validate_script_exists(script_name)

        process = QProcess(self.parent)
        process.setProgram(sys.executable)

        if getattr(sys, "frozen", False):
            args = SixteenSService.build_frozen_args(
                script_name=script_name,
                input_dir=input_dir,
                output_dir=output_dir,
                mode=mode,
                app_base_dir=app_base_dir,
            )
        else:
            script_path = SixteenSService.get_script_path(script_name)
            args = SixteenSService.build_script_args(
                script_path=script_path,
                input_dir=input_dir,
                output_dir=output_dir,
                mode=mode,
                app_base_dir=app_base_dir,
            )

        process.setArguments(args)
        process.setProcessEnvironment(self._build_process_environment())

        return process