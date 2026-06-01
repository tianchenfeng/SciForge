from pathlib import Path
from typing import List

from modules.module_16s.config import SCRIPTS_DIR


class SixteenSService:
    @staticmethod
    def get_script_path(script_name: str) -> Path:
        return SCRIPTS_DIR / script_name

    @staticmethod
    def validate_script_exists(script_name: str) -> None:
        script_path = SixteenSService.get_script_path(script_name)
        if not script_path.exists():
            raise FileNotFoundError(f"找不到脚本文件：{script_path}")

    @staticmethod
    def normalize_selected_scripts(selected_scripts: List[str], mode: str) -> List[str]:
        scripts = list(selected_scripts)

        if mode == "single_factor" and "09_diff_visualization.py" in scripts:
            scripts.remove("09_diff_visualization.py")

        return scripts

    @staticmethod
    def build_script_args(
        script_path: Path,
        input_dir: str,
        output_dir: str,
        mode: str,
        app_base_dir: str,
    ) -> list[str]:
        args = [
            str(script_path),
            "--input_dir", input_dir,
            "--output_dir", output_dir,
            "--mode", mode,
        ]

        # 只有 05_picrust2.py 才额外传 app_base_dir
        if script_path.name == "05_picrust2.py":
            args.extend(["--app_base_dir", app_base_dir])

        return args

    @staticmethod
    def build_frozen_args(
        script_name: str,
        input_dir: str,
        output_dir: str,
        mode: str,
        app_base_dir: str,
    ) -> list[str]:
        args = [
            "--run-script",
            str(Path("modules") / "module_16s" / "scripts" / script_name),
            "--input_dir", input_dir,
            "--output_dir", output_dir,
            "--mode", mode,
        ]

        # 只有 05_picrust2.py 才额外传 app_base_dir
        if script_name == "05_picrust2.py":
            args.extend(["--app_base_dir", app_base_dir])

        return args