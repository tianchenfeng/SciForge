from pathlib import Path
import pandas as pd


SUPPORTED_EXTENSIONS = [
    ".xlsx",
    ".csv",
    ".tsv",
    ".txt"
]


def load_table_auto(file_path, **kwargs):
    """
    自动识别表格格式并读取

    支持：
    - xlsx
    - csv
    - tsv
    - txt

    Parameters
    ----------
    file_path : str or Path

    Returns
    -------
    pandas.DataFrame
    """

    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"文件不存在: {file_path}")

    suffix = file_path.suffix.lower()

    # Excel
    if suffix in [".xlsx"]:
        return pd.read_excel(file_path, **kwargs)

    # CSV
    elif suffix == ".csv":
        return pd.read_csv(file_path, **kwargs)

    # TSV/TXT
    elif suffix in [".tsv", ".txt"]:
        return pd.read_csv(file_path, sep="\t", **kwargs)

    else:
        raise ValueError(
            f"不支持的文件格式: {suffix}\n"
            f"支持格式: {SUPPORTED_EXTENSIONS}"
        )


def find_first_table(folder):
    """
    在文件夹中自动寻找第一个可识别表格
    """

    folder = Path(folder)

    if not folder.exists():
        raise FileNotFoundError(f"文件夹不存在: {folder}")

    for ext in SUPPORTED_EXTENSIONS:
        files = list(folder.glob(f"*{ext}"))

        if files:
            return files[0]

    raise FileNotFoundError(
        f"{folder} 中未找到支持的表格文件"
    )