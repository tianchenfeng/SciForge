from pathlib import Path
import ast
from datetime import datetime


PROJECT_ROOT = Path(r"F:\WorkSpace\toolbox")
OUTPUT_FILE = PROJECT_ROOT / "project_snapshot.md"

IGNORE_DIRS = {
    "__pycache__",
    ".git",
    ".idea",
    ".venv",
    "venv",
    "env",
    "build",
    "dist",
    "output",
    "outputs",
    "logs",
    ".pytest_cache",
}

IGNORE_SUFFIXES = {
    ".pyc",
    ".pyo",
    ".exe",
    ".dll",
    ".so",
    ".pyd",
    ".zip",
    ".rar",
    ".7z",
    ".png",
    ".jpg",
    ".jpeg",
    ".gif",
    ".pdf",
    ".docx",
    ".xlsx",
    ".xls",
    ".csv",
    ".tsv",
}


def should_ignore(path: Path) -> bool:
    parts = set(path.parts)
    if parts & IGNORE_DIRS:
        return True

    if path.suffix.lower() in IGNORE_SUFFIXES:
        return True

    return False


def read_text_safe(path: Path) -> str:
    for enc in ("utf-8", "utf-8-sig", "gbk"):
        try:
            return path.read_text(encoding=enc)
        except UnicodeDecodeError:
            continue
        except Exception:
            return ""

    return ""


def get_tree(root: Path) -> list[str]:
    lines = []

    def walk(current: Path, prefix: str = ""):
        children = [
            p for p in sorted(current.iterdir(), key=lambda x: (x.is_file(), x.name.lower()))
            if not should_ignore(p)
        ]

        for i, child in enumerate(children):
            connector = "└── " if i == len(children) - 1 else "├── "
            lines.append(f"{prefix}{connector}{child.name}")

            if child.is_dir():
                extension = "    " if i == len(children) - 1 else "│   "
                walk(child, prefix + extension)

    lines.append(root.name)
    walk(root)
    return lines


def analyze_py_file(path: Path) -> dict:
    text = read_text_safe(path)
    result = {
        "path": str(path.relative_to(PROJECT_ROOT)),
        "lines": len(text.splitlines()),
        "imports": [],
        "classes": [],
        "functions": [],
        "parse_error": "",
    }

    if not text.strip():
        return result

    try:
        tree = ast.parse(text)
    except SyntaxError as e:
        result["parse_error"] = f"SyntaxError: line {e.lineno}, {e.msg}"
        return result
    except Exception as e:
        result["parse_error"] = repr(e)
        return result

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                result["imports"].append(alias.name)

        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            names = ", ".join(alias.name for alias in node.names)
            result["imports"].append(f"from {module} import {names}")

        elif isinstance(node, ast.ClassDef):
            result["classes"].append(node.name)

        elif isinstance(node, ast.FunctionDef):
            result["functions"].append(node.name)

    result["imports"] = sorted(set(result["imports"]))
    result["classes"] = sorted(set(result["classes"]))
    result["functions"] = sorted(set(result["functions"]))

    return result


def main():
    py_files = [
        p for p in PROJECT_ROOT.rglob("*.py")
        if not should_ignore(p)
    ]

    py_infos = [analyze_py_file(p) for p in sorted(py_files)]

    lines = []
    lines.append("# 田晨丰工具箱项目快照")
    lines.append("")
    lines.append(f"- 扫描时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"- 项目路径：`{PROJECT_ROOT}`")
    lines.append(f"- Python 文件数量：{len(py_infos)}")
    lines.append("")

    lines.append("## 一、目录结构")
    lines.append("")
    lines.append("```text")
    lines.extend(get_tree(PROJECT_ROOT))
    lines.append("```")
    lines.append("")

    lines.append("## 二、Python 文件总览")
    lines.append("")
    lines.append("| 文件 | 行数 | 类 | 函数 | 解析状态 |")
    lines.append("|---|---:|---|---|---|")

    for info in py_infos:
        cls = ", ".join(info["classes"]) if info["classes"] else "-"
        funcs = ", ".join(info["functions"][:12]) if info["functions"] else "-"
        if len(info["functions"]) > 12:
            funcs += f" ... +{len(info['functions']) - 12}"
        status = info["parse_error"] if info["parse_error"] else "OK"

        lines.append(
            f"| `{info['path']}` | {info['lines']} | {cls} | {funcs} | {status} |"
        )

    lines.append("")

    lines.append("## 三、逐文件详情")
    lines.append("")

    for info in py_infos:
        lines.append(f"### `{info['path']}`")
        lines.append("")
        lines.append(f"- 行数：{info['lines']}")

        if info["parse_error"]:
            lines.append(f"- 解析错误：{info['parse_error']}")

        lines.append("- 类：")
        if info["classes"]:
            for item in info["classes"]:
                lines.append(f"  - `{item}`")
        else:
            lines.append("  - 无")

        lines.append("- 函数：")
        if info["functions"]:
            for item in info["functions"]:
                lines.append(f"  - `{item}`")
        else:
            lines.append("  - 无")

        lines.append("- Imports：")
        if info["imports"]:
            for item in info["imports"]:
                lines.append(f"  - `{item}`")
        else:
            lines.append("  - 无")

        lines.append("")

    OUTPUT_FILE.write_text("\n".join(lines), encoding="utf-8")
    print(f"[OK] 项目快照已生成：{OUTPUT_FILE}")


if __name__ == "__main__":
    main()