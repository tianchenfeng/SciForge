from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent
TEST_DIR = BASE_DIR / "tests"
TEST_DIR.mkdir(exist_ok=True)

df = pd.DataFrame({
    "Sample": [
        "S1", "S2", "S3",
        "S4", "S5", "S6",
        "S7", "S8", "S9",
        "S10", "S11", "S12"
    ],
    "A": [
        "C", "C", "C",
        "L", "L", "L",
        "M", "M", "M",
        "E", "E", "E"
    ],
    "B": [
        "3", "7", "15",
        "3", "7", "15",
        "3", "7", "15",
        "3", "7", "15"
    ],
    "value": [
        10.2, 11.1, 12.3,
        13.2, 15.4, 16.8,
        14.1, 17.6, 18.2,
        12.8, 14.6, 15.9
    ]
})

out_path = TEST_DIR / "demo_data.xlsx"
df.to_excel(out_path, index=False)

print(f"[OK] 示例数据已生成: {out_path}")