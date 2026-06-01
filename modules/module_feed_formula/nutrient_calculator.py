from __future__ import annotations

from typing import Dict, List, Any


class NutrientCalculator:
    """
    根据优化器输出的 mix 结果，回算配方整体营养值和成本。
    mix 中每一项应包含：
    - amount
    - cost
    - type
    - cp / starch / ndf / adf / ca / p
    """

    @staticmethod
    def calculate(mix: List[Dict[str, Any]]) -> Dict[str, float]:
        stats = {
            "cost": 0.0,
            "roughage": 0.0,
            "cp": 0.0,
            "starch": 0.0,
            "ndf": 0.0,
            "adf": 0.0,
            "ca": 0.0,
            "p": 0.0,
        }

        for row in mix:
            amount = float(row.get("amount", 0.0))

            stats["cost"] += amount * float(row.get("cost", 0.0))

            if str(row.get("type", "")).strip().lower() == "roughage":
                stats["roughage"] += amount

            stats["cp"] += amount * float(row.get("cp", 0.0)) / 100.0
            stats["starch"] += amount * float(row.get("starch", 0.0)) / 100.0
            stats["ndf"] += amount * float(row.get("ndf", 0.0)) / 100.0
            stats["adf"] += amount * float(row.get("adf", 0.0)) / 100.0
            stats["ca"] += amount * float(row.get("ca", 0.0)) / 100.0
            stats["p"] += amount * float(row.get("p", 0.0)) / 100.0

        return {k: round(v, 6) for k, v in stats.items()}