from __future__ import annotations

from typing import Dict, List, Any

import numpy as np
from scipy.optimize import linprog


class FeedFormulaOptimizer:
    """
    自由配方版线性规划求解器（SciPy 版）

    目标：
    - 最小化总成本
    - 总配方量固定为 100
    - 支持营养指标 target ± error 约束
    - 支持粗饲料占比约束（ratio -> roughage）
    """

    NUTRIENT_KEYS = ("cp", "starch", "ndf", "adf", "ca", "p")

    def solve(
        self,
        ingredients: List[Dict[str, Any]],
        constraints: Dict[str, Dict[str, Any]],
    ) -> Dict[str, Any]:
        clean_ingredients = self._normalize_ingredients(ingredients)
        clean_constraints = self._normalize_constraints(constraints)

        n = len(clean_ingredients)
        if n == 0:
            raise ValueError("原料列表为空。")

        # 目标函数：成本最小
        c = np.array([ing["cost"] for ing in clean_ingredients], dtype=float)

        # 等式约束：总量 = 100
        A_eq = [np.ones(n, dtype=float)]
        b_eq = [100.0]

        # 不等式约束：A_ub x <= b_ub
        A_ub = []
        b_ub = []

        # 营养约束
        for key in self.NUTRIENT_KEYS:
            setting = clean_constraints.get(key)
            if setting and setting["enable"]:
                coeff = np.array([ing[key] / 100.0 for ing in clean_ingredients], dtype=float)
                min_val = max(0.0, setting["target"] - setting["error"])
                max_val = setting["target"] + setting["error"]

                # coeff * x <= max_val
                A_ub.append(coeff)
                b_ub.append(max_val)

                # coeff * x >= min_val  ->  -coeff * x <= -min_val
                A_ub.append(-coeff)
                b_ub.append(-min_val)

        # 粗饲料占比约束
        ratio_setting = clean_constraints.get("ratio")
        if ratio_setting and ratio_setting["enable"]:
            roughage_coeff = np.array(
                [1.0 if ing["type"] == "roughage" else 0.0 for ing in clean_ingredients],
                dtype=float,
            )
            min_ratio = max(0.0, ratio_setting["target"] - ratio_setting["error"])
            max_ratio = min(100.0, ratio_setting["target"] + ratio_setting["error"])

            # roughage <= max_ratio
            A_ub.append(roughage_coeff)
            b_ub.append(max_ratio)

            # roughage >= min_ratio -> -roughage <= -min_ratio
            A_ub.append(-roughage_coeff)
            b_ub.append(-min_ratio)

        A_eq = np.array(A_eq, dtype=float)
        b_eq = np.array(b_eq, dtype=float)

        if A_ub:
            A_ub = np.array(A_ub, dtype=float)
            b_ub = np.array(b_ub, dtype=float)
        else:
            A_ub = None
            b_ub = None

        # 每种原料 0~100
        bounds = [(0.0, 100.0) for _ in range(n)]

        result = linprog(
            c=c,
            A_ub=A_ub,
            b_ub=b_ub,
            A_eq=A_eq,
            b_eq=b_eq,
            bounds=bounds,
            method="highs",
        )

        if not result.success:
            return {
                "success": False,
                "message": self._build_failure_message(result),
                "mix": [],
                "solver_status": result.message,
            }

        final_mix: List[Dict[str, Any]] = []
        x = result.x

        for ing, amount in zip(clean_ingredients, x):
            amount = float(amount)
            if amount > 1e-6:
                row = dict(ing)
                row["amount"] = round(amount, 6)
                final_mix.append(row)

        final_mix.sort(key=lambda item: item["amount"], reverse=True)

        return {
            "success": True,
            "message": "求解成功。",
            "mix": final_mix,
            "solver_status": result.message,
        }

    def _normalize_ingredients(self, ingredients: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        clean_rows: List[Dict[str, Any]] = []

        for i, ing in enumerate(ingredients, start=1):
            name = str(ing.get("name", "")).strip()
            if not name:
                raise ValueError(f"第 {i} 行原料名称为空。")

            row = {
                "name": name,
                "type": str(ing.get("type", "concentrate")).strip().lower(),
                "cost": self._to_non_negative_float(ing.get("cost", 0), f"{name} 的 cost"),
                "dm": self._to_non_negative_float(ing.get("dm", 0), f"{name} 的 dm"),
                "cp": self._to_non_negative_float(ing.get("cp", 0), f"{name} 的 cp"),
                "starch": self._to_non_negative_float(ing.get("starch", 0), f"{name} 的 starch"),
                "ndf": self._to_non_negative_float(ing.get("ndf", 0), f"{name} 的 ndf"),
                "adf": self._to_non_negative_float(ing.get("adf", 0), f"{name} 的 adf"),
                "ca": self._to_non_negative_float(ing.get("ca", 0), f"{name} 的 ca"),
                "p": self._to_non_negative_float(ing.get("p", 0), f"{name} 的 p"),
            }

            if row["type"] not in ("roughage", "concentrate"):
                raise ValueError(f"{name} 的 type 只能是 roughage 或 concentrate。")

            clean_rows.append(row)

        names = [x["name"] for x in clean_rows]
        if len(set(names)) != len(names):
            raise ValueError("原料名称存在重复，请保持唯一。")

        return clean_rows

    def _normalize_constraints(
        self,
        constraints: Dict[str, Dict[str, Any]],
    ) -> Dict[str, Dict[str, Any]]:
        clean: Dict[str, Dict[str, Any]] = {}

        for key, value in constraints.items():
            clean[key] = {
                "enable": bool(value.get("enable", False)),
                "target": self._to_non_negative_float(value.get("target", 0), f"{key}.target"),
                "error": self._to_non_negative_float(value.get("error", 0), f"{key}.error"),
            }

        return clean

    @staticmethod
    def _to_non_negative_float(value: Any, field_name: str) -> float:
        try:
            val = float(value)
        except (TypeError, ValueError):
            raise ValueError(f"{field_name} 不是合法数字。")

        if val < 0:
            raise ValueError(f"{field_name} 不能为负数。")

        return val

    @staticmethod
    def _build_failure_message(result) -> str:
        msg = str(result.message)
        if "infeasible" in msg.lower():
            return "无解。请放宽误差范围，或检查原料数据与约束条件是否冲突。"
        return f"求解失败：{msg}"