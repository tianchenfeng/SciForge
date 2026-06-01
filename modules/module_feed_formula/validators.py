from __future__ import annotations

from typing import List, Dict, Any


class FeedFormulaValidator:
    """
    输入数据校验模块
    负责在进入优化器前拦截所有脏数据
    """

    REQUIRED_FIELDS = ["name", "type", "cost"]

    NUMERIC_FIELDS = [
        "cost", "dm", "cp", "starch", "ndf", "adf", "ca", "p"
    ]

    @staticmethod
    def validate(
        ingredients: List[Dict[str, Any]],
        constraints: Dict[str, Dict[str, Any]],
    ) -> None:
        """
        主入口：统一调用
        出错直接 raise ValueError
        """
        FeedFormulaValidator._validate_ingredients(ingredients)
        FeedFormulaValidator._validate_constraints(constraints)

    # =========================
    # 原料检查
    # =========================
    @staticmethod
    def _validate_ingredients(ingredients: List[Dict[str, Any]]) -> None:
        if not ingredients:
            raise ValueError("原料列表为空，请至少添加一个原料。")

        names = []

        for i, row in enumerate(ingredients, start=1):
            name = str(row.get("name", "")).strip()

            if not name:
                raise ValueError(f"第 {i} 行原料名称为空。")

            names.append(name)

            # type 检查
            type_val = str(row.get("type", "")).strip().lower()
            if type_val not in ("roughage", "concentrate"):
                raise ValueError(f"{name} 的 type 必须是 roughage 或 concentrate。")

            # 数值字段检查
            for field in FeedFormulaValidator.NUMERIC_FIELDS:
                value = row.get(field, 0)

                try:
                    val = float(value)
                except (TypeError, ValueError):
                    raise ValueError(f"{name} 的 {field} 不是合法数字。")

                if val < 0:
                    raise ValueError(f"{name} 的 {field} 不能为负数。")

        # 名称重复检查
        if len(set(names)) != len(names):
            raise ValueError("原料名称存在重复，请保持唯一。")

    # =========================
    # 约束检查
    # =========================
    @staticmethod
    def _validate_constraints(constraints: Dict[str, Dict[str, Any]]) -> None:
        if not constraints:
            raise ValueError("约束设置为空。")

        enabled_count = 0

        for key, setting in constraints.items():
            enable = bool(setting.get("enable", False))

            if enable:
                enabled_count += 1

                target = setting.get("target", 0)
                error = setting.get("error", 0)

                try:
                    target_val = float(target)
                    error_val = float(error)
                except (TypeError, ValueError):
                    raise ValueError(f"{key} 的 target 或 error 不是合法数字。")

                if target_val < 0:
                    raise ValueError(f"{key} 的 target 不能为负数。")

                if error_val < 0:
                    raise ValueError(f"{key} 的 error 不能为负数。")

        if enabled_count == 0:
            raise ValueError("请至少启用一个约束条件。")

    # =========================
    # 粗料约束专项检查
    # =========================
    @staticmethod
    def validate_roughage_presence(
        ingredients: List[Dict[str, Any]],
        constraints: Dict[str, Dict[str, Any]],
    ) -> None:
        ratio_setting = constraints.get("ratio", {})

        if not ratio_setting.get("enable", False):
            return

        has_roughage = any(
            str(row.get("type", "")).strip().lower() == "roughage"
            for row in ingredients
        )

        if not has_roughage:
            raise ValueError("启用了粗料占比约束，但没有任何粗饲料原料。")