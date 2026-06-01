from __future__ import annotations

from typing import Dict, List, Any

from modules.module_feed_formula.validators import FeedFormulaValidator
from modules.module_feed_formula.optimizer import FeedFormulaOptimizer
from modules.module_feed_formula.nutrient_calculator import NutrientCalculator


class FeedFormulaService:
    """
    饲料配方模块总服务层
    职责：
    1. 校验输入
    2. 调用优化器求解
    3. 回算营养结果
    4. 输出统一结果格式
    """

    @staticmethod
    def run_formula(
        ingredients: List[Dict[str, Any]],
        constraints: Dict[str, Dict[str, Any]],
    ) -> Dict[str, Any]:
        """
        返回统一结果结构：
        {
            "success": True/False,
            "message": "...",
            "mix": [...],
            "stats": {...},
            "solver_status": "..."
        }
        """
        # 1. 输入校验
        FeedFormulaValidator.validate(ingredients, constraints)
        FeedFormulaValidator.validate_roughage_presence(ingredients, constraints)

        # 2. 优化求解
        optimizer = FeedFormulaOptimizer()
        result = optimizer.solve(ingredients, constraints)

        # 3. 如果无解，直接返回
        if not result["success"]:
            return {
                "success": False,
                "message": result.get("message", "求解失败。"),
                "mix": [],
                "stats": {},
                "solver_status": result.get("solver_status", ""),
            }

        # 4. 回算营养
        stats = NutrientCalculator.calculate(result["mix"])

        # 5. 输出统一结果
        return {
            "success": True,
            "message": result.get("message", "求解成功。"),
            "mix": result["mix"],
            "stats": stats,
            "solver_status": result.get("solver_status", ""),
        }