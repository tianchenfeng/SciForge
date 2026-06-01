from modules.module_feed_formula.service import FeedFormulaService


def main():
    ingredients = [
        {
            "name": "原料A",
            "type": "roughage",
            "cost": 2.5,
            "dm": 90,
            "cp": 18,
            "starch": 1.2,
            "ndf": 42,
            "adf": 30,
            "ca": 1.4,
            "p": 0.25,
        },
        {
            "name": "原料B",
            "type": "concentrate",
            "cost": 3.2,
            "dm": 89,
            "cp": 46.5,
            "starch": 5.0,
            "ndf": 12,
            "adf": 5,
            "ca": 0.35,
            "p": 0.65,
        },
        {
            "name": "原料C",
            "type": "concentrate",
            "cost": 2.0,
            "dm": 88,
            "cp": 8.5,
            "starch": 70,
            "ndf": 10,
            "adf": 3,
            "ca": 0.03,
            "p": 0.28,
        },
    ]

    constraints = {
        "cp": {"enable": True, "target": 17.5, "error": 1.0},
        "starch": {"enable": True, "target": 24.0, "error": 2.0},
        "ndf": {"enable": False, "target": 30.0, "error": 2.0},
        "adf": {"enable": False, "target": 19.0, "error": 2.0},
        "ca": {"enable": False, "target": 0.6, "error": 0.1},
        "p": {"enable": False, "target": 0.4, "error": 0.05},
        "ratio": {"enable": True, "target": 50.0, "error": 5.0},
    }

    result = FeedFormulaService.run_formula(ingredients, constraints)

    print("\n===== 配方结果 =====\n")

    if not result["success"]:
        print("失败：", result["message"])
        return

    for row in result["mix"]:
        print(f"{row['name']:10s}  {row['amount']:6.2f} %")

    stats = result["stats"]

    print("\n===== 营养验算 =====\n")
    print(f"总成本      : {stats['cost']:.4f}")
    print(f"粗料占比    : {stats['roughage']:.4f} %")
    print(f"粗蛋白 CP   : {stats['cp']:.4f} %")
    print(f"淀粉        : {stats['starch']:.4f} %")
    print(f"NDF         : {stats['ndf']:.4f} %")
    print(f"ADF         : {stats['adf']:.4f} %")
    print(f"Ca          : {stats['ca']:.4f} %")
    print(f"P           : {stats['p']:.4f} %")

    print("\n状态：", result["solver_status"])


if __name__ == "__main__":
    main()