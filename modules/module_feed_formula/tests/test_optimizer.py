from modules.module_feed_formula.optimizer import FeedFormulaOptimizer


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

    optimizer = FeedFormulaOptimizer()
    result = optimizer.solve(ingredients, constraints)

    print(result)


if __name__ == "__main__":
    main()