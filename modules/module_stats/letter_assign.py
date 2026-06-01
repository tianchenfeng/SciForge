import string
import numpy as np


def assign_letters_by_tukey(groups, mean_dict, p_matrix, alpha=0.05):
    """
    根据 Tukey 两两比较结果分配显著性字母。

    规则：
    - 均值越大，优先获得 a
    - 两组差异不显著，可以共享字母
    - 两组差异显著，不能共享字母
    """

    groups = list(groups)

    if len(groups) == 0:
        return {}

    if len(groups) == 1:
        return {groups[0]: ""}

    # 按均值从大到小排序
    groups_sorted = sorted(
        groups,
        key=lambda g: mean_dict.get(g, -np.inf),
        reverse=True
    )

    letters = list(string.ascii_lowercase)
    result = {g: "" for g in groups_sorted}
    letter_groups = []

    for g in groups_sorted:
        placed = False

        for i, existing_groups in enumerate(letter_groups):
            can_share = True

            for eg in existing_groups:
                p = p_matrix.get((g, eg), p_matrix.get((eg, g), np.nan))

                # 如果 p 缺失，为稳妥起见，不共享字母
                if np.isnan(p) or p < alpha:
                    can_share = False
                    break

            if can_share:
                result[g] += letters[i]
                existing_groups.append(g)
                placed = True
                break

        if not placed:
            idx = len(letter_groups)
            if idx >= len(letters):
                letter = f"x{idx + 1}"
            else:
                letter = letters[idx]

            result[g] += letter
            letter_groups.append([g])

    return result