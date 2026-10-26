from PySide6.QtCore import Qt, QSize
from PySide6.QtWidgets import QFrame, QGridLayout, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget


class HomePage(QWidget):
    MODULES = [
        ("stats", "统计分析", "ANOVA · 多重比较 · 批量报告", "Aa", "#007aff", "#eaf3ff"),
        ("16s", "16S 微生物组", "多样性 · 物种组成 · 功能预测", "16", "#20a38b", "#e6f7f1"),
        ("metabolomics", "代谢组分析", "质量控制 · 多变量 · 差异分析", "M", "#9461d4", "#f1eafa"),
        ("rnaseq", "RNA-Seq 分析", "差异基因 · 富集分析 · 候选筛选", "R", "#4977c1", "#eaf0fb"),
        ("gas_fit", "体外产气拟合", "动力学模型 · 参数估计 · 曲线拟合", "↗", "#db8b35", "#fff3e6"),
        ("silage_eval", "青贮综合评价", "隶属函数 · 熵权赋权 · 综合排序", "S", "#509e65", "#eaf5eb"),
        ("feed", "智能饲料配方", "营养约束 · 原料配比 · 成本优化", "F", "#cc7184", "#fceef2"),
    ]

    def __init__(self):
        super().__init__()
        self._columns = 0
        for key, *_ in self.MODULES:
            setattr(self, f"on_open_{key}", None)
        root = QVBoxLayout(self)
        root.setContentsMargins(4, 8, 4, 8)
        root.setSpacing(22)
        hero = QFrame()
        hero.setObjectName("hero")
        hero_layout = QVBoxLayout(hero)
        hero_layout.setContentsMargins(30, 28, 30, 30)
        hero_layout.setSpacing(12)
        for text, name in [
            ("SCIFORGE  /  RESEARCH WORKSPACE", "eyebrow"),
            ("让数据，成为发现。", "heroTitle"),
            ("从实验数据到清晰结论，让专注回到科研本身。", "heroSubtitle"),
        ]:
            label = QLabel(text)
            label.setObjectName(name)
            label.setWordWrap(True)
            hero_layout.addWidget(label)
        root.addWidget(hero)
        heading = QHBoxLayout()
        title = QLabel("开始你的探索")
        title.setObjectName("sectionTitle")
        heading.addWidget(title)
        heading.addStretch()
        hint = QLabel("7 个分析模块 · 一个工作空间")
        hint.setObjectName("mutedLabel")
        heading.addWidget(hint)
        root.addLayout(heading)
        self.grid = QGridLayout()
        self.grid.setSpacing(14)
        self.cards = []
        for key, title, description, symbol, color, background in self.MODULES:
            card = QPushButton()
            card.setObjectName("moduleCard")
            card.setAccessibleName(title)
            card.setMinimumSize(240, 142)
            card.setCursor(Qt.PointingHandCursor)
            card.clicked.connect(lambda checked=False, k=key: self._open(k))
            layout = QVBoxLayout(card)
            layout.setContentsMargins(20, 18, 20, 18)
            layout.setSpacing(8)
            top = QHBoxLayout()
            badge = QLabel(symbol)
            badge.setFixedSize(32, 32)
            badge.setAlignment(Qt.AlignCenter)
            badge.setStyleSheet(f"color: {color}; background: {background}; border-radius: 9px; font-weight: 600; font-size: 16px;")
            top.addWidget(badge)
            top.addStretch()
            arrow = QLabel("↗")
            arrow.setObjectName("cardArrow")
            top.addWidget(arrow)
            layout.addLayout(top)
            for text, name in [(title, "cardTitle"), (description, "cardDescription")]:
                label = QLabel(text)
                label.setWordWrap(True)
                label.setObjectName(name)
                layout.addWidget(label)
            for child in card.findChildren(QWidget):
                child.setAttribute(Qt.WA_TransparentForMouseEvents)
            setattr(self, f"btn_{key}", card)
            self.cards.append(card)
        root.addLayout(self.grid)
        note = QLabel("本地处理数据  ·  结果可追溯  ·  报告随时导出")
        note.setObjectName("mutedLabel")
        note.setWordWrap(True)
        root.addWidget(note)
        root.addStretch()
        self._reflow(3)

    def _open(self, key):
        callback = getattr(self, f"on_open_{key}", None)
        if callable(callback):
            callback()

    def _reflow(self, columns):
        if columns == self._columns:
            return
        self._columns = columns
        for card in self.cards:
            self.grid.removeWidget(card)
        for index, card in enumerate(self.cards):
            self.grid.addWidget(card, index // columns, index % columns)
        for column in range(3):
            self.grid.setColumnStretch(column, 1 if column < columns else 0)

    def resizeEvent(self, event):
        self._reflow(max(1, min(3, (event.size().width() + 14) // 270)))
        super().resizeEvent(event)

    def minimumSizeHint(self):
        return QSize(260, self.layout().minimumSize().height())
