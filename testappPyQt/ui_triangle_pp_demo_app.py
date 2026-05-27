# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'TrianglePPDemoApp.ui'
##
## Created by: Qt User Interface Compiler version 6.11.1
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QBrush, QColor, QConicalGradient, QCursor,
    QFont, QFontDatabase, QGradient, QIcon,
    QImage, QKeySequence, QLinearGradient, QPainter,
    QPalette, QPixmap, QRadialGradient, QTransform)
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QHBoxLayout,
    QLabel, QMainWindow, QMenuBar, QPushButton,
    QScrollArea, QSizePolicy, QSpacerItem, QStatusBar,
    QToolBar, QToolButton, QVBoxLayout, QWidget)

from drawing_area import DrawingArea
import TrianglePPDemoApp_rc

class Ui_TrianglePPDemoAppClass(object):
    def setupUi(self, TrianglePPDemoAppClass):
        if not TrianglePPDemoAppClass.objectName():
            TrianglePPDemoAppClass.setObjectName(u"TrianglePPDemoAppClass")
        TrianglePPDemoAppClass.resize(857, 576)
        icon = QIcon()
        icon.addFile(u":/TrianglePPDemo/triangle-PP-ico.png", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        TrianglePPDemoAppClass.setWindowIcon(icon)
        self.centralWidget = QWidget(TrianglePPDemoAppClass)
        self.centralWidget.setObjectName(u"centralWidget")
        self.verticalLayout = QVBoxLayout(self.centralWidget)
        self.verticalLayout.setSpacing(6)
        self.verticalLayout.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.horizontalLayout = QHBoxLayout()
        self.horizontalLayout.setSpacing(6)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.horizontalLayout.setContentsMargins(9, -1, -1, -1)
        self.label = QLabel(self.centralWidget)
        self.label.setObjectName(u"label")

        self.horizontalLayout.addWidget(self.label)

        self.pointModeComboBox = QComboBox(self.centralWidget)
        self.pointModeComboBox.addItem("")
        self.pointModeComboBox.addItem("")
        self.pointModeComboBox.addItem("")
        self.pointModeComboBox.addItem("")
        self.pointModeComboBox.addItem("")
        self.pointModeComboBox.addItem("")
        self.pointModeComboBox.setObjectName(u"pointModeComboBox")

        self.horizontalLayout.addWidget(self.pointModeComboBox)

        self.generatePointsPushButton = QPushButton(self.centralWidget)
        self.generatePointsPushButton.setObjectName(u"generatePointsPushButton")

        self.horizontalLayout.addWidget(self.generatePointsPushButton)

        self.horizontalSpacer = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout.addItem(self.horizontalSpacer)

        self.triangualtePointsPushButton = QPushButton(self.centralWidget)
        self.triangualtePointsPushButton.setObjectName(u"triangualtePointsPushButton")

        self.horizontalLayout.addWidget(self.triangualtePointsPushButton)

        self.verticalLayout_2 = QVBoxLayout()
        self.verticalLayout_2.setSpacing(2)
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.useConstraintsCheckBox = QCheckBox(self.centralWidget)
        self.useConstraintsCheckBox.setObjectName(u"useConstraintsCheckBox")

        self.verticalLayout_2.addWidget(self.useConstraintsCheckBox)

        self.hideMarkersCheckBox = QCheckBox(self.centralWidget)
        self.hideMarkersCheckBox.setObjectName(u"hideMarkersCheckBox")

        self.verticalLayout_2.addWidget(self.hideMarkersCheckBox)


        self.horizontalLayout.addLayout(self.verticalLayout_2)

        self.horizontalSpacer_3 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout.addItem(self.horizontalSpacer_3)

        self.tesselatePointsPushButton = QPushButton(self.centralWidget)
        self.tesselatePointsPushButton.setObjectName(u"tesselatePointsPushButton")

        self.horizontalLayout.addWidget(self.tesselatePointsPushButton)

        self.horizontalSpacer_2 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout.addItem(self.horizontalSpacer_2)

        self.optionsToolButton = QToolButton(self.centralWidget)
        self.optionsToolButton.setObjectName(u"optionsToolButton")
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.optionsToolButton.sizePolicy().hasHeightForWidth())
        self.optionsToolButton.setSizePolicy(sizePolicy)
        self.optionsToolButton.setAutoRaise(True)

        self.horizontalLayout.addWidget(self.optionsToolButton)


        self.verticalLayout.addLayout(self.horizontalLayout)

        self.scrollArea = QScrollArea(self.centralWidget)
        self.scrollArea.setObjectName(u"scrollArea")
        self.scrollArea.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.scrollArea.setWidgetResizable(True)
        self.drawAreaWidget = DrawingArea()
        self.drawAreaWidget.setObjectName(u"drawAreaWidget")
        self.drawAreaWidget.setGeometry(QRect(0, 0, 837, 457))
        sizePolicy1 = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Preferred)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.drawAreaWidget.sizePolicy().hasHeightForWidth())
        self.drawAreaWidget.setSizePolicy(sizePolicy1)
        self.scrollArea.setWidget(self.drawAreaWidget)

        self.verticalLayout.addWidget(self.scrollArea)

        TrianglePPDemoAppClass.setCentralWidget(self.centralWidget)
        self.menuBar = QMenuBar(TrianglePPDemoAppClass)
        self.menuBar.setObjectName(u"menuBar")
        self.menuBar.setGeometry(QRect(0, 0, 857, 21))
        TrianglePPDemoAppClass.setMenuBar(self.menuBar)
        self.mainToolBar = QToolBar(TrianglePPDemoAppClass)
        self.mainToolBar.setObjectName(u"mainToolBar")
        TrianglePPDemoAppClass.addToolBar(Qt.ToolBarArea.TopToolBarArea, self.mainToolBar)
        self.statusBar = QStatusBar(TrianglePPDemoAppClass)
        self.statusBar.setObjectName(u"statusBar")
        TrianglePPDemoAppClass.setStatusBar(self.statusBar)

        self.retranslateUi(TrianglePPDemoAppClass)

        QMetaObject.connectSlotsByName(TrianglePPDemoAppClass)
    # setupUi

    def retranslateUi(self, TrianglePPDemoAppClass):
        TrianglePPDemoAppClass.setWindowTitle(QCoreApplication.translate("TrianglePPDemoAppClass", u"Triangle++ Demo", None))
        self.label.setText(QCoreApplication.translate("TrianglePPDemoAppClass", u"Point Generation", None))
        self.pointModeComboBox.setItemText(0, QCoreApplication.translate("TrianglePPDemoAppClass", u"Manual", None))
        self.pointModeComboBox.setItemText(1, QCoreApplication.translate("TrianglePPDemoAppClass", u"Automatic", None))
        self.pointModeComboBox.setItemText(2, QCoreApplication.translate("TrianglePPDemoAppClass", u"From Image", None))
        self.pointModeComboBox.setItemText(3, QCoreApplication.translate("TrianglePPDemoAppClass", u"From File", None))
        self.pointModeComboBox.setItemText(4, QCoreApplication.translate("TrianglePPDemoAppClass", u"Example 1", None))
        self.pointModeComboBox.setItemText(5, QCoreApplication.translate("TrianglePPDemoAppClass", u"Example 2", None))

        self.generatePointsPushButton.setText(QCoreApplication.translate("TrianglePPDemoAppClass", u"Generate Points", None))
        self.triangualtePointsPushButton.setText(QCoreApplication.translate("TrianglePPDemoAppClass", u"Triangulate", None))
        self.useConstraintsCheckBox.setText(QCoreApplication.translate("TrianglePPDemoAppClass", u"Quality constraints", None))
        self.hideMarkersCheckBox.setText(QCoreApplication.translate("TrianglePPDemoAppClass", u"Hide markers", None))
        self.tesselatePointsPushButton.setText(QCoreApplication.translate("TrianglePPDemoAppClass", u"Tesselate", None))
#if QT_CONFIG(tooltip)
        self.optionsToolButton.setToolTip(QCoreApplication.translate("TrianglePPDemoAppClass", u"Settings", None))
#endif // QT_CONFIG(tooltip)
        self.optionsToolButton.setText(QCoreApplication.translate("TrianglePPDemoAppClass", u"...", None))
    # retranslateUi

