# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'TrianglePpOptions.ui'
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
from PySide6.QtWidgets import (QAbstractButton, QApplication, QCheckBox, QDialog,
    QDialogButtonBox, QGridLayout, QLabel, QLineEdit,
    QPushButton, QSizePolicy, QSpacerItem, QVBoxLayout,
    QWidget)

class Ui_TrianglePpOptions(object):
    def setupUi(self, TrianglePpOptions):
        if not TrianglePpOptions.objectName():
            TrianglePpOptions.setObjectName(u"TrianglePpOptions")
        TrianglePpOptions.resize(436, 552)
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(TrianglePpOptions.sizePolicy().hasHeightForWidth())
        TrianglePpOptions.setSizePolicy(sizePolicy)
        self.verticalLayout = QVBoxLayout(TrianglePpOptions)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.gridLayout = QGridLayout()
        self.gridLayout.setObjectName(u"gridLayout")
        self.label_17 = QLabel(TrianglePpOptions)
        self.label_17.setObjectName(u"label_17")

        self.gridLayout.addWidget(self.label_17, 19, 1, 1, 1)

        self.label_18 = QLabel(TrianglePpOptions)
        self.label_18.setObjectName(u"label_18")

        self.gridLayout.addWidget(self.label_18, 20, 1, 1, 1)

        self.segmentColorButton = QPushButton(TrianglePpOptions)
        self.segmentColorButton.setObjectName(u"segmentColorButton")
        sizePolicy.setHeightForWidth(self.segmentColorButton.sizePolicy().hasHeightForWidth())
        self.segmentColorButton.setSizePolicy(sizePolicy)

        self.gridLayout.addWidget(self.segmentColorButton, 19, 2, 1, 1)

        self.label_12 = QLabel(TrianglePpOptions)
        self.label_12.setObjectName(u"label_12")
        self.label_12.setIndent(-1)

        self.gridLayout.addWidget(self.label_12, 14, 1, 1, 1)

        self.maxAreaLineEdit = QLineEdit(TrianglePpOptions)
        self.maxAreaLineEdit.setObjectName(u"maxAreaLineEdit")
        sizePolicy1 = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.maxAreaLineEdit.sizePolicy().hasHeightForWidth())
        self.maxAreaLineEdit.setSizePolicy(sizePolicy1)

        self.gridLayout.addWidget(self.maxAreaLineEdit, 3, 2, 1, 1)

        self.constrainedDelaunayCheckBox = QCheckBox(TrianglePpOptions)
        self.constrainedDelaunayCheckBox.setObjectName(u"constrainedDelaunayCheckBox")
        self.constrainedDelaunayCheckBox.setChecked(True)

        self.gridLayout.addWidget(self.constrainedDelaunayCheckBox, 11, 2, 1, 1)

        self.maxPointCountLineEdit = QLineEdit(TrianglePpOptions)
        self.maxPointCountLineEdit.setObjectName(u"maxPointCountLineEdit")
        sizePolicy1.setHeightForWidth(self.maxPointCountLineEdit.sizePolicy().hasHeightForWidth())
        self.maxPointCountLineEdit.setSizePolicy(sizePolicy1)

        self.gridLayout.addWidget(self.maxPointCountLineEdit, 9, 2, 1, 1)

        self.applyQualityConstraintsCheckBox = QCheckBox(TrianglePpOptions)
        self.applyQualityConstraintsCheckBox.setObjectName(u"applyQualityConstraintsCheckBox")

        self.gridLayout.addWidget(self.applyQualityConstraintsCheckBox, 15, 2, 1, 1)

        self.label_11 = QLabel(TrianglePpOptions)
        self.label_11.setObjectName(u"label_11")
        self.label_11.setWordWrap(True)

        self.gridLayout.addWidget(self.label_11, 4, 1, 1, 1)

        self.label_8 = QLabel(TrianglePpOptions)
        self.label_8.setObjectName(u"label_8")
        font = QFont()
        font.setBold(True)
        self.label_8.setFont(font)

        self.gridLayout.addWidget(self.label_8, 10, 0, 1, 2)

        self.horizontalSpacer_2 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.gridLayout.addItem(self.horizontalSpacer_2, 3, 4, 1, 1)

        self.label_2 = QLabel(TrianglePpOptions)
        self.label_2.setObjectName(u"label_2")

        self.gridLayout.addWidget(self.label_2, 8, 1, 1, 1)

        self.horizontalSpacer = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.gridLayout.addItem(self.horizontalSpacer, 2, 4, 1, 1)

        self.minAngleLineEdit = QLineEdit(TrianglePpOptions)
        self.minAngleLineEdit.setObjectName(u"minAngleLineEdit")
        sizePolicy1.setHeightForWidth(self.minAngleLineEdit.sizePolicy().hasHeightForWidth())
        self.minAngleLineEdit.setSizePolicy(sizePolicy1)

        self.gridLayout.addWidget(self.minAngleLineEdit, 2, 2, 1, 1)

        self.voronoiColorButton = QPushButton(TrianglePpOptions)
        self.voronoiColorButton.setObjectName(u"voronoiColorButton")
        sizePolicy.setHeightForWidth(self.voronoiColorButton.sizePolicy().hasHeightForWidth())
        self.voronoiColorButton.setSizePolicy(sizePolicy)

        self.gridLayout.addWidget(self.voronoiColorButton, 20, 2, 1, 1)

        self.label_6 = QLabel(TrianglePpOptions)
        self.label_6.setObjectName(u"label_6")
        self.label_6.setFont(font)

        self.gridLayout.addWidget(self.label_6, 1, 0, 1, 2)

        self.minPointCountLineEdit = QLineEdit(TrianglePpOptions)
        self.minPointCountLineEdit.setObjectName(u"minPointCountLineEdit")
        sizePolicy1.setHeightForWidth(self.minPointCountLineEdit.sizePolicy().hasHeightForWidth())
        self.minPointCountLineEdit.setSizePolicy(sizePolicy1)

        self.gridLayout.addWidget(self.minPointCountLineEdit, 8, 2, 1, 1)

        self.label_4 = QLabel(TrianglePpOptions)
        self.label_4.setObjectName(u"label_4")

        self.gridLayout.addWidget(self.label_4, 9, 1, 1, 1)

        self.label_16 = QLabel(TrianglePpOptions)
        self.label_16.setObjectName(u"label_16")

        self.gridLayout.addWidget(self.label_16, 15, 1, 1, 1)

        self.label_10 = QLabel(TrianglePpOptions)
        self.label_10.setObjectName(u"label_10")

        self.gridLayout.addWidget(self.label_10, 11, 1, 1, 1)

        self.horizontalSpacer_5 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.gridLayout.addItem(self.horizontalSpacer_5, 8, 0, 1, 1)

        self.vertexColorButton = QPushButton(TrianglePpOptions)
        self.vertexColorButton.setObjectName(u"vertexColorButton")
        sizePolicy.setHeightForWidth(self.vertexColorButton.sizePolicy().hasHeightForWidth())
        self.vertexColorButton.setSizePolicy(sizePolicy)

        self.gridLayout.addWidget(self.vertexColorButton, 18, 2, 1, 1)

        self.conformingDelaunayCheckBox = QCheckBox(TrianglePpOptions)
        self.conformingDelaunayCheckBox.setObjectName(u"conformingDelaunayCheckBox")

        self.gridLayout.addWidget(self.conformingDelaunayCheckBox, 13, 2, 1, 1)

        self.label_14 = QLabel(TrianglePpOptions)
        self.label_14.setObjectName(u"label_14")

        self.gridLayout.addWidget(self.label_14, 31, 1, 1, 1)

        self.removeConcavitiesCheckBox = QCheckBox(TrianglePpOptions)
        self.removeConcavitiesCheckBox.setObjectName(u"removeConcavitiesCheckBox")
        self.removeConcavitiesCheckBox.setChecked(False)

        self.gridLayout.addWidget(self.removeConcavitiesCheckBox, 14, 2, 1, 1)

        self.label_7 = QLabel(TrianglePpOptions)
        self.label_7.setObjectName(u"label_7")

        self.gridLayout.addWidget(self.label_7, 2, 3, 1, 1)

        self.horizontalSpacer_3 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.gridLayout.addItem(self.horizontalSpacer_3, 8, 4, 1, 1)

        self.segmentPointsLineEdit = QLineEdit(TrianglePpOptions)
        self.segmentPointsLineEdit.setObjectName(u"segmentPointsLineEdit")

        self.gridLayout.addWidget(self.segmentPointsLineEdit, 4, 2, 1, 2)

        self.label_5 = QLabel(TrianglePpOptions)
        self.label_5.setObjectName(u"label_5")
        self.label_5.setFont(font)

        self.gridLayout.addWidget(self.label_5, 7, 0, 1, 2)

        self.horizontalSpacer_4 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.gridLayout.addItem(self.horizontalSpacer_4, 9, 4, 1, 1)

        self.label_9 = QLabel(TrianglePpOptions)
        self.label_9.setObjectName(u"label_9")

        self.gridLayout.addWidget(self.label_9, 13, 1, 1, 1)

        self.label_15 = QLabel(TrianglePpOptions)
        self.label_15.setObjectName(u"label_15")

        self.gridLayout.addWidget(self.label_15, 18, 1, 1, 1)

        self.label_3 = QLabel(TrianglePpOptions)
        self.label_3.setObjectName(u"label_3")

        self.gridLayout.addWidget(self.label_3, 3, 1, 1, 1)

        self.seperateSegmentColorCheckBox = QCheckBox(TrianglePpOptions)
        self.seperateSegmentColorCheckBox.setObjectName(u"seperateSegmentColorCheckBox")
        self.seperateSegmentColorCheckBox.setChecked(True)

        self.gridLayout.addWidget(self.seperateSegmentColorCheckBox, 31, 2, 1, 1)

        self.label_13 = QLabel(TrianglePpOptions)
        self.label_13.setObjectName(u"label_13")
        self.label_13.setFont(font)

        self.gridLayout.addWidget(self.label_13, 16, 0, 1, 2)

        self.label = QLabel(TrianglePpOptions)
        self.label.setObjectName(u"label")

        self.gridLayout.addWidget(self.label, 2, 1, 1, 1)

        self.restoreColorsButton = QPushButton(TrianglePpOptions)
        self.restoreColorsButton.setObjectName(u"restoreColorsButton")
        sizePolicy2 = QSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.restoreColorsButton.sizePolicy().hasHeightForWidth())
        self.restoreColorsButton.setSizePolicy(sizePolicy2)
        self.restoreColorsButton.setMinimumSize(QSize(0, 0))
        self.restoreColorsButton.setLayoutDirection(Qt.LayoutDirection.RightToLeft)

        self.gridLayout.addWidget(self.restoreColorsButton, 21, 2, 1, 2)


        self.verticalLayout.addLayout(self.gridLayout)

        self.verticalSpacer = QSpacerItem(20, 44, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.verticalLayout.addItem(self.verticalSpacer)

        self.buttonBox = QDialogButtonBox(TrianglePpOptions)
        self.buttonBox.setObjectName(u"buttonBox")
        self.buttonBox.setOrientation(Qt.Orientation.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.StandardButton.Cancel|QDialogButtonBox.StandardButton.Ok)

        self.verticalLayout.addWidget(self.buttonBox)


        self.retranslateUi(TrianglePpOptions)
        self.buttonBox.accepted.connect(TrianglePpOptions.accept)
        self.buttonBox.rejected.connect(TrianglePpOptions.reject)

        QMetaObject.connectSlotsByName(TrianglePpOptions)
    # setupUi

    def retranslateUi(self, TrianglePpOptions):
        TrianglePpOptions.setWindowTitle(QCoreApplication.translate("TrianglePpOptions", u"Triangle++ Options", None))
        self.label_17.setText(QCoreApplication.translate("TrianglePpOptions", u"Segment color", None))
        self.label_18.setText(QCoreApplication.translate("TrianglePpOptions", u"Voronoi color", None))
        self.segmentColorButton.setText("")
        self.label_12.setText(QCoreApplication.translate("TrianglePpOptions", u"Remove Concavities:", None))
        self.constrainedDelaunayCheckBox.setText("")
        self.applyQualityConstraintsCheckBox.setText("")
        self.label_11.setText(QCoreApplication.translate("TrianglePpOptions", u"Segments (by their start/end point indexes)", None))
        self.label_8.setText(QCoreApplication.translate("TrianglePpOptions", u"Triangulation Type:", None))
        self.label_2.setText(QCoreApplication.translate("TrianglePpOptions", u"Min. Point Count: ", None))
        self.voronoiColorButton.setText("")
        self.label_6.setText(QCoreApplication.translate("TrianglePpOptions", u"Quality and Constraints:", None))
        self.label_4.setText(QCoreApplication.translate("TrianglePpOptions", u"Max. Point Count: ", None))
        self.label_16.setText(QCoreApplication.translate("TrianglePpOptions", u"Apply Quality Constraints", None))
        self.label_10.setText(QCoreApplication.translate("TrianglePpOptions", u"Constrained Delaunay: ", None))
        self.vertexColorButton.setText("")
        self.conformingDelaunayCheckBox.setText("")
        self.label_14.setText(QCoreApplication.translate("TrianglePpOptions", u"Different color for segments:", None))
        self.removeConcavitiesCheckBox.setText("")
        self.label_7.setText(QCoreApplication.translate("TrianglePpOptions", u"[deg]", None))
        self.label_5.setText(QCoreApplication.translate("TrianglePpOptions", u"Point Generation:", None))
        self.label_9.setText(QCoreApplication.translate("TrianglePpOptions", u"Conforming Delaunay: ", None))
        self.label_15.setText(QCoreApplication.translate("TrianglePpOptions", u"Point and vertex color", None))
        self.label_3.setText(QCoreApplication.translate("TrianglePpOptions", u"Max. Area: ", None))
        self.seperateSegmentColorCheckBox.setText("")
        self.label_13.setText(QCoreApplication.translate("TrianglePpOptions", u"View Options:", None))
        self.label.setText(QCoreApplication.translate("TrianglePpOptions", u"Min. Angle: ", None))
        self.restoreColorsButton.setText(QCoreApplication.translate("TrianglePpOptions", u"Restore Default Colors", None))
    # retranslateUi

