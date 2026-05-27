from PySide6.QtWidgets import QDialog, QColorDialog, QMessageBox 
from PySide6.QtGui import QPalette, QColor, QDoubleValidator
from PySide6.QtCore import Qt, QRegularExpression
from ui_triangle_pp_options import Ui_TrianglePpOptions


class TrianglePpOptions(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = Ui_TrianglePpOptions()
        self.ui.setupUi(self)

        self.vertex_color = QColor(Qt.GlobalColor.blue)
        self.segment_color = QColor("limegreen")
        self.voronoi_color = QColor(Qt.GlobalColor.red)

        self.min_angle_ok = 0.0
        self.min_angle_warning = 0.0
        self.normal_palette = self.palette()

        # Validators
        self.ui.minAngleLineEdit.setValidator(
            QDoubleValidator(0.0, 60.0, 2, self)
        )

        # Connect signals
        self.ui.constrainedDelaunayCheckBox.clicked.connect(self.on_constrained_delaunay_clicked)
        self.ui.conformingDelaunayCheckBox.clicked.connect(self.on_conforming_delaunay_clicked)
        self.ui.applyQualityConstraintsCheckBox.clicked.connect(self.on_apply_quality_clicked)
        self.ui.minAngleLineEdit.textChanged.connect(self.on_min_angle_text_changed)
        self.ui.seperateSegmentColorCheckBox.clicked.connect(self.on_seperate_segment_color_clicked)

        self.ui.vertexColorButton.clicked.connect(self.on_vertex_color_clicked)
        self.ui.segmentColorButton.clicked.connect(self.on_segment_color_clicked)
        self.ui.voronoiColorButton.clicked.connect(self.on_voronoi_color_clicked)
        self.ui.restoreColorsButton.clicked.connect(self.on_restore_colors_clicked)

    # ==================== Public API ====================

    def fillContents(self, min_angle, max_area=-1, min_points=-1, max_points=-1,
                     quality_constr=False, conf_delaunay=False, convex_hull=True,
                     diff_color_segments=True):
        self.ui.minAngleLineEdit.setText(str(min_angle) if min_angle >= 0 else "")
        self.ui.maxAreaLineEdit.setText(str(max_area) if max_area >= 0 else "")
        self.ui.minPointCountLineEdit.setText(str(min_points) if min_points >= 0 else "")
        self.ui.maxPointCountLineEdit.setText(str(max_points) if max_points >= 0 else "")

        self.ui.applyQualityConstraintsCheckBox.setChecked(quality_constr)
        self.ui.conformingDelaunayCheckBox.setChecked(conf_delaunay)
        self.ui.constrainedDelaunayCheckBox.setChecked(not conf_delaunay)
        self._enable_min_max(not conf_delaunay)

        self.ui.removeConcavitiesCheckBox.setChecked(not convex_hull)
        self.ui.seperateSegmentColorCheckBox.setChecked(diff_color_segments)

    def fillColors(self, vertex_color, segment_color, voronoi_color):
        self.vertex_color = vertex_color
        self.segment_color = segment_color
        self.voronoi_color = voronoi_color

        self.ui.vertexColorButton.setPalette(QPalette(vertex_color))

        shown_seg = segment_color if self.ui.seperateSegmentColorCheckBox.isChecked() else vertex_color
        self.ui.segmentColorButton.setPalette(QPalette(shown_seg))
        self.ui.segmentColorButton.setEnabled(self.ui.seperateSegmentColorCheckBox.isChecked())

        self.ui.voronoiColorButton.setPalette(QPalette(voronoi_color))

    def setMinAngleBoundaries(self, max_ok, max_warning):
        self.min_angle_ok = max_ok
        self.min_angle_warning = max_warning

    def setSegmentPointIndexes(self, indexes):
        self.ui.segmentPointsLineEdit.setText(", ".join(map(str, indexes)))

    def getMinAngle(self) -> int:
        return int(self.ui.minAngleLineEdit.text() or 0)

    def getMaxArea(self) -> int:
        return int(self.ui.maxAreaLineEdit.text() or 0)

    def getMinPointCount(self) -> int:
        return int(self.ui.minPointCountLineEdit.text() or 0)

    def getMaxPointCount(self) -> int:
        return int(self.ui.maxPointCountLineEdit.text() or 0)

    def useConformingDelaunay(self) -> bool:
        return self.ui.conformingDelaunayCheckBox.isChecked()

    def applyQualityConstraints(self) -> bool:
        return self.ui.applyQualityConstraintsCheckBox.isChecked()

    def includeConvexHull(self) -> bool:
        return not self.ui.removeConcavitiesCheckBox.isChecked()

    def seperateSegmentColor(self) -> bool:
        return self.ui.seperateSegmentColorCheckBox.isChecked()

    def getSegmentPointIndexes(self):
        text = self.ui.segmentPointsLineEdit.text().strip()
        if not text:
            return []
        return [int(x.strip()) for x in text.split(",")]

    def getDelaunayColors(self):
        return self.vertex_color, self.segment_color, self.voronoi_color

    # ==================== Slots ====================

    def on_constrained_delaunay_clicked(self, checked):
        self.ui.conformingDelaunayCheckBox.setChecked(not checked)
        self._enable_min_max(checked)

    def on_conforming_delaunay_clicked(self, checked):
        self.ui.constrainedDelaunayCheckBox.setChecked(not checked)
        self._enable_min_max(not checked)

    def on_apply_quality_clicked(self, checked):
        self.ui.conformingDelaunayCheckBox.setEnabled(not checked)
        if checked:
            self.ui.conformingDelaunayCheckBox.setChecked(False)
        self._enable_min_max(not checked)

    def on_min_angle_text_changed(self):
        ok = self.ui.minAngleLineEdit.hasAcceptableInput()
        pal = QPalette()
        if not ok:
            pal.setColor(QPalette.ColorRole.Base, Qt.GlobalColor.red)
            pal.setColor(QPalette.ColorRole.Text, Qt.GlobalColor.white)
        else:
            pal = self.normal_palette
        self.ui.minAngleLineEdit.setPalette(pal)
        self.ui.buttonBox.button(QDialog.DialogCode.Accepted).setEnabled(ok)

    def on_seperate_segment_color_clicked(self, _):
        self.fillColors(self.vertex_color, self.segment_color, self.voronoi_color)

    def on_vertex_color_clicked(self):
        color = QColorDialog.getColor(self.vertex_color, self)
        if color.isValid():
            self.vertex_color = color
            self.ui.vertexColorButton.setPalette(QPalette(color))

    def on_segment_color_clicked(self):
        color = QColorDialog.getColor(self.segment_color, self)
        if color.isValid():
            self.segment_color = color
            self.ui.segmentColorButton.setPalette(QPalette(color))

    def on_voronoi_color_clicked(self):
        color = QColorDialog.getColor(self.voronoi_color, self)
        if color.isValid():
            self.voronoi_color = color
            self.ui.voronoiColorButton.setPalette(QPalette(color))

    def on_restore_colors_clicked(self):
        self.vertex_color = QColor(Qt.GlobalColor.blue)
        self.segment_color = QColor("limegreen")
        self.voronoi_color = QColor(Qt.GlobalColor.red)
        self.fillColors(self.vertex_color, self.segment_color, self.voronoi_color)
        self.ui.seperateSegmentColorCheckBox.setChecked(True)

    def _enable_min_max(self, enable):
        self.ui.minAngleLineEdit.setEnabled(enable)
        self.ui.maxAreaLineEdit.setEnabled(enable)
