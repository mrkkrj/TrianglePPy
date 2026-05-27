from PySide6.QtWidgets import (QMainWindow, QMessageBox, QFileDialog, QMenu)
from PySide6.QtGui import QKeySequence, QAction, QFont, QColor
from PySide6.QtCore import Qt, QPointF, QSettings, QFileInfo, QSignalBlocker
import random
import sys

try:
    import triangle_ppy as tpp  # Your Triangle++ Python binding
except ImportError:
    tpp = None
    print("Warning: tpp module not available. Triangulation disabled.")

from drawing_area import DrawingArea
from ui_triangle_pp_demo_app import Ui_TrianglePPDemoAppClass
from triangle_pp_options import TrianglePpOptions


class TrianglePPDemoApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.ui = Ui_TrianglePPDemoAppClass()
        self.ui.setupUi(self)

        self.ui.drawAreaWidget.setDrawMode(DrawingArea.DrawMode.DrawPoints)
        self.ui.optionsToolButton.setText("☰")

        # State variables
        self.mode = 1  # AutomaticMode
        self.use_quality_constr = False
        self.triangulated = False
        self.tesselated = False
        self.min_angle = -1
        self.max_area = -1
        self.min_points = -1
        self.max_points = -1
        self.use_conforming_delaunay = False
        self.include_convex_hull = True
        self.seperate_segment_color = True
        self.last_file_dir = "."

        self.voronoi_points = []
        self.segment_endpoint_indexes = []
        self.hole_points = []
        self.region_points = []
        self.region_max_areas = []

        self.read_from_file = False
        self.scale_factor = 1.0
        self.offset_x = 0.0
        self.offset_y = 0.0
        self.zoom_factor = 1.0
        self.original_size = None

        self.vertex_points_orig = []
        self.hole_points_orig = []
        self.region_points_orig = []
        self.region_max_areas_orig = []

        self.triangle_color = QColor(Qt.GlobalColor.blue)
        self.voronoi_color = QColor(Qt.GlobalColor.red)
        self.segment_color = QColor("limegreen")
        self.hole_marker_color = QColor("limegreen")

        self.read_from_settings()
        self._add_ui_shortcuts()
        self._connect_signals()

        self.ui.pointModeComboBox.setCurrentIndex(self.mode)
        self._set_generate_button_text()

    def _connect_signals(self):
        self.ui.generatePointsPushButton.clicked.connect(self.on_generate_points)
        self.ui.triangualtePointsPushButton.clicked.connect(self.on_triangulate)
        self.ui.tesselatePointsPushButton.clicked.connect(self.on_tesselate)
        self.ui.pointModeComboBox.currentIndexChanged.connect(self.on_point_mode_changed)
        self.ui.useConstraintsCheckBox.toggled.connect(self.on_use_constraints_toggled)
        self.ui.hideMarkersCheckBox.toggled.connect(self.on_hide_markers_toggled)
        self.ui.optionsToolButton.clicked.connect(self.on_options_clicked)

        da = self.ui.drawAreaWidget
        da.pointDeleted.connect(self.on_triangulation_point_deleted)
        da.linePointsSelected.connect(self.on_segment_endpoints_selected)
        da.pointChangedToHoleMarker.connect(self.on_point_changed_to_hole_marker)
        da.pointMoved.connect(self.on_triangulation_point_moved)

    # ====================== Public Slots ======================

    def on_generate_points(self):
        self.clear_display()
        self.reset_zoom()
        self.read_from_file = False

        if self.mode == 0:    # Manual
            pass
        elif self.mode == 1:  # Automatic
            self._generate_random_points()
        elif self.mode == 3:  # From File
            self.read_from_file()
        elif self.mode == 4:  # Example 1
            self.show_example1()
        elif self.mode == 5:  # Example 2
            self.show_example2()

    def on_triangulate(self):
        self.clear_voronoi_points()
        drawn_points = self.ui.drawAreaWidget.getPointCoordinates()
        self.ui.drawAreaWidget.clearImage()
        self.ui.drawAreaWidget.setDrawColor(self.triangle_color)

        for p in drawn_points:
            self.ui.drawAreaWidget.drawPoint(p)

        if len(drawn_points) < 3:
            QMessageBox.critical(self, "Triangle++", "Not enough points to triangulate!")
            return

        if tpp is None:
            QMessageBox.critical(self, "Error", "Triangle++ binding not available.")
            return

        input_points = [(p.x(), p.y()) for p in drawn_points]

        try:
            delaunay = tpp.Delaunay(input_points)
            self._config_delaunay(delaunay)

            if self.use_conforming_delaunay:
                delaunay.TriangulateConf(self.use_quality_constr)
            else:
                delaunay.Triangulate(self.use_quality_constr)

            self.triangulated = True
            self.tesselated = False
            self._draw_triangulation(delaunay, drawn_points)
        except Exception as e:
            QMessageBox.critical(self, "Triangle++", f"Triangulation error:\n{str(e)}")

    def on_tesselate(self):
        self.clear_voronoi_points()
        drawn_points = self.ui.drawAreaWidget.getPointCoordinates()

        if len(drawn_points) < 3:
            QMessageBox.critical(self, "Triangle++", "Not enough points!")
            return

        if tpp is None:
            return

        try:
            delaunay = tpp.Delaunay([(p.x(), p.y()) for p in drawn_points])
            delaunay.Tesselate(self.use_conforming_delaunay)
            self.tesselated = True
            self._draw_voronoi_tesselation(delaunay)
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))

    def on_point_mode_changed(self, index):
        self.mode = index
        self._set_generate_button_text()
        if self.mode in (3, 4, 5):
            self.on_generate_points()

    def on_use_constraints_toggled(self, checked):
        self.use_quality_constr = checked
        if self.triangulated:
            self.on_triangulate()

    def on_hide_markers_toggled(self, checked):
        color = Qt.GlobalColor.white if checked else self.hole_marker_color
        for hole in self.hole_points:
            self._draw_marker_point(hole, color, "H")
        self.ui.drawAreaWidget.setDrawColor(self.triangle_color)

    def on_options_clicked(self):
        menu = QMenu(self)
        menu.addAction("Save to Image (Ctrl+I)", self.save_to_image)
        menu.addAction("Save to File (Ctrl+S)", self.write_to_file)
        menu.addAction("Options (Ctrl+O)", self.show_triangulation_options)
        menu.addAction("Info", self.show_info)
        menu.addAction("Zoom In", self.zoom_in)
        menu.addAction("Zoom Out", self.zoom_out)
        menu.addAction("Undo Last Point (Ctrl+Z)", self.undo_point_creation)
        menu.addAction("Quit", self.close)
        menu.exec(self.ui.optionsToolButton.mapToGlobal(
            self.ui.optionsToolButton.rect().bottomLeft()))

    # ====================== Drawing Methods ======================

    def _draw_triangulation(self, delaunay, points_on_screen):
        """Draw all triangles"""
        self.ui.drawAreaWidget.setDrawColor(self.triangle_color)

        for face in delaunay.faces():  # Assuming iterator interface
            p1 = self._get_result_point(delaunay, face.org())
            p2 = self._get_result_point(delaunay, face.dest())
            p3 = self._get_result_point(delaunay, face.apex())

            self.ui.drawAreaWidget.drawLine(p1, p2)
            self.ui.drawAreaWidget.drawLine(p2, p3)
            self.ui.drawAreaWidget.drawLine(p3, p1)

        # Draw constraint segments
        self.ui.drawAreaWidget.setDrawColor(self.segment_color()())
        for i in range(0, len(self.segment_endpoint_indexes), 2):
            start = self.segment_endpoint_indexes[i]
            end = self.segment_endpoint_indexes[i + 1]
            if start < len(points_on_screen) and end < len(points_on_screen):
                self.ui.drawAreaWidget.drawLine(points_on_screen[start], points_on_screen[end])

        # Draw markers
        self.on_hide_markers_toggled(self.ui.hideMarkersCheckBox.isChecked())

        self.statusBar().showMessage(f"Created {delaunay.triangleCount()} triangles")

    def _draw_voronoi_tesselation(self, delaunay):
        """Draw Voronoi diagram"""
        self.ui.drawAreaWidget.setDrawColor(self.voronoi_color)

        # Voronoi vertices
        for v in delaunay.voronoi_vertices():
            pt = QPointF(v[0], v[1])
            self.ui.drawAreaWidget.drawPoint(pt)
            self.voronoi_points.append(pt)

        # Voronoi edges
        for edge in delaunay.voronoi_edges():
            p1 = QPointF(edge.org()[0], edge.org()[1])
            if edge.is_finite():
                p2 = QPointF(edge.dest()[0], edge.dest()[1])
                self.ui.drawAreaWidget.drawLine(p1, p2)
            else:
                # Infinite ray - draw long line
                dx = edge.dest()[0]
                dy = edge.dest()[1]
                p2 = QPointF(p1.x() + dx * 100, p1.y() + dy * 100)
                self.ui.drawAreaWidget.drawLine(p1, p2)

        self.ui.drawAreaWidget.setDrawColor(self.triangle_color)

    def _config_delaunay(self, delaunay):
        """Configure constraints"""
        if self.use_quality_constr:
            if self.min_angle > 0:
                delaunay.setMinAngle(self.min_angle)
            if self.max_area > 0:
                delaunay.setMaxArea(self.max_area)

        delaunay.useConvexHullWithSegments(self.include_convex_hull)

        if self.segment_endpoint_indexes:
            delaunay.setSegmentConstraint(self.segment_endpoint_indexes)

        if self.hole_points:
            holes = [(p.x(), p.y()) for p in (self.hole_points_orig if self.read_from_file else self.hole_points)]
            delaunay.setHolesConstraint(holes)

    def _get_result_point(self, delaunay, idx):
        """Helper for Steiner points"""
        if idx == -1:  # Steiner point
            return QPointF(0, 0)  # Replace with actual logic from your binding
        x, y = delaunay.pointAtVertexId(idx)
        return QPointF(x, y)

    # ====================== Event Handlers ======================

    def on_triangulation_point_deleted(self, pos: QPointF):
        if pos in self.hole_points:
            self.hole_points.remove(pos)
        if self.triangulated:
            self.on_triangulate()

    def on_segment_endpoints_selected(self, start_idx: int, end_idx: int):
        self.ui.drawAreaWidget.setDrawColor(self.segment_color())
        points = self.ui.drawAreaWidget.getPointCoordinates()
        if start_idx < len(points) and end_idx < len(points):
            self.ui.drawAreaWidget.drawLine(points[start_idx], points[end_idx])
            self.segment_endpoint_indexes.extend([start_idx, end_idx])
        self.ui.drawAreaWidget.setDrawColor(self.triangle_color)

    def on_point_changed_to_hole_marker(self, point_idx: int, pos: QPointF):
        self._draw_marker_point(pos, self.hole_marker_color, "H")
        self.hole_points.append(pos)
        self.ui.hideMarkersCheckBox.show()
        self.ui.drawAreaWidget.setDrawColor(self.triangle_color)

    def on_triangulation_point_moved(self, pos1: QPointF, pos2: QPointF):
        if self.triangulated:
            self.on_triangulate()
        if self.tesselated:
            self.on_tesselate()

    def _draw_marker_point(self, pos: QPointF, color: QColor, text: str):
        old_mode = self.ui.drawAreaWidget.getDrawMode()
        self.ui.drawAreaWidget.setDrawColor(color)
        font = QFont()
        font.setPixelSize(22)
        self.ui.drawAreaWidget.drawText(pos * 0.96, text, font)
        self.ui.drawAreaWidget.setDrawMode(DrawingArea.DrawMode.DrawHoleMarker)
        self.ui.drawAreaWidget.drawPoint(pos)
        self.ui.drawAreaWidget.setDrawMode(old_mode)

    # ====================== Utility Methods ======================

    def _generate_random_points(self):
        min_p = max(3, self.min_points)
        max_p = max(100, self.max_points)
        count = random.randint(min_p, max_p)

        w = self.ui.drawAreaWidget.width()
        h = self.ui.drawAreaWidget.height()
        ps = self.ui.drawAreaWidget.getPointSize()

        for _ in range(count):
            x = random.randint(ps // 2, w - ps // 2)
            y = random.randint(ps // 2, h - ps // 2)
            self.ui.drawAreaWidget.drawPoint(QPointF(x, y))

    def show_example1(self):
        pts = [(0,0),(0,1),(0,3),(2,0),(4,1.25),(4,3),(6,0),(8,1.25),(9,0),(9,0.75),(9,3)]
        for x, y in pts:
            self.ui.drawAreaWidget.drawPoint(QPointF(x * 60 + 30, y * 60 + 30))

    def show_example2(self):
        # Add your own example points here
        pass

    def show_triangulation_options(self):
        dlg = TrianglePpOptions(self)
        dlg.fillContents(
            self.min_angle if self.min_angle >= 0 else 20,
            self.max_area,
            self.min_points if self.min_points >= 0 else 3,
            self.max_points if self.max_points >= 0 else 100,
            self.use_quality_constr, self.use_conforming_delaunay,
            self.include_convex_hull, self.seperate_segment_color
        )
        dlg.fillColors(self.triangle_color, self.segment_color, self.voronoi_color)

        if dlg.exec() == 1:  # Accepted
            self.min_angle = dlg.getMinAngle()
            self.max_area = dlg.getMaxArea()
            self.min_points = dlg.getMinPointCount()
            self.max_points = dlg.getMaxPointCount()
            self.use_quality_constr = dlg.applyQualityConstraints()
            self.use_conforming_delaunay = dlg.useConformingDelaunay()
            self.include_convex_hull = dlg.includeConvexHull()
            self.seperate_segment_color = dlg.seperateSegmentColor()
            self.triangle_color, self.segment_color, self.voronoi_color = dlg.getDelaunayColors()
            self.hole_marker_color = self.segment_color
            self.write_to_settings()
            self.ui.useConstraintsCheckBox.setChecked(self.use_quality_constr)

    def show_info(self):
        QMessageBox.about(self, "Triangle++ Demo (Python)",
                          "PySide6 port of Marek Krajewski's Triangle++ Demo\n\n"
                          "Supports Delaunay, CDT, and Voronoi.")

    def clear_display(self):
        self.ui.drawAreaWidget.clearImage()
        self.segment_endpoint_indexes.clear()
        self.hole_points.clear()
        self.region_points.clear()
        self.triangulated = self.tesselated = False

    def clear_voronoi_points(self):
        for p in self.voronoi_points:
            self.ui.drawAreaWidget.clearPoint(p)
        self.voronoi_points.clear()

    def segment_color(self):
        return self.segment_color if self.seperate_segment_color else self.triangle_color

    def zoom_in(self):
        self.zoom_factor += 0.25
        self._zoom_points(1.25)

    def zoom_out(self):
        self.zoom_factor = max(0.1, self.zoom_factor - 0.25)
        self._zoom_points(0.75)

    def undo_point_creation(self):
        self.ui.drawAreaWidget.clearLastPoint()

    def save_to_image(self):
        fn, _ = QFileDialog.getSaveFileName(self, "Save Image", self.last_file_dir,
                                            "PNG (*.png);;JPEG (*.jpg)")
        if fn:
            self.ui.drawAreaWidget.saveImage(fn)

    def write_to_file(self):
        QMessageBox.information(self, "Info", "Save to .node/.poly coming soon.")

    def read_from_file(self):
        fn, _ = QFileDialog.getOpenFileName(self, "Open Triangle File", self.last_file_dir,
                                            "Triangle Files (*.node *.poly)")
        if fn:
            QMessageBox.information(self, "Info", f"Loading {fn} - full support requires tpp binding.")

    def reset_zoom(self):
        if self.zoom_factor == 1.0:
            return
        self.zoom_factor = 1.0
        if self.original_size:
            self.ui.drawAreaWidget.resize(self.original_size)

    def _zoom_points(self, factor):
        # Implementation similar to C++ version
        pass  # Expand if needed

    def _set_generate_button_text(self):
        texts = ["Clear Points", "Generate Points", "Find Points", "Open File",
                 "Read Example", "Read Example"]
        self.ui.generatePointsPushButton.setText(texts[self.mode])

    def _add_ui_shortcuts(self):
        for action in [
            QAction("Zoom In", self, shortcut=QKeySequence.ZoomIn, triggered=self.zoom_in),
            QAction("Zoom Out", self, shortcut=QKeySequence.ZoomOut, triggered=self.zoom_out),
            QAction("Save", self, shortcut=QKeySequence.Save, triggered=self.write_to_file),
            QAction("Undo", self, shortcut=QKeySequence.Undo, triggered=self.undo_point_creation),
        ]:
            self.addAction(action)

    def read_from_settings(self):
        s = QSettings("IB-Krajewski", "Triangle++ Demo")
        self.min_angle = s.value("minAngle", -1, int)

    def write_to_settings(self):
        s = QSettings("IB-Krajewski", "Triangle++ Demo")
        s.setValue("minAngle", self.min_angle)

