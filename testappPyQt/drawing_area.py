from PySide6.QtWidgets import QWidget, QMenu
from PySide6.QtGui import QPainter, QPen, QBrush, QColor, QImage, QFont, QMouseEvent, QPaintEvent, QResizeEvent, QAction
from PySide6.QtCore import Qt, QPointF, QSize, Signal, QRectF

class DrawingArea(QWidget):
    pointAdded = Signal(QPointF)
    pointDeleted = Signal(QPointF)
    holeMarkerDeleted = Signal(QPointF)
    lineAdded = Signal(QPointF, QPointF)
    lineDeleted = Signal(QPointF, QPointF)
    linePointsSelected = Signal(int, int)
    pointChangedToHoleMarker = Signal(int, QPointF)
    pointMoved = Signal(QPointF, QPointF)
    pointPositionChanged = Signal(QPointF, QPointF)
    drawingAreaResized = Signal()

    class DrawMode:
        DrawPoints = 0
        DrawLines = 1
        DrawHoleMarker = 2
        MovePoint = 3

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAttribute(Qt.WidgetAttribute.WA_StaticContents)

        self.mode = self.DrawMode.DrawLines
        self.pen_width = 1
        self.point_size = 6
        self.pen_color = QColor(Qt.GlobalColor.blue)

        self.line_started = False
        self.start_pos = QPointF()
        self.last_moving_pos = QPointF()
        self.start_pos_index = -1
        self.points = []          # QList<QPointF> -> list
        self.hole_marker_points = []

        self.line_start_point_idx = -1

        self.img = QImage()
        self.img_dirty = False

        self._resize_image(self.img, QSize(800, 600))  # initial size

    # ==================== Public API ====================

    def setDrawMode(self, mode: int):
        self.mode = mode

    def setDrawColor(self, color: QColor):
        self.pen_color = color

    def setPointSize(self, diameter: int):
        self.point_size = diameter

    def setLineWidth(self, width: int):
        self.pen_width = width

    def getDrawMode(self) -> int:
        return self.mode

    def getDrawColor(self) -> QColor:
        return self.pen_color

    def getPointSize(self) -> int:
        return self.point_size

    def getLineWidth(self) -> int:
        return self.pen_width

    def isModified(self) -> bool:
        return self.img_dirty

    def hasPoints(self) -> bool:
        return len(self.points) > 0

    def drawPoint(self, pos: QPointF):
        self._draw_point_at(pos)
        if self.mode == self.DrawMode.DrawPoints:
            self.points.append(pos)
        elif self.mode == self.DrawMode.DrawHoleMarker:
            self.hole_marker_points.append(pos)

    def clearPoint(self, pos: QPointF):
        if self._remove_point(pos, self.points):
            self.pointDeleted.emit(pos)

    def clearLastPoint(self):
        self._delete_point_at_last_pos()

    def drawLine(self, start: QPointF, end: QPointF):
        self.start_pos = start
        self._draw_line_to(end)

    def drawText(self, pos: QPointF, text: str, font: QFont = None):
        painter = QPainter(self.img)
        painter.setPen(QPen(self.pen_color))
        if font:
            painter.setFont(font)
        painter.drawText(pos, text)
        self.update()
        self.img_dirty = True

    def getPointCoordinates(self) -> list[QPointF]:
        return self.points[:]

    def getHoleMarkerCoordinates(self) -> list[QPointF]:
        return self.hole_marker_points[:]

    def openImage(self, fileName: str) -> bool:
        img = QImage()
        if not img.load(fileName):
            return False
        sz = img.size().expandedTo(self.size())
        self._resize_image(img, sz)
        self.img = img
        self.img_dirty = False
        self.update()
        return True

    def saveImage(self, fileName: str, fmt: str = None) -> bool:
        img = self.img.copy()
        self._resize_image(img, self.size())
        ok = img.save(fileName, fmt)
        if ok:
            self.img_dirty = False
        return ok

    def clearImage(self):
        self.img.fill(Qt.GlobalColor.white)
        if self.mode == self.DrawMode.DrawPoints:
            self.points.clear()
            self.hole_marker_points.clear()
        self.line_start_point_idx = -1
        self.img_dirty = True
        self.update()

    # ==================== Protected Events ====================

    def mousePressEvent(self, ev: QMouseEvent):
        if ev.button() == Qt.MouseButton.LeftButton:
            self.start_pos = ev.position()
            if self.mode == self.DrawMode.DrawLines:
                self.line_started = True
            else:
                idx = -1
                if self._point_clicked(ev.position(), idx):
                    self.start_pos = self.points[idx]
                    self.last_moving_pos = self.start_pos
                    self.start_pos_index = idx
                    self.mode = self.DrawMode.MovePoint
                    return

        elif ev.button() == Qt.MouseButton.RightButton:
            if self.mode == self.DrawMode.DrawPoints:
                idx = -1
                if self._point_clicked(ev.position(), idx):
                    self.start_pos = self.points[idx]
                    self._show_point_ctx_menu(ev.position())
                    return

    def mouseMoveEvent(self, ev: QMouseEvent):
        if (ev.buttons() & Qt.MouseButton.LeftButton) and self.line_started:
            self._draw_line_to(ev.position())
        elif (ev.buttons() & Qt.MouseButton.LeftButton) and self.mode == self.DrawMode.MovePoint:
            # erase old
            old_color = self.pen_color
            self.pen_color = QColor(Qt.GlobalColor.white)
            self._draw_point_at(self.last_moving_pos)
            self.pen_color = old_color

            self._draw_point_at(ev.position())
            self.points[self.start_pos_index] = ev.position()
            self.pointMoved.emit(self.last_moving_pos, ev.position())
            self.last_moving_pos = ev.position()

    def mouseReleaseEvent(self, ev: QMouseEvent):
        if ev.button() == Qt.MouseButton.LeftButton:
            if self.mode == self.DrawMode.DrawLines and self.line_started:
                self._draw_line_to(ev.position())
                self.line_started = False
            elif self.mode == self.DrawMode.DrawPoints:
                self._draw_point_at(ev.position())
                self.points.append(ev.position())
            elif self.mode == self.DrawMode.MovePoint:
                self.points[self.start_pos_index] = ev.position()
                self.mode = self.DrawMode.DrawPoints
                self.start_pos_index = -1
                self.pointPositionChanged.emit(self.start_pos, ev.position())

    def paintEvent(self, ev: QPaintEvent):
        painter = QPainter(self)
        dirty = ev.rect()
        painter.drawImage(dirty, self.img, dirty)

    def resizeEvent(self, ev: QResizeEvent):
        if self.width() > self.img.width() or self.height() > self.img.height():
            new_w = max(self.width() + 128, self.img.width())
            new_h = max(self.height() + 128, self.img.height())
            self._resize_image(self.img, QSize(new_w, new_h))
            self.update()
        self.drawingAreaResized.emit()

    # ==================== Private Helpers ====================

    def _draw_line_to(self, end_pos: QPointF):
        painter = QPainter(self.img)
        pen = QPen(self.pen_color, self.pen_width, Qt.PenStyle.SolidLine,
                   Qt.PenCapStyle.RoundCap, Qt.PenJoinStyle.RoundJoin)
        painter.setPen(pen)
        painter.drawLine(self.start_pos, end_pos)

        rad = (self.pen_width // 2) + 2
        update_rect = QRectF(self.start_pos, end_pos).toRect().normalized().adjusted(-rad, -rad, rad, rad)
        self.update(update_rect)

        self.start_pos = end_pos
        self.img_dirty = True

    def _draw_point_at(self, pos: QPointF):
        painter = QPainter(self.img)
        painter.setPen(QPen(self.pen_color, self.pen_width))

        if self.point_size <= 1:
            painter.drawPoint(pos)
        else:
            painter.setBrush(QBrush(self.pen_color))
            r = self.point_size / 2.0
            painter.drawEllipse(pos, r, r)
        self.update()
        self.img_dirty = True

    def _remove_point(self, pos: QPointF, points_list: list) -> bool:
        try:
            idx = points_list.index(pos)
        except ValueError:
            return False

        del points_list[idx]

        # overpaint with white
        old_color = self.pen_color
        self.pen_color = QColor(Qt.GlobalColor.white)
        self._draw_point_at(pos)
        self.pen_color = old_color

        self.img_dirty = True
        return True

    def _resize_image(self, img: QImage, new_size: QSize):
        if img.size() == new_size:
            return
        resized = QImage(new_size, QImage.Format.Format_RGB32)
        resized.fill(Qt.GlobalColor.white)
        painter = QPainter(resized)
        painter.drawImage(QPointF(0, 0), img)
        img.swap(resized)   # modify in place

    def _point_clicked(self, click_pos: QPointF, out_index: list) -> bool:
        """out_index must be a mutable list with one element"""
        half = self.point_size / 2.0
        for i, pt in enumerate(self.points):
            if abs(click_pos.x() - pt.x()) <= half and abs(click_pos.y() - pt.y()) <= half:
                out_index[0] = i
                return True
        return False

    def _show_point_ctx_menu(self, pos: QPointF):
        menu = QMenu(self)
        is_hole = pos in self.hole_marker_points

        act_delete = QAction("Delete HoleMarker" if is_hole else "Delete Point", self)
        act_delete.triggered.connect(self._delete_point_at_last_pos)
        menu.addAction(act_delete)

        # You can extend other actions (Start Segment, etc.) similarly

        menu.exec(self.mapToGlobal(self.start_pos.toPoint()))

    def _delete_point_at_last_pos(self):
        self.clearPoint(self.start_pos)
