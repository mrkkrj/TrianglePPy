import sys
from PySide6.QtWidgets import QApplication
from triangle_pp_demo_app import TrianglePPDemoApp

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = TrianglePPDemoApp()
    window.show()
    sys.exit(app.exec())

