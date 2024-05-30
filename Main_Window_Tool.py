import sys
import os
import subprocess
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QVBoxLayout, QLabel, QFrame, QHBoxLayout
from PyQt5.QtGui import QFont  
from PyQt5.QtCore import Qt


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Sequence Alignment Tool")
        self.setGeometry(300, 300, 800, 600)

        vbox = QVBoxLayout()
        self.setLayout(vbox)

        # Heading
        heading_label = QLabel("Sequence Alignment Tool")
        heading_label.setStyleSheet("font-size: 22pt; font-weight: bold")
        heading_label.setAlignment(Qt.AlignCenter)
        vbox.addWidget(heading_label)

        # Line
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        vbox.addWidget(line)

        # Horizontal layout for buttons
        hbox = QHBoxLayout()

        self.pairwise_button = QPushButton("Pairwise Sequence Alignment", self)
        self.pairwise_button.setStyleSheet("background-color: #D8CFC4")  # Set color
        self.pairwise_button.setFont(QFont("Arial", 12))
        self.pairwise_button.setFixedSize(450, 150)  # Set button size
        self.pairwise_button.setCursor(Qt.PointingHandCursor)  # Set cursor
        hbox.addWidget(self.pairwise_button)

        self.multiple_button = QPushButton("Multiple Sequence Alignment", self)
        self.multiple_button.setStyleSheet("background-color: #CCB08F")  # Set color
        self.multiple_button.setFont(QFont("Arial", 12))
        self.multiple_button.setFixedSize(450, 150)  # Set button size
        self.multiple_button.setCursor(Qt.PointingHandCursor)  # Set cursor
        hbox.addWidget(self.multiple_button)

        vbox.addLayout(hbox)  # Add horizontal layout to vertical layout

        # Connect button clicks to respective functions
        self.pairwise_button.clicked.connect(self.open_pairwise_window)
        self.multiple_button.clicked.connect(self.open_multiple_window)

    def open_pairwise_window(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        pairwise_script = os.path.join(script_dir, "Pairwise_Sequence_Alignement_Tool.py")
        subprocess.Popen(["python", pairwise_script])

    def open_multiple_window(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        multiple_script = os.path.join(script_dir, "Multiple_Sequence_Alignement_Tool.py")
        subprocess.Popen(["python", multiple_script])


if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())
