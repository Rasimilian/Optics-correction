import sys
import numpy as np
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtWidgets import QApplication, QWidget
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import QObject, pyqtSignal
import pyqtgraph as pg
import paramiko
import epics
import plotly
from plotly.graph_objects import Figure, Scatter
# import PyQtWebEngine
from PyQt5.QtWebEngineWidgets import QWebEngineView
import plotly.offline as po
import plotly.graph_objs as go
from typing import List, Tuple, Dict
import sys
from datetime import datetime

import design2
from GUI.first_tab.ORM import *
from GUI.output_writer import EmittingStream


class Application(QtWidgets.QMainWindow, design2.Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.initUI()
        self.is_structure_loaded = False

    def initUI(self):
        self.setWindowTitle('LOCO')
        self.resize(800, 600)
        self.setWindowIcon(QIcon('GUI/binp.gif'))
        self.show()
        sys.stdout = EmittingStream(textWritten=self.redirect_output_text)
        print('Hello!')

        # Load Structures
        self.actionOpen.triggered.connect(self.load_structure)
        # self.actionOpen_2.triggered.connect(self.load_structure)

        # Plot ideal ORM
        self.pushButton.clicked.connect(lambda: self.plot_ORM(self.ideal_structure,
                                                              "madx/elements/quads.txt",
                                                              None,
                                                              "madx/correctors/correctors.txt",
                                                              1e-6,
                                                              self.gridLayout_3,
                                                              (0, 0, 1, 1)))

        # Plot real ORM
        self.pushButton_4.clicked.connect(lambda: self.plot_ORM("madx/structures/VEPP4M_full1_all_grads_errors.txt",
                                                                "madx/elements/quads.txt",
                                                                None,
                                                                "madx/correctors/correctors.txt",
                                                                1e-6,
                                                                self.gridLayout_3,
                                                                (0, 2, 1, 1)))


        # Plot BPM data
        # self.ydata = epics.caget('VEPP4:NEP0:turns_x-I')
        # self.xdata = list(range(len(self.ydata)))
        # self.canvas = MplCanvas(width=5, height=4, dpi=100)
        # self.gridLayout.addWidget(self.canvas)
        # self.setLayout(self.gridLayout)
        #
        # self.timer = QtCore.QTimer()
        # self.timer.setInterval(100)
        # self.timer.timeout.connect(self.update_plot)
        # self.timer.start()

    def update_plot(self):
        # self.ydata = epics.caget('VEPP4:NEP0:turns_x-I') + np.random.random(len(self.ydata))
        self.ydata = epics.caget('VEPP4:NEP0:turns_x-I')
        self.canvas.axes.cla()
        self.canvas.axes.plot(self.xdata, self.ydata, 'r')
        self.canvas.draw()

    def plot_ORM(self,
                 structure: str,
                 elements_to_vary: List[str],
                 accumulative_param_additive: np.ndarray,
                 correctors: List[str],
                 corrector_step: float,
                 layout: QtWidgets.QGridLayout,
                 position: Tuple):
        print('aa')
        matrix = OrbitResponseMatrix(structure,
                                     elements_to_vary,
                                     accumulative_param_additive,
                                     correctors,
                                     corrector_step)
        matrix.plot(layout, position)

    def load_structure(self):
        structure_name = QtWidgets.QFileDialog.getOpenFileName(self,'Open file','C:/Users/r_mam/PycharmProjects/optics_correction/madx/structures')[0]

        if self.is_structure_loaded == False:
            self.ideal_structure = structure_name
            self.is_structure_loaded = True
        else:
            print("Structure is already loaded!")

    def __del__(self):
        # Restore sys.stdout
        sys.stdout = sys.__stdout__

    def redirect_output_text(self, text):
        """Append text to the Text Widget."""
        self.textBrowser.append(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+': '+text)
        # cursor = self.textBrowser.textCursor()
        # cursor.movePosition(QtGui.QTextCursor.End)
        # cursor.insertText(str(datetime.now().isoformat())+' '+text)
        # self.textBrowser.setTextCursor(cursor)
        # self.textBrowser.ensureCursorVisible()


def main():
    app = QtWidgets.QApplication(sys.argv)
    window = Application()
    window.show()
    app.exec_()


if __name__ == '__main__':
    main()
