from typing import List, Tuple
import matplotlib
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import plotly.graph_objs as go
import plotly.offline as po
from PyQt5.QtWebEngineWidgets import QWebEngineView

from madx.madx_tool import Structure
from element_parser.data_parser import describe_elements, describe_correctors

matplotlib.use('Qt5Agg')


class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super().__init__(fig)


class OrbitResponseMatrix:
    def __init__(self,
                 structure: str,
                 elements_to_vary: List[str],
                 accumulative_param_additive: np.ndarray,
                 correctors: List[str],
                 corrector_step: float):
        self.structure = structure
        self.elements_to_vary, _= describe_elements(self.structure, elements_to_vary)

        if accumulative_param_additive:
            self.accumulative_param_additive = accumulative_param_additive
        else:
            self.accumulative_param_additive = np.zeros(len(self.elements_to_vary))
        self.correctors, _ = describe_correctors(self.structure, correctors)
        self.corrector_step = corrector_step
        self.matrix = Structure.calculate_response_matrix(self.structure,
                                                          self.elements_to_vary,
                                                          self.accumulative_param_additive,
                                                          self.correctors,
                                                          self.corrector_step)

    def show_qt(self, fig: go.Figure):
        raw_html = '<html><head><meta charset="utf-8" />'
        raw_html += '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script></head>'
        raw_html += '<body>'
        raw_html += po.plot(fig, include_plotlyjs=False, output_type='div')
        raw_html += '</body></html>'

        fig_view = QWebEngineView()
        # setHtml has a 2MB size limit, need to switch to setUrl on tmp file
        # for large figures.

        fig_view.setHtml(raw_html)
        fig_view.show()
        fig_view.raise_()
        return fig_view

    def plot(self, layout, position):
        fig = go.Figure(data=[go.Surface(z=self.matrix)])
        fig_view = self.show_qt(fig)
        layout.addWidget(fig_view, *position)
