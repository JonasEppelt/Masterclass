from matplotlib.pyplot import sca
import pandas as pd
from copy import deepcopy
import numpy as np


from matplotlib.colors import to_rgba_array
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

class Crystal:
    def __init__(self, x, y, content, edge_length, line_width):
        self._x = x
        self._y = y
        self._content = content
        self._edge_length = edge_length - line_width/ # account for line thickness for better visbility
        self._line_width = line_width

    def get_patch(self):
        return Rectangle((self._x, self._y), self._edge_length, self._edge_length, linewidth = self._line_width)

class ECal:
    def __init__(self, crystal_edge_width = 0.5, noise_rate = 0, linewidth = 1.75) -> None:
        self._crystal_edge_width = crystal_edge_width
        self._noise_rate = noise_rate
        self._linewidth = linewidth
        self._barrel_rows = 46
        self._barrel_columns = 144
        self._barrel_crystal_patches = []
        for c in range(self._barrel_columns):
            for r in range(self._barrel_rows):
                


    @property
    def get