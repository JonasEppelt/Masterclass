import matplotlib.pyplot as plt
from itertools import product
import numpy as np

class B_field:
    def __init__(self, B_field, lim) -> None:
        self._lim = lim
        self._B_field = B_field
        self._grid = self.get_grid()
    
    @property
    def field_sign(self):
        if self._B_field == 0:
            return 0
        return self._B_field/abs(self._B_field)

    @property
    def field_marker(self):
        if self.field_sign == 0:
            return
        elif self.field_sign ==-1:
            return "x"
        return "o"
    
    @property
    def marker_size(self):
        return abs(self._B_field)*100

    @property 
    def marker_sizes(self):
        return np.array([self.marker_size]*len(self._grid))

    def get_grid(self):
        grid_x = np.linspace(-self._lim, self._lim,10)
        grid_y = np.linspace(-self._lim, self._lim,10)
        grid = np.array(list(product(grid_x, grid_y))).T
        return grid
    
    def plot_field(self,ax):
        return ax.scatter(self._grid[0], self._grid[1], marker = self.field_marker, s = self.marker_size, color = np.array([1,0,0,0.5]))