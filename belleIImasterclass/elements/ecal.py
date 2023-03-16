from matplotlib.pyplot import sca
import pandas as pd
from copy import deepcopy
import numpy as np


from matplotlib.colors import to_rgba_array
from matplotlib.patches import Rectangle

class ECal:
    def __init__(self, crystal_size = 4, linewidth = 1.75) -> None:
        self._crystal_size = crystal_size
        self._linewidth = linewidth

        self._barrel_rows = 46
        self._barrel_columns = 144
        self._barrel_center = np.array([0,-100])
        barrel_coordinates, barrel_angles, barrel_offsets = self.get_barrel_patches

        self._forward_center = np.array([-300,200])
        self._forward_number_of_rings = 13
        self._forward_rings_sizes = np.array([48,48,64,64,64,96,96,96,96,96,96,144,144])
        self._forward_ring_radiuses = np.linspace(49.62,125.82,13)
        forward_coordinates, forward_angles, forward_offsets = self.get_caps_patches(self._forward_center, self._forward_number_of_rings, self._forward_rings_sizes, self._forward_ring_radiuses, )


        self._backward_center = np.array([300,200])
        self._backward_number_of_rings = 10
        self._backward_rings_sizes =  np.array([144,144,96,96,96,96,96,64,64,64])
        self._backward_ring_radiuses = np.linspace(129.2,58.51,10)
        backward_coordinates, backward_angles, backward_offsets = self.get_caps_patches(self._backward_center, self._backward_number_of_rings, self._backward_rings_sizes, self._backward_ring_radiuses, )

        self._patch_coordinates = np.append(np.append(forward_coordinates, barrel_coordinates, axis = 1), backward_coordinates, axis=1)
        self._patch_angles = np.append(np.append(forward_angles, barrel_angles), backward_angles)
        self._patch_offsets = np.append(np.append(forward_offsets, barrel_offsets, axis = 1), backward_offsets, axis = 1)
        
        self._n_patches = 8736
        self.patches = [Rectangle(
            (self._patch_coordinates[:,n]+self._patch_offsets[:,n]),
            width = self._crystal_size,
            height = self._crystal_size,
            angle = self._patch_angles[n]*180/np.pi,
            linewidth = 10,
        ) for n in range(self._n_patches)]
    
    @property
    def get_barrel_patches(self):
        x_positions = np.linspace(-np.pi*139.1, np.pi*139.1, self._barrel_columns)
        y_positions = np.linspace(-153.6,153.6, 46)
        coordinates = np.empty([2,0])
        angles = np.empty([0])
        offsets = np.empty([2,0])
        for row in range(self._barrel_rows):
            coordinates = np.append(coordinates, np.full((self._barrel_columns,2), self._barrel_center).T+np.array([x_positions, np.ones(self._barrel_columns)*y_positions[row]]), axis = 1)
            angles = np.append(angles, np.zeros(144))
            offsets = np.append(offsets, -np.ones([2,self._barrel_columns])*self._crystal_size, axis = 1)
        return coordinates, angles, offsets
    
    def get_caps_patches(self, cap_center, number_of_rings, ring_sizes, ring_radiuses):

        coordinates = np.empty([2,0])
        angles = np.empty([0])
        offsets = np.empty([2,0])
        for n_ring in range(number_of_rings):
            phi = np.linspace(0,2*np.pi, ring_sizes[n_ring], endpoint=False)  #- np.pi/2
            angles = np.append(angles, phi)
            coordinates = np.append(coordinates, np.full((ring_sizes[n_ring], 2), cap_center).T+ring_radiuses[n_ring]*np.array([np.cos(phi), np.sin(phi)]), axis = 1)
            offsets = np.append(offsets, -np.sqrt(2)*self._crystal_size*np.array([np.sin(np.pi/4-phi), np.cos(np.pi/4-phi)])/2, axis = 1)
        return coordinates, angles, offsets
