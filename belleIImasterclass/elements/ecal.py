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

        self.thetavalues=np.array([13.2,14.7,16.2,17.8,19.2,20.8,22.3,23.7,25.2,26.6,28,29.4,30.8,32.9,34.3,35.7,37.1,38.6,40.2,
                                   41.7,43.4,45.1,46.9,48.6,50.5,52.4,54.4,56.4,58.5,60.6,62.8,65,67.2,69.5,71.9,74.2,76.6,79,
                                   81.6,84,86.5,88.9,91.2,93.5,96,98.5,101,103.4,105.8,108.2,110.5,112.8,115.1,117.3,119.4,
                                   121.5,123.6,125.6,127.6,131.4,133.6,135.8,138.1,140.5,142.9,145.5,148.2,150.7,153.4])*np.pi/180

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
    
    def get_cell_coordinates(self,theta,phi):
        cell_ids=[]
        for i in range(len(theta)):
            _phi=phi[i]
            if _phi < 0:
                _phi=_phi+2*np.pi
            theta_id=np.argmin(abs(theta[i]-self.thetavalues))
            if theta_id<self._forward_number_of_rings:
                id=self._forward_rings_sizes[np.arange(theta_id)].sum()
                id+=int(self._forward_rings_sizes[theta_id]*_phi/(2*np.pi))

            elif theta_id>=self._forward_number_of_rings and theta_id<self._forward_number_of_rings+self._barrel_rows:
                id=self._forward_rings_sizes.sum()
                id+=self._barrel_columns*(theta_id-self._forward_number_of_rings)
                id+=int(self._barrel_columns*_phi/(2*np.pi))
            else:
                id=self._forward_rings_sizes.sum()
                id+=self._barrel_columns*self._barrel_rows
                id+=self._backward_rings_sizes[np.arange(theta_id-self._barrel_rows-self._forward_number_of_rings)].sum()
                id+=int(self._backward_rings_sizes[theta_id-self._barrel_rows-self._forward_number_of_rings]*_phi/(2*np.pi))

            cell_ids.append(id)    
        return self._patch_coordinates[:,cell_ids]+self._patch_offsets[:,cell_ids]
    