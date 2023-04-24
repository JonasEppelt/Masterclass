from matplotlib.collections import LineCollection
import pandas as pd
import numpy as np
from belleIImasterclass.particlesmanager import ParticlesManager
from typing import Tuple

import time

class Tracker:
    '''
    Class to design a tracking detector
    '''
    def __init__(self, 
    n_layers: int = 50,
    n_segments: int = 8, 
    add_segments: int = 3, 
    noise_ratio: float = 0.2, 
    granularity: int =100, 
    linewidth: float = 2.5,
    B_field: float = None,
    ) -> None:
        self._n_layers = n_layers
        self._noise_ratio = noise_ratio
        self._granularity = granularity
        self._n_segments = n_segments
        self._add_segments = add_segments
        self._linewidth = linewidth        
        self._B_field = 0.5/self._n_layers if not B_field else B_field
        self._segment_df = pd.DataFrame(
            index = np.arange(0,self.n_seg_total, dtype = int),
            columns = ["begin", "end","radius", "noiseflag","hitflag"],
            dtype = float
            )
        
        self._segment_df.hitflag = False

        segments_counter = 0
        for l in range(1,self._n_layers+1):
            len_segment = 2*np.pi/(self._n_segments+self._add_segments*l) # length of segment in layer l in rad
            offset = l%2
            for i in range(n_segments+l*self._add_segments):
                self._segment_df.loc[segments_counter, "begin"] = len_segment*i+0.1/(l+1)+offset
                self._segment_df.loc[segments_counter, "end"] = len_segment*(i+1)-0.1/(l+1)+offset
                self._segment_df.loc[segments_counter, "radius"] = l
                segments_counter += 1

        theta = np.linspace(self._segment_df["begin"], self._segment_df["end"], self._granularity, dtype = float)
        self.tracker_lines = (self._segment_df["radius"].to_numpy()[None,None,:] * np.array([np.cos(theta),np.sin(theta)])).T

    def get_hits(self, particle, selected=True):
        tracker_selected = "tracker_" if selected else ""
        particle_radius = (particle[f"{tracker_selected}pt"] /abs(self._B_field)).astype(float)
        if particle_radius <= 0.0 or particle[f"{tracker_selected}charge"] == 0:
            return self._segment_df.radius != self._segment_df.radius
        relative_charge = particle[f"{tracker_selected}charge"]*(self._B_field/abs(self._B_field))

        mask0=self._segment_df.radius<=(2*particle_radius)
        theta=relative_charge*np.pi/2-particle[f"{tracker_selected}phi"]+np.arccos(self._segment_df.radius[mask0]/(2*particle_radius))*relative_charge-np.pi/2

        mask1 = theta < 0
        theta[mask1] += 2*np.pi
        mask1 = theta < 0
        theta[mask1] += 2*np.pi
        mask1 = theta > 2*np.pi
        theta[mask1] -= 2*np.pi

        theta = np.append(theta,10*np.ones(len(mask0)-mask0.sum()))
        if not selected:
            self._segment_df.hitflag = self._segment_df.hitflag | (theta > self._segment_df.begin) & (theta < self._segment_df.end) & (self._segment_df.radius <= abs(2*particle_radius))

        return (theta > self._segment_df.begin) & (theta < self._segment_df.end) & (self._segment_df.radius <= abs(2*particle_radius))

    def get_hits_and_misses(self,particles_dataframe,particle_index): #compares hits of the true particles and hits of the simulated particles, not needed anymore
        particle=particles_dataframe.loc[particle_index,:]
        hits=np.logical_and(self.get_hits(particle),(self._segment_df.hitflag | self._segment_df.noiseflag))
        misses=np.logical_and(self.get_hits(particle),np.logical_not(self._segment_df.hitflag | self._segment_df.noiseflag))
        return [hits.sum(),misses.sum()]

    def get_hit_segments(self, particles_dataframe, particle_index: int = None, selected = True) -> Tuple[np.array, np.array]:
        hit_segments = np.empty((0,self._granularity, 2))
        colors = []
        for i in range(len(particles_dataframe)):
            hits = self.get_hits(particles_dataframe.loc[i,:], selected=selected)
            hit_segments = np.append(hit_segments, self.tracker_lines[hits], axis = 0)            
            if(i == particle_index):
                colors.extend(["blue"]*hits.sum())
            else:
                colors.extend(["yellow"]*hits.sum())
        return hit_segments, colors

    def get_arrowangle(self,particles: ParticlesManager, particle_index: int):                                          #returns angle for the arrow (angle of the outer most segment)
        last_hit_segment = (self.tracker_lines[self.get_hits(particles[particle_index], False),:,:])[-1,:,:]
        x,y = last_hit_segment.T
        phi = np.arctan2(x,y)
        return -np.mean(phi)+np.pi/2

    @property
    def tracker_colors(self) -> np.array:
        colors = np.array(["lightslategray"]*self.n_seg_total)
        colors[np.random.random(len(colors))<self._noise_ratio] = "red"
        return colors
    
    @property
    def linewidths(self) -> np.array:
        return np.array([self._linewidth]*self.n_seg_total)

    @property
    def n_seg_total(self) -> int:
        return int(self._n_layers*(self._n_segments + (self._n_layers+1)/2*self._add_segments))

    @property
    def tracker_base_collection(self) -> LineCollection:
        return LineCollection(self.tracker_lines, colors = self.tracker_colors, linewidths = self.linewidths)
    
    @property
    def n_layers(self) -> int:
        return self._n_layers
