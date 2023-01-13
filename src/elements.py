import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

class Tracks:
    '''
    Class representing particle tracks inside the tracking detector.
    
    Attributes
        - pt: float
            transverse momentum of particle in GeV
        - phi: float
            polar angle of particle in rad, between -pi and pi
        - charge: int
            electric charge of particle, either -1 or 1
        - granularity: int
            granularity of points to return for drawing

    Properties
        

    Functions
        - get_trace_array
            Returns array for Track in the detector depending on the B-field
        - track_radius: float
            radius of the particles track curve depending on the b-field
        - track_center_x: float
            x coordinate of the center point of the particles track curve
        - track_center_y: float
            y coordinate of the center point of the particles track curve
    '''
    def __init__(self, pt, phi, charge, granularity=100) -> None:
        self._pt = pt
        self._phi = phi
        self._charge = charge
        self._granularity = granularity
    
    def get_trace_array(self, B):
        track_thetas = np.linspace(-np.pi/2+self.phi,np.pi/2+self.phi, self._granularity) # angle interval to draw the half circle
        track_trace_x = abs(self.track_radius(B))*np.sin(track_thetas) + self.track_center_x
        track_trace_y = abs(self.track_radius(B))*np.cos(track_thetas) + self.track_center_y
        return np.array([track_trace_x, track_trace_y])

    def get_track_radius(self, B):
        return self._pt / B

    def track_center_x(self, B):
        return  self.track_radius(B)*np.sin((self.phi + self._charge * np.pi/2))
        
    def track_center_y(self, B):
        return self.track_radius(B)*np.cos((self.phi + self._charge * np.pi/2))
    
    @property
    def pt(self):
        return self._pt
    @pt.setter
    def pt(self, value):
        self._pt = value
    
    @property
    def phi(self):
        return self._phi
    @phi.setter
    def phi(self, value):
        if value > -np.pi and value < np.pi:
            self._phi = value
        else:
            raise Exception
    
    @property
    def charge(self):
        return self._charge
    @charge.setter
    def charge(self, value):
        if value == -1 or value == 1:
            self._charge = value
        else:
            raise Exception

class Tracker:
    '''
    Class to design a tracking detector

    Attributes:
        - n_layers: int
            number of tracking layers
        - n_segments: int
            number of tracking segments in the first layer
        - add_segments: int
            number of additional segments per new layer
        - granularity: int
            controlls the number of points used to draw the tracker. More points look better, but take longer to draw
        - noise_ration: float
            percentage of all segments with noise. Must be between 0 and 1
        - segment_df: pd.Dataframe
            Dataframe with beginning and end angles, radius and noiseflags of all tracker segments
        
    Properties
        - n_seg_total: int
            total number of segments
        
    Funktions:
        - n_segments_in_layer
    '''
    def __init__(self, n_layers, n_segments, add_segments, noise_ratio, granularity, ) -> None:
        self._n_layers = n_layers
        self._noise_ratio = noise_ratio
        self._granularity = granularity
        self._n_segments = n_segments
        self._add_segments = add_segments
        self._segment_df = pd.DataFrame(
            index = np.arange(0,self.n_seg_total),
            columns = ["begin", "end","raidus", "noiseflag"]
            )
        for l in range(1,self._n_layers):
            len_segment = 2*np.pi/(self._n_segments+k*l) # length of segment in layer l in rad
            for i in range(n_segments )

    @property
    def n_seg_total(self):
        return self._n_layers*(self._n_segments + (self._n_layers+1)/2*self._add_segments)


