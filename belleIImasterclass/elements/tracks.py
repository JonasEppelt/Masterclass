import numpy as np

class Tracks:
    '''
    Class representing particle tracks inside the tracking detector.
    '''
    def __init__(self, pt, phi, charge, B, granularity=100) -> None:
        self._pt = pt
        self._phi = phi
        self._charge = charge
        self._B = B
        self._granularity = granularity
    
    def get_trace_array(self):
        track_thetas = np.linspace(-np.pi/2+self.phi,np.pi/2+self.phi, self._granularity) # angle interval to draw the half circle
        track_trace_x = abs(self.track_radius)*np.sin(track_thetas) + self.track_center_x
        track_trace_y = abs(self.track_radius)*np.cos(track_thetas) + self.track_center_y
        return np.array([track_trace_x, track_trace_y])
    
        
    @property
    def track_radius(self):
        return self._pt / self._B

    @property
    def track_center_x(self):
        return  self.track_radius*np.sin((self.phi + self._charge * np.pi/2))
    @property
    def track_center_y(self):
        return self.track_radius*np.cos((self.phi + self._charge * np.pi/2))

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