import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
import ipywidgets as widgets
from ipywidgets import *
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
    '''
    def __init__(self, 
    n_layers: int,
    n_segments: int, 
    add_segments: int, 
    noise_ratio: float, 
    granularity: int, 
    linewidth: float = 8,
    ) -> None:
        self._n_layers = n_layers
        self._noise_ratio = noise_ratio
        self._granularity = granularity
        self._n_segments = n_segments
        self._add_segments = add_segments
        self._linewidth = linewidth
        self._segment_df = pd.DataFrame(
            index = np.arange(0,self.n_seg_total, dtype = int),
            columns = ["begin", "end","radius", "noiseflag"]
            )
        
        segments_counter = 0
        for l in range(1,self._n_layers+1):
            len_segment = 2*np.pi/(self._n_segments+l) # length of segment in layer l in rad
            for i in range(n_segments+l*self._add_segments):
                self._segment_df.loc[segments_counter, "begin"] = len_segment*i+0.2/(l+1)
                self._segment_df.loc[segments_counter, "end"] = len_segment*(i+1)-0.2/(l+1)
                self._segment_df.loc[segments_counter, "radius"] = l
                segments_counter += 1

    @property
    def tracker_lines(self) -> np.array:
        '''
        returns for each segments lines in the given granularity, in kartesian coordinates
        '''
        theta = np.linspace(self._segment_df["begin"], self._segment_df["end"], self._granularity, dtype = float)
        return self._segment_df["radius"].to_numpy()[None,None,:] * np.array([np.cos(theta),np.sin(theta)])

    @property
    def tracker_colors(self) -> np.array:
        return np.array(["black"]*self.n_seg_total)
    
    @property
    def linewidths(self) -> np.array:
        return np.array([self._linewidth]*self.n_seg_total)

    @property
    def n_seg_total(self) -> int:
        return self._n_layers*(self._n_segments + (self._n_layers+1)/2*self._add_segments)

    @property
    def tracker_base_collection(self) -> LineCollection:
        return LineCollection(self.tracker_lines, colors = self.tracker_colors, linewidths = self.linewidths)


class BlitManager: 
    '''
    managing the blitting for the interactive widgets: only changed parts will be drawn in each iteration
    '''
    def __init__(self, canvas, artist, artist2=None):
        """copy from matplotlib website (blitting)"""
        self._canvas = canvas
        self._background = None
        self._artist = artist
        if artist2 is not None:
            self._artist2=artist2
            self._twoartists=True
        else:
            self._twoartists=False
        
        # grab the background on every draw
        self.cid = canvas.mpl_connect("draw_event", self.on_draw)

    def on_draw(self, event) -> None:
        canvas = self._canvas
        if event is not None:
            if event.canvas != canvas:
                raise RuntimeError
        self._background = canvas.copy_from_bbox(canvas.figure.bbox)
        self._draw_animated()


    def _draw_animated(self) -> None:
        fig = self._canvas.figure
        fig.draw_artist(self.artist)
        if self.twoartists:
            fig.draw_artist(self.artist2)

    def update(self) -> None:
        cv = self.canvas
        fig = cv.figure
        if self._background is None:
            self.on_draw(None)
        else:
            cv.restore_region(self._bg)
            self._draw_animated()
            cv.blit(fig.bbox)
        cv.flush_events()

class TrackingWidget:
    '''
    Widget displaying the tracker and manipulating it
    '''
    def __init__(self, particle_df: pd.DataFrame, n_layers: int, n_segments: int, noise_rate: float = 0, continuous_update = True) -> None:
        self._particle_df = particle_df
        self._widget_df = pd.DataFrame(index = particle_df.index)
        self._continuous_update = continuous_update
        self._widget_df["hits_counter_widget"] = self.generate_widget_per_particle(Label, value="0 hits & 0 misses")
        self._widget_df["pt_slider_widget"] = self.generate_widget_per_particle(FloatSlider, 
            value = 0, min = 0,  max = 5, step = 0.01, description = "$p_T$", continuous_update = self._continuous_update)
        self._widget_df["pt_fineslider_widget"] = self.generate_widget_per_particle(FloatSlider,
            value = 0, min = 0, max = 1, step = 0.001, description = "$p_T$, fein", continuous_update = self._continuous_update)
        self._widget_df["phi_slider_widget"] = self.generate_widget_per_particle(FloatSlider,
            value = 0, min = -np.pi, max = np.pi, step = 0.01, description = "$\phi$", continuous_update = self._continuous_update)
        self._widget_df["phi_fineslider_widget"] = self.generate_widget_per_particle(FloatSlider,
            value = 0, min = 0, max = 0.15, step = 0.001, description = "$\phi$ fein", continuous_update = self.continuous_update)]
        self._widget_df["charge_widget"] = self.generate_widget_per_particle(RadioButtons,
            options = ["positiv", "negativ"], description = "elektrische Ladung:")

        for widget_column_name in ["pt_slider_widget", "pt_fineslider_widget", "phi_slider_widget", "phi_fineslider_widget", "charge_widget"]:
            self._widget_df[widget_column_name].apply(lambda x: x.observe(self.update, names = "value"))

        self._widget_df["widget_box"] = self._widget_df.apply(lambda x:
            VBox([x["hits_counter_widget"], x["pt_slider_widget"], x["pt_fineslider_widget"], 
                x["phi_slider_widget"], x["phi_fineslider_widget"], x["charge_widget"]]))
        
        self.particle_selector = Accordion(children=self._widget_df["widget_box"].to_list(), titles = list(range(1,self.n_particles+1)))
        self.particle_selector.observe(self.change_particle, names = "selected_index")
        
        particle_box = HBox(children=[self.particle_selector])
        particle_box.layout = Layout(
                                border='solid 1px black',
                                margin='0px 10px 10px 0px',
                                padding='5px 5px 5px 5px',
                                height = "750px ",
                                width = "500px"
                            )
        plot_box = HBox([self.out])
        plot_box.layout = Layout(
                                border='solid 1px black',
                                margin='0px 10px 10px 0px',
                                padding='5px 5px 5px 5px',
                                height = "750px ",
                                width = "750px"
                            )
        self.final_box = HBox([  self.tabs_box,self.plot_box])



    def generate_widget_per_particle(self, widget, kwargs) -> pd.Series:
        return self._widget_df.apply(lambda x: widget(**kwargs))

    @property
    def n_particles(self) -> int:
        return len(self._particle_df)

    def show(self):
        self._out = widgets.Output()
        with self.out:
            fig, ax = plt.subplots(figsize=(7,7), constrained_layout = True)
        limit = self._tracker.layers +3
        ax.set_xlim(-limit,limit)
        ax.set_ylim(-limit,limit)
        ax.set_axis_off()
        artist = ax.add_collection(LineCollection([]))
        artist.set_animated(True)
        fig.canvas.toolbar_position = "left"
        ax.add_collection(self._tracker.tracker_base_collection)
        blitmanager = BlitManager(fig.canvas, artist)
        with self._out:
            plt.show()
        fig.canvas.draw()
        display(self.final_box)
        self.update(True)