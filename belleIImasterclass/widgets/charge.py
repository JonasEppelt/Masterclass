from belleIImasterclass.elements.B_field import B_field
from belleIImasterclass.elements.tracker import Tracker
from belleIImasterclass.elements.tracks import Tracks
from belleIImasterclass.widgets.blitmanager import BlitManager
from ipywidgets import FloatSlider, RadioButtons, Output, Layout, HBox, VBox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import pandas as pd

class ChargeWidget:
    def __init__(self, granularity = 100, continuous_update = True):
        self._granularity = granularity
        self._tracker = Tracker(granularity=granularity, noise_ratio=0)
        self._lim = self._tracker.n_layers +3
        self._b_field_vis = B_field(0, lim = self._lim)
        self._continuous_update = continuous_update
        self.pt_slider = FloatSlider(
            value=3, 
            min = 0, 
            max = 5, 
            step = 0.01, 
            description = "Impuls", 
            continuous_update = self._continuous_update,
        )
        self.pt_slider.observe(self.update, names = "value") 
        self.phi_slider = FloatSlider(
            value = 0,
            min=-np.pi,
            max=np.pi,
            step = 0.01,
            description="Winkel $\phi$",
        )
        self.phi_slider.observe(self.update, names = "value")
        self.charge_buttons = RadioButtons(
            options = ["positiv", "negativ", "neutral"],
            description="elektrisch Ladung",
        )
        self.charge_buttons.observe(self.update, names = "value")

        self.B_slider = FloatSlider(
            value = 0,
            min = -1,
            max = 1,
            step = 0.01,
            description = "Magnetfeldstärke",
            continuous_update = self._continuous_update,
        )
        self.B_slider.observe(self.update, names = "value")
        self._out = Output()
        particle_box = VBox(children=[self.pt_slider, self.phi_slider, self.charge_buttons, self.B_slider])
        particle_box.layout = Layout(
                                border='solid 1px black',
                                margin='0px 10px 10px 0px',
                                padding='5px 5px 5px 5px',
                                height = "750px ",
                                width = "500px"
                            )
        plot_box = HBox(children=[self._out])
        plot_box.layout = Layout(
                                border='solid 1px black',
                                margin='0px 10px 10px 0px',
                                padding='5px 5px 5px 5px',
                                height = "750px ",
                                width = "750px"
                            )
        self.final_box = HBox(children=[particle_box, plot_box])

    def show(self):
        with self._out:
            self._fig, self._ax = plt.subplots(figsize = (7,7), constrained_layout=True)
        self._ax.set_xlim(-self._lim, self._lim)
        self._ax.set_ylim(-self._lim, self._lim)
        self._ax.set_axis_off()

        self.tracker_col = self._ax.add_collection(self._tracker.tracker_base_collection)


        trace_collection = self._ax.add_collection(LineCollection([]))
        self._trace_collection = trace_collection
        self._trace_collection.set_animated(True)
        self._fig.canvas.toolbar_position = "left"
        self._trace_blitmanager = BlitManager(self._fig.canvas, self._trace_collection,)

        self._field_collection = self._b_field_vis.plot_field(self._ax)

        with self._out:
            plt.show()
        self._fig.canvas.draw()
        display(self.final_box)
        self.update()
    
    def update(self, change = None):
        pt = self.pt_slider.value
        phi = self.phi_slider.value
        charge_str = self.charge_buttons.value
        B = self.B_slider.value
        #B = B if B!=0 else 1e-9
        if charge_str == "negativ":
            charge = -1
        elif charge_str == "positiv":
            charge = 1
        else:
            charge = 0
        self._b_field_vis._B_field = B
        
        self._ax.clear()
        self._ax.set_xlim(-self._lim, self._lim)
        self._ax.set_ylim(-self._lim, self._lim)
        self._ax.set_axis_off()
        self._field_collection = self._b_field_vis.plot_field(self._ax)

        trace = Tracks(pt=pt, phi = phi, charge=charge, B = B, granularity=self._granularity)
        segments = np.array([trace.get_trace_array().T])
        colors = ["blue"]

        
        self._trace_collection.set_segments(segments)
        self._trace_collection.set_colors(colors)
        self._trace_blitmanager.update()
        