from ipywidgets import Label, FloatSlider, RadioButtons, Accordion, Output, Layout, HBox, VBox, Text
import pandas as pd
from belleIImasterclass.particlesmanager import ParticlesManager
from belleIImasterclass.widgets.blitmanager import BlitManager
from belleIImasterclass.elements.tracker import Tracker
from belleIImasterclass.elements.tracks import Tracks
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

import time 

class TrackingWidget:
    '''
    Widget displaying the tracker and manipulating it
    '''
    def __init__(self, particles_manager: ParticlesManager, continuous_update = True, granularity = 100,truthvalues=False,show_hitcounter=False, **kwargs) -> None:
        self._helper = None
        self._particles_manager = particles_manager
        self._granularity = granularity
        self._tracker = Tracker(**kwargs, granularity=granularity)
        #self._arrows = []
        #for i in range(len(self._particles_manager)):
        #    self._arrows.append(self.get_arrow(i))

        self._widget_df = pd.DataFrame(index = np.arange(particles_manager.n_particles))
        self._continuous_update = continuous_update
        self._widget_df["hits_counter_widget"] = self.generate_widget_per_particle(Text,
            disabled = True, description = "hit counter:",placeholder = "")
        self._widget_df["pt_slider_widget"] = self.generate_widget_per_particle(FloatSlider, 
            value = 0, min = 0,  max = 7, step = 0.01, description = "$p_T$", continuous_update = self._continuous_update)
        self._widget_df["pt_fineslider_widget"] = self.generate_widget_per_particle(FloatSlider,
            value = 0, min = 0, max = 0.5, step = 0.001, description = "$p_T$, fein", continuous_update = self._continuous_update)
        self._widget_df["phi_slider_widget"] = self.generate_widget_per_particle(FloatSlider,
            value = 0, min = -np.pi, max = np.pi, step = 0.01, description = "$\phi$", continuous_update = self._continuous_update)
        self._widget_df["phi_fineslider_widget"] = self.generate_widget_per_particle(FloatSlider,
            value = 0, min = -0.10, max = 0.10, step = 0.001, description = "$\phi$, fein", continuous_update = self._continuous_update)
        self._widget_df["charge_widget"] = self.generate_widget_per_particle(RadioButtons,
            options = ["positiv", "negativ","neutral"], description = "elektrische Ladung:")

        if truthvalues:
            for i in range(len(self._particles_manager)):
                self._widget_df.loc[i, "pt_slider_widget"].value = self._particles_manager._df.loc[i,"pt"]
                self._widget_df.loc[i, "phi_slider_widget"].value = self._particles_manager._df.loc[i,"phi"]
                charge="positiv" *(self._particles_manager._df.loc[i,"charge"]==1) +"negativ"*(self._particles_manager._df.loc[i,"charge"]==-1)+"neutral"*(self._particles_manager._df.loc[i,"charge"]==0)
                self._widget_df.loc[i, "charge_widget"].value = charge
        else:
            for i in range(len(self._particles_manager)):
                self._widget_df.loc[i, "pt_slider_widget"].value = self._particles_manager._df.loc[i,"tracker_pt"]
                self._widget_df.loc[i, "phi_slider_widget"].value = self._particles_manager._df.loc[i,"tracker_phi"]
                charge="positiv" *(self._particles_manager._df.loc[i,"tracker_charge"]==1) +"negativ"*(self._particles_manager._df.loc[i,"tracker_charge"]==-1)+"neutral"*(self._particles_manager._df.loc[i,"tracker_charge"]==0)
                self._widget_df.loc[i, "charge_widget"].value = charge

        for widget_column_name in ["pt_slider_widget", "pt_fineslider_widget", "phi_slider_widget", "phi_fineslider_widget", "charge_widget"]:
            self._widget_df.apply(lambda x: x[widget_column_name].observe(self.update, names = "value"), 1)

        if show_hitcounter:
            self._widget_df["widget_box"] = self._widget_df.apply(lambda x: 
                VBox([x["hits_counter_widget"], x["pt_slider_widget"], x["pt_fineslider_widget"], 
                    x["phi_slider_widget"], x["phi_fineslider_widget"], x["charge_widget"]]),1)
        else:
             self._widget_df["widget_box"] = self._widget_df.apply(lambda x: 
                VBox([x["pt_slider_widget"], x["pt_fineslider_widget"], 
                    x["phi_slider_widget"], x["phi_fineslider_widget"], x["charge_widget"]]),1)           
        
        self.particle_selector = Accordion(children=self._widget_df["widget_box"].to_list(), titles =[f"Teilchen {str(i)}" for i in list(range(self.n_particles))])
        self.particle_selector.observe(self.change_particle, names = "selected_index")
        self._current_particle_index = 0
        
        self._out = Output()
        
        particle_box = HBox(children=[self.particle_selector])
        particle_box.layout = Layout(
                                border='solid 1px black',
                                margin='3px 3px 3px 3px',
                                padding='5px 5px 5px 5px',
                                height = "750px ",
                                width = "355px"
                            )
        plot_box = HBox([self._out])
        plot_box.layout = Layout(
                                border='solid 1px black',
                                margin='3px 3px 3px 3px',
                                padding='5px 5px 5px 5px',
                                height = "750px ",
                                width = "750px"
                            )
        self.final_box = HBox([particle_box, plot_box])

    def generate_widget_per_particle(self, widget, **kwargs) -> pd.Series:
        return self._widget_df.apply(lambda x: widget(**kwargs),1)

    @property
    def n_particles(self) -> int:
        return self._particles_manager.n_particles

    def show(self):        
        with self._out:
            self._fig, ax = plt.subplots(figsize=(7,7), constrained_layout = True)
        limit = self._tracker.n_layers +3
        ax.set_xlim(-limit,limit)
        ax.set_ylim(-limit,limit)
        ax.set_axis_off()
        selection_hit_collection = ax.add_collection(LineCollection([]))
        self._selection_hit_collection = selection_hit_collection
        self._selection_hit_collection.set_animated(True)
        self._fig.canvas.toolbar_position = "left"
        ax.add_collection(self._tracker.tracker_base_collection)
        for i in range(len(self._particles_manager)):
            segments, colors = self._tracker.get_hit_segments(self._particles_manager._df, i, selected = False)
            colors = ["red"]*len(colors)
            ax.add_collection(LineCollection(segments, color = colors, linewidth = 2.5))
        self._selection_blitmanager = BlitManager(self._fig.canvas, self._selection_hit_collection)
        with self._out:
            plt.show()
        self._fig.canvas.draw()
        display(self.final_box)
        self.change_particle(1)
        self.update(True)
    
    def update(self, update):
        if self._current_particle_index != None:
            pt = self._widget_df.loc[self._current_particle_index, ["pt_slider_widget", "pt_fineslider_widget"]].apply(lambda x: x.value).sum()
            phi = np.clip(self._widget_df.loc[self._current_particle_index, ["phi_slider_widget", "phi_fineslider_widget"]].apply(lambda x: x.value).sum(),-np.pi,np.pi)
            charge = (self._widget_df.loc[self._current_particle_index, "charge_widget"].value == "positiv") + (-1)*(self._widget_df.loc[self._current_particle_index, "charge_widget"].value == "negativ")
            self._particles_manager.tracker_measurement(index = self._current_particle_index, pt = pt if charge!=0 else 0, phi = phi , charge = charge)
        df = self._particles_manager._df
        segments, colors = self._tracker.get_hit_segments(df, self._current_particle_index)
        if self._current_particle_index != None:
            trace = Tracks(pt = pt, phi = phi, charge = charge, B = self._tracker._B_field, granularity = self._granularity)
            trace_segments = trace.get_trace_array()
            segments = np.append(segments, [trace_segments.T], axis = 0)
            colors = np.append(colors, ["blue"])
            arrow_segments = self.get_arrow(self._current_particle_index)
            segments = np.append(segments, [arrow_segments], axis = 0)
            colors = np.append(colors, ["green"])
            hits,misses=self._tracker.get_hits_and_misses(df,self._current_particle_index)
            self._widget_df.loc[self._current_particle_index,"hits_counter_widget"].value=str(hits)+" hits & "+str(misses)+" misses"
        self._selection_hit_collection.set_segments(segments)
        self._selection_hit_collection.set_colors(colors)
        self._selection_blitmanager.update()

    def change_particle(self, update):
        self._current_particle_index = self.particle_selector.selected_index
        self.update(1)
    
    def get_arrow(self, i) -> np.array:
        if self._particles_manager._df.loc[i,"charge"]!=0:
            size = 0.18
            phi = self._tracker.get_arrowangle(self._particles_manager, i)
            posx = np.cos(phi)*(self._tracker.n_layers+2.8)
            posy = np.sin(phi)*(self._tracker.n_layers+2.8)
            x=(np.array([0,0,-2,2,0,0])*np.cos(phi+np.pi/2)-np.array([-3,3,0,0,3,-3])*np.sin(phi+np.pi/2))*size
            y=(np.array([0,0,-2,2,0,0])*np.sin(phi+np.pi/2)+np.array([-3,3,0,0,3,-3])*np.cos(phi+np.pi/2))*size
            arrow = np.append(np.array([x,y]).T,np.zeros((self._granularity-6,2)),axis=0)
            arrow = arrow + [posx,posy]
        else:
            arrow=1000*np.ones((self._granularity,2))
        return arrow
