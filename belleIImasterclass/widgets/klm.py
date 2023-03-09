from faulthandler import disable
import ipywidgets as widgets
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

from copy import deepcopy

from matplotlib.collections import LineCollection

from belleIImasterclass.elements.klm_detector import klm_detector
from belleIImasterclass.particlesmanager import ParticlesManager
from belleIImasterclass.widgets.blitmanager import BlitManager
from copy import deepcopy
from matplotlib.collections import LineCollection



class KLMWidget():
    def __init__(self,particles_manager: ParticlesManager,always_hit=False,true_particles=False,B=0.05):
        self.always_hit=always_hit
        self._particles_manager = particles_manager
        self.B=B
        self.klmsegments=18
        self.segmentwidth=4
        self.klmradius=19
        self.truepart=true_particles
        
        self.klm_detector=klm_detector(self._particles_manager,self.klmsegments,self.segmentwidth,self.klmradius,self.B,self.always_hit)

        self.klm_collection=self.klm_detector.make_klm_collection() 
        self.hit_collection=self.klm_detector.make_hit_collection()
        self.tracker_collection=self.klm_detector.make_tracker_collection()
        self.ecl_collection=self.klm_detector.make_ecl_collection()
    
    def update(self, change):
        self.index=self.tabs.selected_index if self.tabs.selected_index is not None else self.index

        charge=self._particles_manager._df.iloc[self.index]["tracker_charge"] if self.truepart==False else self._particles_manager._df.iloc[self.index]["charge"]
        phi_0=self._particles_manager._df.iloc[self.index]["tracker_phi"] if self.truepart==False else self._particles_manager._df.iloc[self.index]["phi"]
        R_0=self._particles_manager._df.iloc[self.index]["tracker_pt"]/self.B if self.truepart==False else self._particles_manager._df.iloc[self.index]["pt"]/self.B

        trace = self.klm_detector.make_trace(charge,phi_0,R_0)
        self.lineartist.set_segments([trace])
        self.lineartist.set_colors(["red"])
        self.bm.update()        
        self._particles_manager.Klm_measurement(self.index, (self.tickbox[self.index].value == "ja"))

    def show(self):

        self.out = widgets.Output()
        with self.out:
            fig, ax = plt.subplots(figsize=(7,7),constrained_layout=True)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_ylim(-28,28)
        ax.set_xlim(-28,28)
        self.lineartist = ax.add_collection(LineCollection([]))
        self.lineartist.set_animated(True)
        ax.add_collection(self.tracker_collection)
        ax.add_collection(self.ecl_collection)
        ax.add_collection(self.klm_collection)
        ax.add_collection(self.hit_collection)
        self.bm = BlitManager(fig.canvas ,self.lineartist)

        self.tabs = widgets.Accordion()
        self.tabs.observe(self.update, names = "selected_index")
        self.tickbox = []
        self.box_list = []
        self.boxtext=widgets.Text(value = "Wurde hier ein Teilchen erkannt?", disabled = True)
        for i in range(self._particles_manager.n_particles): 
            self.tabs.set_title(i,f"Teilchen {i}")
            self.tickbox.append(widgets.RadioButtons(options=['ja', 'nein']))
            self.tickbox[i].observe(self.update, names = "value")
            self.box_list.append(widgets.HBox([self.boxtext,self.tickbox[i]]))
        self.tabs.children = self.box_list
        self.final_box = widgets.VBox(children=[self.tabs, self.out])
        with self.out:
            plt.show()
        display(self.final_box)
        self.update(0)