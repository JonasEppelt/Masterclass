from faulthandler import disable
import ipywidgets as widgets
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

from copy import deepcopy

from matplotlib.collections import LineCollection


from belleIImasterclass.particlesmanager import ParticlesManager
from belleIImasterclass.widgets.blitmanager import BlitManager
from copy import deepcopy
from matplotlib.collections import LineCollection



class KLMWidget():
    def __init__(self,particles_manager: ParticlesManager,always_hit=False,B=0.1):
        self.always_hit=always_hit
        self._particles_manager = particles_manager
        self.B=B
        self.klmsegments=18
        self.segmentwidth=4
        self.klmradius=19

        for i in range(self._particles_manager.n_particles): 
            if (self._particles_manager._df.iloc[i]["pdg"]==13 or self._particles_manager._df.iloc[i]["pdg"]==-13 or self.always_hit):
                charge=self._particles_manager._df.iloc[i]["charge"]
                phi_0=self._particles_manager._df.iloc[i]["phi"]
                R_0=self._particles_manager._df.iloc[i]["pt"]/self.B
                trace = self.make_trace(charge,phi_0,R_0)
                for l in range(50,80):
                    R=np.sqrt(trace[l,0]**2+trace[l,1]**2)
                    if abs(R-self.klmradius)<0.2:
                        inner_phi=np.arctan2(trace[l,1],trace[l,0])
                    elif abs(R-self.klmradius-self.segmentwidth)<0.2:
                        outer_phi=np.arctan2(trace[l,1],trace[l,0])


        #self.tracker=Tracker(layers = 13, n_segments = 2, ecl_segments=14, k=2,dist=0.2, noise = 0, linewidth = 2, ignore_noise = True,granularity=100,trackercolor="gray")
        self.make_klm_collection()

        self.out = widgets.Output()
        with self.out:
            fig, ax = plt.subplots(figsize=(7,7),constrained_layout=True)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_ylim(-28,28)
        ax.set_xlim(-28,28)
        self.lineartist = ax.add_collection(LineCollection([]))
        self.lineartist.set_animated(True)
        #ax.add_collection(self.tracker.get_tracker_collection())

        self.ecl_collection=LineCollection([16*np.array([np.cos(np.linspace(0,6.3)),np.sin(np.linspace(0,6.3))]).T], color = np.array([1,0,0,0.6]), linewidths = 5)
        ax.add_collection(self.ecl_collection)

        ax.add_collection(self.klm_collection)
        self.bm = BlitManager(fig.canvas ,self.lineartist)


    def make_klm_collection(self):
        self.segments_coords=np.zeros((self.klmsegments,100,2))
        self.segments_angle=np.zeros((self.klmsegments,2))
        for i in range(self.klmsegments):
            points=np.zeros((100,2))
            self.segments_angle[i]=np.array([i*2*np.pi/self.klmsegments+0.015,(i+1)*2*np.pi/self.klmsegments-0.015])
            t = np.linspace(self.segments_angle[i,0],self.segments_angle[i,1], 25)   
            t_rev = np.linspace(self.segments_angle[i,1],self.segments_angle[i,0], 25) 
            points[np.arange(0,25)]=self.klmradius*np.array([np.sin(t),np.cos(t)]).T
            points[np.arange(50,75)]=(self.segmentwidth+self.klmradius)*np.array([np.sin(t_rev),np.cos(t_rev)]).T
            points[np.arange(25,50)]=np.array([np.linspace(points[24,0],points[50,0],25),np.linspace(points[24,1],points[50,1],25)]).T
            points[np.arange(75,100)]=np.array([np.linspace(points[74,0],points[0,0],25),np.linspace(points[74,1],points[0,1],25)]).T
            self.segments_coords[i]=points
        self.klm_collection=LineCollection(self.segments_coords, color = np.array([0,0,1,0.9]), linewidths = 3)

    def make_trace(self,charge,phi_0,R_0):
        magnetradius=17.5
        if 2*R_0 > magnetradius:
            r=np.linspace(0,magnetradius,50)
            theta=-phi_0+np.arccos(r/(2*R_0))*charge+(charge-1)*np.pi/2
            phi_2=-np.arctan2(r[48]*np.cos(theta[48])-r[1+48]*np.cos(theta[1+48]),r[48]*np.sin(theta[48])-r[1+48]*np.sin(theta[1+48]))
            theta2=phi_2-np.arccos(r/(2*R_0))*charge+np.pi*(charge/2-0.5)        
            trace=np.append(np.array([r*np.cos(theta),r*np.sin(theta)]).T,
                            np.array([r[-1]*np.cos(theta[-1])+r*np.cos(theta2),r[-1]*np.sin(theta[-1])+r*np.sin(theta2)]).T,axis=0)
        elif R_0==0:
            trace=np.zeros((100,2))
        else:
            r=np.linspace(0,2*R_0,100)
            theta=-phi_0+np.arccos(r/(2*R_0))*charge+(charge-1)*np.pi/2  
            trace=np.array([r*np.cos(theta),r*np.sin(theta)]).T     
        return trace

    def update(self, change):
        self.index=self.tabs.selected_index if self.tabs.selected_index is not None else self.index

        charge=self._particles_manager._df.iloc[self.index]["tracker_charge"]
        phi_0=self._particles_manager._df.iloc[self.index]["tracker_phi"]
        R_0=self._particles_manager._df.iloc[self.index]["tracker_pt"]/self.B

        trace = self.make_trace(charge,phi_0,R_0)
        self.lineartist.set_segments([trace])
        self.lineartist.set_colors(["red"])
        self.bm.update()        

    def show(self):

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

    @property
    def KLM_hit(self):
        hit = []
        for i in range(len(self.tickbox)):
            hit.append(1 if (self.tickbox[i].value == "ja") else 0)
        return hit