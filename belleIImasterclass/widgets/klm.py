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
    def __init__(self,particles_manager: ParticlesManager,always_hit=False,true_particles=False,B=0.05):
        self.always_hit=always_hit
        self._particles_manager = particles_manager
        self.B=B
        self.klmsegments=18
        self.segmentwidth=4
        self.klmradius=19
        self.truepart=true_particles
        self.index = None

        self.make_klm_collection()
        self.make_hit_collection()
        self.make_tracker_ecl_collection()
        
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
    
    def make_tracker_ecl_collection(self):
        self.ecl_collection=LineCollection([16*np.array([np.cos(np.linspace(0,6.3)),np.sin(np.linspace(0,6.3))]).T], color = np.array([1,0,0,0.6]), linewidths = 5)
        lines=[]
        for l in range(1,14+1):
            len_segment = 2*np.pi/(3+2*l) # length of segment in layer l in rad
            offset = l%2
            for i in range(14+l*2):
                phi=np.linspace(len_segment*i+0.1/(l+1)+offset,len_segment*(i+1)-0.1/(l+1)+offset)
                lines.append(l*np.array([np.cos(phi),np.sin(phi)]).T)
        self.tracker_collection=LineCollection(lines, color = "gray", linewidths = 6)

    def make_klm_collection(self):
        self.segments_coords=np.zeros((self.klmsegments,100,2))
        self.segments_angle=np.zeros((self.klmsegments,2))
        for i in range(self.klmsegments):
            points=np.zeros((100,2))
            self.segments_angle[i]=np.array([i*2*np.pi/self.klmsegments+0.015,(i+1)*2*np.pi/self.klmsegments-0.015])
            t = np.linspace(self.segments_angle[i,0],self.segments_angle[i,1], 25)   
            t_rev = np.linspace(self.segments_angle[i,1],self.segments_angle[i,0], 25) 
            points[np.arange(0,25)]=self.klmradius*np.array([np.cos(t),np.sin(t)]).T
            points[np.arange(50,75)]=(self.segmentwidth+self.klmradius)*np.array([np.cos(t_rev),np.sin(t_rev)]).T
            points[np.arange(25,50)]=np.array([np.linspace(points[24,0],points[50,0],25),np.linspace(points[24,1],points[50,1],25)]).T
            points[np.arange(75,100)]=np.array([np.linspace(points[74,0],points[0,0],25),np.linspace(points[74,1],points[0,1],25)]).T
            self.segments_coords[i]=points
        self.klm_collection=LineCollection(self.segments_coords, color = np.array([0,0,1,0.9]), linewidths = 3.6)

    def make_hit_collection(self):
        self.hits=np.full((self.klmsegments), False)
        for i in range(self._particles_manager.n_particles): 
            if (self._particles_manager._df.iloc[i]["pdg"]==13 or self._particles_manager._df.iloc[i]["pdg"]==-13 or self.always_hit):
                charge=self._particles_manager._df.iloc[i]["charge"]
                phi_0=self._particles_manager._df.iloc[i]["phi"]
                R_0=self._particles_manager._df.iloc[i]["pt"]/self.B
                trace = self.make_trace(charge,phi_0,R_0)
                inner_phi=0
                outer_phi=0
                for l in range(0,100):
                    R=np.sqrt(trace[l,0]**2+trace[l,1]**2)
                    if abs(R-self.klmradius)<1:
                        inner_phi=np.arctan2(trace[l,1],trace[l,0])
                    elif abs(R-self.klmradius-self.segmentwidth)<1:
                        outer_phi=np.arctan2(trace[l,1],trace[l,0])
                    if inner_phi<0: 
                        inner_phi=inner_phi+2*np.pi
                    if outer_phi<0: 
                        outer_phi=outer_phi+2*np.pi
                self.hits=np.logical_or(self.hits,(np.logical_and(self.segments_angle[:,0]<inner_phi,self.segments_angle[:,1]>inner_phi)))
                self.hits=np.logical_or(self.hits,(np.logical_and(self.segments_angle[:,0]<outer_phi,self.segments_angle[:,1]>outer_phi)))
        self.hit_collection=LineCollection(self.segments_coords[self.hits], color = np.array([1,0,0,0.9]), linewidths = 2)

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
        if not self.index is None:
            charge=self._particles_manager._df.iloc[self.index]["tracker_charge"] if self.truepart==False else self._particles_manager._df.iloc[self.index]["charge"]
            phi_0=self._particles_manager._df.iloc[self.index]["tracker_phi"] if self.truepart==False else self._particles_manager._df.iloc[self.index]["phi"]
            R_0=self._particles_manager._df.iloc[self.index]["tracker_pt"]/self.B if self.truepart==False else self._particles_manager._df.iloc[self.index]["pt"]/self.B

            trace = self.make_trace(charge,phi_0,R_0)
            self.lineartist.set_segments([trace])
            self.lineartist.set_colors(["red"])
            self.bm.update()        
            self._particles_manager.Klong_measurement(self.index, (self.tickbox[self.index].value == "ja"))

    def show(self):
        self.tickbox = []
        self.box_list = []
        self.boxtext=widgets.Text(value = "Wurde hier ein Teilchen erkannt?", disabled = True)
        for i in range(self._particles_manager.n_particles): 
            self.tickbox.append(widgets.RadioButtons(options=['ja', 'nein']))
            self.tickbox[i].observe(self.update, names = "value")
            self.box_list.append(widgets.HBox([self.boxtext,self.tickbox[i]]))
        self.tabs = widgets.Accordion(
            titles = [f"Teilchen {str(i)}" for i in range(self._particles_manager.n_particles)],
            children=self.box_list,)
        self.tabs.observe(self.update, names = "selected_index")
        self.final_box = widgets.VBox(children=[self.tabs, self.out])
        with self.out:
            plt.show()
        display(self.final_box)
        self.update(0)