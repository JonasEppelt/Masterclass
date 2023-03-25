from faulthandler import disable
import ipywidgets as widgets
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

from copy import deepcopy
from matplotlib.patches import FancyArrow
from matplotlib.collections import PatchCollection

#from belleIImasterclass.elements.missing_energy import mis
from belleIImasterclass.particlesmanager import ParticlesManager
from belleIImasterclass.widgets.blitmanager import BlitManager

class EnergyWidget():
    def __init__(self,particles_manager: ParticlesManager,total_energy=10.45,true_particles=False):
        self._total_energy=total_energy #in GeV
        self.true_particles=true_particles
        self._particles_manager = particles_manager
        self.max_pt=np.amax(self._particles_manager._df.loc[:,"pt"])+0.01# if self.true_particles else self._particles_manager._df.loc[:,"tracker_pt"])+0.01

    def show(self):
        self.out = widgets.Output()
        with self.out:
            fig, ax = plt.subplots(figsize=(7,7),constrained_layout=True)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_ylim(-self.max_pt*1.5,self.max_pt*1.5)
        ax.set_xlim(-self.max_pt*1.5,self.max_pt*1.5)
        ax.scatter(0,0,s=320, marker='*',color="red")
        self.patchartist = ax.add_collection(PatchCollection([]))
        self.patchartist.set_animated(True)
        self.bm = BlitManager(fig.canvas ,self.patchartist)

        self.px_slider=widgets.FloatSlider( 0,min = -self.max_pt*1.3, max = self.max_pt*1.3, step = 0.01, description = "$p_x$")
        self.px_slider.observe(self.update, names = "value")
        self.py_slider=widgets.FloatSlider( 0,min = -self.max_pt*1.3, max = self.max_pt*1.3, step = 0.01, description = "$p_y$")
        self.py_slider.observe(self.update, names = "value")
        self.E_slider=widgets.FloatSlider( 0,min = 0, max = 10, step = 0.01, description = "$Energie$")
        self.E_slider.observe(self.update, names = "value")
        self.charge_button=widgets.RadioButtons(options=['positive Ladung', 'negative Ladung',"ungeladen"])
        self.charge_button.observe(self.update, names = "value")

        self.px_text=widgets.Text(description = "px:", value = "0", disabled=True)
        self.py_text=widgets.Text(description = "py:", value = "0", disabled=True)
        self.pt_text=widgets.Text(description = "pt:", value = "0", disabled=True)
        self.energy_text=widgets.Text(description = "Energie:", value = "0", disabled=True)
        self.mass_text=widgets.Text(description = "Masse:", value = "0", disabled=True)
        self.missing_energy_text=widgets.Text(description = "fehlende Energie:", value = "0", disabled=True)
        self.charge_text=widgets.Text(description = "Ladung:", value = "0", disabled=True)
        self.update_button = widgets.Button(description='Update!',disabled=False,tooltip='Update',icon='rotate-right')
        self.update_button.on_click(self.update)
        
        self.box=widgets.VBox(children=[self.px_slider,self.py_slider,self.E_slider,self.charge_button,self.px_text,self.py_text,self.pt_text,self.mass_text,self.energy_text,self.missing_energy_text,self.charge_text,self.update_button])
        self.final_box = widgets.HBox(children=[self.box, self.out])
        with self.out:
            plt.show()
        display(self.final_box)
        self.update(0)

    def update(self,change):
        arrows=[]
        colors=[]
        totalcharge=0
        totalpx=0
        totalpy=0
        totalpz=0
        totalenergy=0
        for index in range(self._particles_manager.n_particles):
            totalcharge+=self._particles_manager._df.loc[index,"charge"] if self.true_particles else self._particles_manager._df.loc[index,"tracker_charge"]
            pt=self._particles_manager._df.loc[index,"pt"] if self.true_particles else self._particles_manager._df.loc[index,"tracker_pt"]
            totalpz+=self._particles_manager._df.loc[index,"pz"]
            phi=self._particles_manager._df.loc[index,"phi"] if self.true_particles else self._particles_manager._df.loc[index,"tracker_phi"]
            totalenergy+=self._particles_manager._df.loc[index,"energy"] if self.true_particles else self._particles_manager._df.loc[index,"ecl_energy"]
            py=pt*np.cos(phi)
            totalpy += py
            px=pt*np.sin(phi)
            totalpx += px
            if px==0:
                px=0.000001
            if py==0:
                py=0.000001            
            arrows.append(FancyArrow(0,0,px,py,width=0.05))
            colors.append("blue")
        
        px=self.px_slider.value
        totalpx += px
        py=self.py_slider.value
        totalpy += py
        pz=-totalpz
        totalcharge+= 1*(self.charge_button.value=='positive Ladung')-1*(self.charge_button.value=='negative Ladung')
        totalenergy+=self.E_slider.value
        mass=np.sqrt(abs(self.E_slider.value**2 - (px**2+py**2+pz**2)))  
        self._particles_manager.missing_particle_measurement(px,py,pz,self.E_slider.value,mass,1*(self.charge_button.value=='positive Ladung')-1*(self.charge_button.value=='negative Ladung'))          
        arrows.append(FancyArrow(0,0,px,py,width=0.05))
        colors.append("red")    
        self.pt_text.value=str(np.sqrt(totalpx**2+totalpy**2))
        self.px_text.value=str(totalpx)
        self.py_text.value=str(totalpy)
        self.mass_text.value=str(mass)
        self.energy_text.value=str(totalenergy)+"GeV"
        self.missing_energy_text.value=str(self._total_energy-totalenergy)+"GeV"
        self.charge_text.value=str(totalcharge)
        self.patchartist.set_paths(arrows)
        self.patchartist.set_color(colors)
        
        self.bm.update()


