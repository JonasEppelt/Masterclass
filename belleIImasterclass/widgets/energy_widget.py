from faulthandler import disable
import ipywidgets as widgets
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

from copy import deepcopy
from matplotlib.patches import FancyArrow,Rectangle
from matplotlib.collections import PatchCollection

#from belleIImasterclass.elements.missing_energy import mis
from belleIImasterclass.particlesmanager import ParticlesManager
from belleIImasterclass.widgets.blitmanager import BlitManager

class EnergyWidget():
    def __init__(self,particles_manager: ParticlesManager,total_energy=10.58,true_particles=False,px_py_sliders=False):
        self._total_energy=total_energy #in GeV
        self.true_particles=true_particles
        self._particles_manager = particles_manager
        self.px_py_sliders=px_py_sliders
        #für die Skalierung der Vektoren wird der größte Impuls als vergleich verwendet. Es wird der wahre Wert und nicht der Wert aus dem Tracker verwendet.
        #evtl wäre es besser den wert aus dem Tracker zu nehmen
        self.max_pt=np.amax(self._particles_manager._df.loc[:,"pt"])+0.0001
        #Liste für die Farben der Teilchen, das Darkmatterteilchen ist immer Rot
        self.colors=["green", "blue", "purple", "magenta", "darkorange", "brown", "darkkhaki", "aqua", "navy", "teal", "orchid", "peru", "lawngreen", "slateblue", "crimson"] *10
        #zusatz skalierungsfaktor (noch nicht implementiert)
        self.scale_factor=1

    def show(self):
        self.out = widgets.Output()
        with self.out:
            fig, self.ax = plt.subplots(figsize=(7,7),constrained_layout=True)
        self.ax.set_yticklabels([]) #Zahlen auf den Achsen wegmachen
        self.ax.set_xticklabels([])
        self.ax.set_ylim(-self.max_pt*1.4,self.max_pt*1.4) #Bereiche festlegen (asymetrisch wegen Energiebalken)
        self.ax.set_xlim(-self.max_pt*1.8,self.max_pt*1.4) 
        self.ax.scatter(0,0,s=320, marker='*',color="red") #Stern in der Mitte
        self.ax.plot([-self.max_pt*1.7,-self.max_pt*1.4],[self.max_pt*1.16,self.max_pt*1.16],color="black") #oberer Balken der Energie Säule
        self.ax.plot([-self.max_pt*1.7,-self.max_pt*1.4],[-self.max_pt*1.16,-self.max_pt*1.16],color="black") #unterer Balken der Energie Säule

        #dummy plots sind nötig um die Legende richtig zu machen
        dummyplots=[]
        for i in range(self._particles_manager.n_particles):
            dummyplot, = self.ax.plot([10000],[10000],label=("Teilchen "+str(i)))
            #beim label könnte man noch zusätzliche teilcheninfos reinmachen evtl
            dummyplots.append(dummyplot)
        dummyplot, = self.ax.plot([10000],[10000],label="Dark Matter Teilchen", color = "black")
        dummyplots.append(dummyplot)
        self.ax.legend(handles=dummyplots)

        self.patchartist = self.ax.add_collection(PatchCollection([])) #artist für Energiesäule UND Pfeile
        self.patchartist.set_animated(True)
        self.bm = BlitManager(fig.canvas ,self.patchartist)

        #Anzeigen für das Dark matter Teilchen
        self.dark_matter_text=widgets.Text(description = "", value = "Dark Matter Teilchen:", disabled=True)
        self.E_slider=widgets.FloatSlider( 0,min = 0, max = 8, step = 0.01, description = "Energie:")
        self.E_slider.observe(self.update, names = "value")
        self.charge_button=widgets.RadioButtons(options=['positive Ladung', 'negative Ladung',"ungeladen"])
        self.charge_button.observe(self.update, names = "value")
        self.mass_text=widgets.Text(description = "Masse:", value = "0", disabled=True)
        if self.px_py_sliders:  #px und py slider und pt und phi anzeige
            self.px_slider=widgets.FloatSlider( 0,min = -self.max_pt*1.25, max = self.max_pt*1.25, step = 0.01, description = "$p_x$:")
            self.px_slider.observe(self.update, names = "value")
            self.py_slider=widgets.FloatSlider( 0,min = -self.max_pt*1.25, max = self.max_pt*1.25, step = 0.01, description = "$p_y$:")
            self.py_slider.observe(self.update, names = "value")
            self.dm_pt_text=widgets.Text(description = "$p_t$", value = "0", disabled=True)
            self.dm_phi_text=widgets.Text(description = "$phi$", value = "0", disabled=True)
            self.dark_matter_box=widgets.VBox(children=[self.dark_matter_text,self.px_slider,self.py_slider,self.dm_pt_text,self.dm_phi_text,
                                                        self.E_slider,self.mass_text,self.charge_button])
            self.dark_matter_box.layout = widgets.Layout(border='solid 1px black',padding='5px 5px 5px 5px',margin='3px 3px 3px 3px',width = "318px")  
        else:                   #pt und phi slider und px und py anzeige
            self.pt_slider=widgets.FloatSlider( 0,min = 0, max = self.max_pt*1.25, step = 0.01, description = "$p_t$:")
            self.pt_slider.observe(self.update, names = "value")
            self.phi_slider=widgets.FloatSlider( 0,min = -np.pi, max = np.pi, step = 0.01, description = "$phi$:")
            self.phi_slider.observe(self.update, names = "value")
            self.dm_px_text=widgets.Text(description = "$p_x$", value = "0", disabled=True)
            self.dm_py_text=widgets.Text(description = "$p_y$", value = "0", disabled=True)
            self.dark_matter_box=widgets.VBox(children=[self.dark_matter_text,self.pt_slider,self.phi_slider,self.dm_px_text,self.dm_py_text,
                                                        self.E_slider,self.mass_text,self.charge_button])
            self.dark_matter_box.layout = widgets.Layout(border='solid 1px black',padding='5px 5px 5px 5px',margin='3px 3px 3px 3px',width = "318px")             

        #Anzeigen für das Gesamtsystem
        self.system_text=widgets.Text(description = "", value = "Gesamtsystem:", disabled=True)
        self.px_text=widgets.Text(description = "$\Sigma \, p_x$:", value = "0", disabled=True)
        self.py_text=widgets.Text(description = "$\Sigma \, p_y$:", value = "0", disabled=True)
        self.pt_text=widgets.Text(description = "$\Sigma \, p_t$:", value = "0", disabled=True)
        self.energy_text=widgets.Text(description = "$\Sigma$ Energie:", value = "0", disabled=True) 
        self.charge_text=widgets.Text(description = "$\Sigma$ Ladung:", value = "0", disabled=True)
        self.system_box=widgets.VBox(children=[self.system_text,self.px_text,self.py_text,self.pt_text,self.energy_text,self.charge_text])
        self.system_box.layout = widgets.Layout(border='solid 1px black',padding='5px 5px 5px 5px',margin='3px 3px 3px 3px',width = "318px")  

        #Updatebutton
        self.update_button = widgets.Button(description='Update!',disabled=False,tooltip='Update',icon='rotate-right')
        self.update_button.layout = widgets.Layout(width = "318px",margin='3px 3px 3px 3px')  
        self.update_button.on_click(self.update)
        
        self.box=widgets.VBox(children=[self.update_button,self.dark_matter_box,self.system_box])
        self.final_box = widgets.HBox(children=[self.box, self.out])
        with self.out:
            plt.show()
        display(self.final_box)
        self.update(0)

    def update(self,change):
        arrows=[]
        bars=[]
        colors=[]
        totalcharge=0
        totalpx=0
        totalpy=0
        totalpz=0
        totalenergy=0

        #Sichtbare Teilchen:
        for index in range(self._particles_manager.n_particles):
            #Ladungen:
            totalcharge+=self._particles_manager._df.loc[index,"charge"] if self.true_particles else self._particles_manager._df.loc[index,"tracker_charge"]
            #Impulse:
            pt=self._particles_manager._df.loc[index,"pt"] if self.true_particles else self._particles_manager._df.loc[index,"tracker_pt"]
            phi=self._particles_manager._df.loc[index,"phi"] if self.true_particles else self._particles_manager._df.loc[index,"tracker_phi"]  
            py=pt*np.cos(phi)
            totalpy += py
            px=pt*np.sin(phi)
            totalpx += px
            if px==0 and py==0:
                px=0.000000001
            totalpz+=self._particles_manager._df.loc[index,"pz"]                          
            #Energien:
            energy=self._particles_manager._df.loc[index,"energy"] if self.true_particles else self._particles_manager._df.loc[index,"ecl_energy"]
            totalenergy+=energy
            #Pfeile und Balken:
            bars.append(Rectangle(xy=(-self.max_pt*1.6,(totalenergy-energy)*(2.3*self.max_pt/self._total_energy)-self.max_pt*1.15),
                                  width=self.max_pt*0.1,height=energy*(2.3*self.max_pt/self._total_energy)))
            colors.append(self.colors[index])       
            arrows.append(FancyArrow(0,0,px,py,width=0.07))

        #Dark matter Teilchen:
        #Ladnung:    
        totalcharge+= 1*(self.charge_button.value=='positive Ladung')-1*(self.charge_button.value=='negative Ladung')
        #Impuls:
        if self.px_py_sliders:
            px=self.px_slider.value
            py=self.py_slider.value
        else:
            px=self.pt_slider.value*np.cos(self.phi_slider.value)
            py=self.pt_slider.value*np.sin(self.phi_slider.value)
        totalpx += px
        totalpy += py
        pz=-totalpz
        #Energie:    
        totalenergy+=self.E_slider.value
        #Masse aus Energie und Impuls:
        mass=np.sqrt(self.E_slider.value**2 - (px**2+py**2+pz**2))  if self.E_slider.value**2 >= (px**2+py**2+pz**2) else None
        #Pfeil und Balken:
        bars.append(Rectangle(xy=(-self.max_pt*1.6,(totalenergy-self.E_slider.value)*(2.3*self.max_pt/self._total_energy)-self.max_pt*1.15),
                              width=self.max_pt*0.1,height=self.E_slider.value*(2.3*self.max_pt/self._total_energy)))      
        arrows.append(FancyArrow(0,0,px,py,width=0.07))
        colors.append("black")

        #Textanzeigen:
        if self.px_py_sliders:
            self.dm_pt_text.value=str(np.round(np.sqrt(self.px_slider.value**2+self.py_slider.value**2),2))
            self.dm_phi_text.value=str(np.round(np.arctan2(self.py_slider.value,self.px_slider.value),2))
        else:
            self.dm_px_text.value=str(np.round(self.pt_slider.value*np.cos(self.phi_slider.value),2))
            self.dm_py_text.value=str(np.round(self.pt_slider.value*np.sin(self.phi_slider.value),2))
        self.pt_text.value=str(np.round(np.sqrt(totalpx**2+totalpy**2),2))
        self.px_text.value=str(np.round(totalpx,2))
        self.py_text.value=str(np.round(totalpy,2))
        self.mass_text.value=str(np.round(mass,2)) if mass is not None else "Fehler: Energie < Gesamtimpuls"
        self.energy_text.value=str(np.round(totalenergy,2))+"GeV"
        self.charge_text.value=str(totalcharge)
        #Werte merken:
        self._particles_manager.missing_particle_measurement(0,px,py,pz if np.sqrt(totalpx**2+totalpy**2)<0.1 else 0,self.E_slider.value,mass if mass is not None else 0,1*(self.charge_button.value=='positive Ladung')-1*(self.charge_button.value=='negative Ladung'))
        #Zeichnen:
        self.patchartist.set_paths(arrows+bars)
        self.patchartist.set_color(colors+colors)
        self.bm.update()


