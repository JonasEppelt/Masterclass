import ipywidgets as widgets
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PatchCollection
from belleIImasterclass.particlesmanager import ParticlesManager


class IDWidget:
    def __init__(self,particles_manager: ParticlesManager,cheat_mode=True, p_cheating_threshhold = 1e-2, e_cheating_threshhold=1e-2) -> None:
        self.pm=particles_manager
        self.cheat_mode=cheat_mode
        self.p_threshhold=p_cheating_threshhold
        self.e_threshhold=e_cheating_threshhold
        self.make_true_particle_data()
        self.patch_collections=[0]*self.pm.n_particles

    def get_values(self): 
        self.energies=self.pm._df.loc[:,"ecl_energy"].to_numpy()
        self.px = self.pm._df.loc[:,"tracker_pt"].to_numpy()*np.cos(self.pm._df.loc[:,"tracker_phi"].to_numpy())
        self.py = self.pm._df.loc[:,"tracker_pt"].to_numpy()*np.sin(self.pm._df.loc[:,"tracker_phi"].to_numpy())
        self.pz = self.pm._df.loc[:,"pz"].to_numpy()
        self.KLM_hits = self.pm._df.loc[:,"klm_detect"].to_numpy()
        self.charge=self.pm._df["tracker_charge"].to_numpy()  

        if(self.cheat_mode):
            self.px_cheat_mask = (abs(self.px - self.pm._df["px"].to_numpy()) < self.p_threshhold)
            self.py_cheat_mask = (abs(self.py - self.pm._df["py"].to_numpy()) < self.p_threshhold)
            self.energies_cheat_mask = (abs(self.energies - self.pm._df["energy"].to_numpy()) < self.e_threshhold)
            self.px[self.px_cheat_mask] =  self.pm._df.loc[self.px_cheat_mask,"px"].to_numpy()
            self.py[self.py_cheat_mask] =  self.pm._df.loc[self.py_cheat_mask,"py"].to_numpy()
            self.energies[self.energies_cheat_mask] = self.pm._df.loc[self.energies_cheat_mask,"energy"] 

        self.energies=self.energies[self.pm.index]
        self.px=self.px[self.pm.index]
        self.py=self.py[self.pm.index]
        self.pz=self.pz[self.pm.index]
        self.KLM_hits=self.KLM_hits[self.pm.index]
        self.charge=self.charge[self.pm.index]

        for i in range(self.pm.n_particles):
            self.patch_collections[i]=PatchCollection(self.pm._df.loc[i,"patches"],color=self.pm._df.loc[i,"patch_facecolors"])  
        self.impuls = np.sqrt(self.px**2 + self.py**2 + self.pz**2)
        self.masse = np.sqrt(abs(self.energies**2 - self.impuls**2))            


    def update(self, change = 0):
        self.get_values()
        sele_index = self.tabs.selected_index

        self.sel_charge[sele_index].value = str(self.truth_particles.loc[self.part_ids[sele_index].value, "el. Ladung"])
        self.sel_mass[sele_index].value = str(self.truth_particles.loc[self.part_ids[sele_index].value, "Masse"])
        self.sel_image[sele_index].value = self.truth_particles.loc[self.part_ids[sele_index].value, "Image"]
        self.sel_label[sele_index].value = "So sieht ein typisches "+self.part_ids[sele_index].value + " Teilchen im Ecal aus:"
        self.sel_KL0[sele_index].value = self.truth_particles.loc[self.part_ids[sele_index].value, "K_L0"] + " im KLM Detektor"
        self.sel_E_p[sele_index].value = self.truth_particles.loc[self.part_ids[sele_index].value, "E_p"]

        self.KL0_txt[sele_index].value = "Hit im KLM Detektor" if self.KLM_hits[sele_index] else "kein Hit im KLM Detektor"
        self.energy_txt[sele_index].value = str(self.energies[sele_index])
        self.charge_txt[sele_index].value = str(self.charge[sele_index])
        self.moment_txt[sele_index].value = str(self.impuls[sele_index])
        self.invmas_txt[sele_index].value = str(self.masse[sele_index])
        self.E_p_txt[sele_index].value = str(self.energies[sele_index]/self.impuls[sele_index])
        self.px_txt[sele_index].value = str(self.px[sele_index])
        self.py_txt[sele_index].value = str(self.py[sele_index])
        self.pz_txt[sele_index].value = str(self.pz[sele_index])
        self.ax.cla()
        self.ax.set_ylim(-28,28)
        self.ax.set_xlim(-28,28)
        self.ax.set_yticklabels([])
        self.ax.set_xticklabels([])
        patchartist=self.ax.add_collection(self.patch_collections[sele_index])
        patchartist.set_edgecolors(self.pm._df.loc[sele_index,"patch_edgecolors"])

    def show(self):
        boxes = []
        self.energy_txt = []
        self.px_txt = []
        self.py_txt = []
        self.pz_txt = []
        self.charge_txt = []
        self.moment_txt = []
        self.invmas_txt = []
        self.E_p_txt = []
        self.sel_mass = []
        self.sel_charge = []
        self.part_ids = []
        self.KL0_txt=[]
        self.sel_KL0=[]
        self.sel_E_p=[]
        self.sel_label=[]
        self.sel_image=[]
        self.out = widgets.Output()
        with self.out:
            self.fig, self.ax = plt.subplots(figsize=(3.3,3.3),constrained_layout=True)
        self.ax.set_ylim(-28,28)
        self.ax.set_xlim(-28,28)
        self.label1=widgets.Text(value = "Resultate", disabled = True)
        self.label2=widgets.Text(value = "Bekannte Teilchen zum Vergleichen", disabled = True)
        self.update_button = widgets.Button(description='Update!',disabled=False,tooltip='Update',icon='rotate-right')
        
        for i in range(self.pm.n_particles):
            self.px_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "$p_x$", disabled = True))
            self.py_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "$p_y$", disabled = True))
            self.pz_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "$p_z$", disabled = True))
            self.energy_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "Energie", disabled = True))
            self.charge_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "el. Ladung", disabled = True))
            self.moment_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "Impuls", disabled = True))
            self.invmas_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "Masse", disabled = True))
            self.E_p_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "$E/p$", disabled = True))
            self.KL0_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "KLM", disabled = True))   
            self.res_box = widgets.VBox(children=[widgets.HBox([self.label1,self.update_button]), self.energy_txt[i], self.charge_txt[i], self.moment_txt[i], self.invmas_txt[i], self.E_p_txt[i],self.px_txt[i],self.py_txt[i],self.pz_txt[i],self.KL0_txt[i],self.out])
            self.res_box.layout = widgets.Layout(border='solid 1px black',margin='0px 10px 10px 0px',padding='5px 5px 5px 5px',height = "705px ",width = "370px")            

            self.part_ids.append(widgets.Select(options = self.truth_particles.index, value = "e+", description = "Teilchen"))
            self.part_ids[i].observe(self.update, "value")
            self.sel_mass.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "Masse", disabled = True))
            self.sel_charge.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "el. Ladung", disabled = True))
            self.sel_E_p.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "$E/p$", disabled = True))
            self.sel_KL0.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "KLM", disabled = True))
            self.sel_label.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", disabled = True))
            self.sel_image.append(widgets.Image(value=self.truth_particles.loc["e+", "Image"],format='png',width=320,height=320))
            self.sel_box = widgets.VBox(children=[self.label2, self.part_ids[i], self.sel_mass[i], self.sel_charge[i], self.sel_KL0[i],self.sel_E_p[i],self.sel_label[i],self.sel_image[i]])
            self.sel_box.layout = widgets.Layout(border='solid 1px black',margin='0px 10px 10px 0px',padding='5px 5px 5px 5px',height = "705px ",width = "370px")  

            box = widgets.HBox(children=[self.res_box, self.sel_box])
            boxes.append(box)
        self.tabs = widgets.Tab(children=boxes)
        self.tabs.observe(self.update, "selected_index")
        self.update_button.on_click(self.update)
        for i in range(self.pm.n_particles):
            self.tabs.set_title(i,f"Teilchen {i}")
        self.update()
        with self.out:
            plt.show()
        display(self.tabs)

    def make_true_particle_data(self):

        true_particle_data = [[0.511, 1],
                      [0.511, -1],
                      [105., +1],
                     [105., -1],
                     [1776., +1],
                     [1776., -1],
                     [938.3, +1],
                     [938.3, -1],
                     [939.6, 0],
                      [135, 0],
                     [139.6, +1],
                      [139.6, -1],
                      [497.6, 0],
                     [493.7, +1],
                      [493.7, -1]]
        true_particle_names = ["e+", "e-", "mu+", "mu-", "tau+", "tau-", "Proton", "Antiproton", "Neutron", "pi0", "pi+", "pi-", "K0", "K+", "K-"]
        for n in range(len(true_particle_names)):
            file=open("Ecal_images/"+true_particle_names[n]+".png", "rb")
            true_particle_data[n].append(file.read())
            if true_particle_names[n]=="mu-" or true_particle_names[n]=="mu+":
                true_particle_data[n].append("Hit")
            else:
                true_particle_data[n].append("kein Hit")
            if true_particle_names[n]=="e+" or true_particle_names[n]=="e-":
                true_particle_data[n].append("≈1")
            else:
                true_particle_data[n].append("≠1")    
        self.truth_particles = pd.DataFrame(columns = ["Masse", "el. Ladung","Image","K_L0","E_p"], data=true_particle_data, index=true_particle_names)
        self.truth_particles.loc[:, "Masse"] = self.truth_particles["Masse"]*10**(-3)