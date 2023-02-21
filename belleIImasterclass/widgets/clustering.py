from matplotlib.path import Path
from matplotlib.widgets import LassoSelector
from belleIImasterclass.particlesmanager import ParticlesManager
from belleIImasterclass.elements.ecal import ECal
from belleIImasterclass.widgets.blitmanager import BlitManager
from ipywidgets import Output, Accordion, Text, HBox, VBox
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
from matplotlib.collections import PatchCollection, LineCollection

class ECLWidget:
    def __init__(self,particles_manager: ParticlesManager, noise_rate = 0.05) -> None:
        self._particles_manager = particles_manager
        self._noise_rate = noise_rate
        self._ecal = ECal()        
        self._selected_crystalls = [np.array([], dtype = int)]*self._particles_manager.n_particles

        self._sel_particle = None

        self._crystall_colors = np.zeros((self._ecal._n_patches,4))
        self._crystall_colors[:,0] = 1
        for n_particle in range(len(self._particles_manager._df)):
            self._crystall_colors[:,3] += self._particles_manager.get_crystall_content(n_particle)
        self.noise = np.clip(np.random.normal(0,self._noise_rate, self._ecal._n_patches),0,10000)

        self._crystall_colors[:,3] += self.noise
        
        self._crystall_content = deepcopy(self._crystall_colors[:,3])
        self._crystall_colors[:,3] = np.clip(2*np.sqrt(self._crystall_colors[:,3]),0,1)
        self._crystall_colors[np.where(self._crystall_colors[:,3]<=0.22),:]=[0,0,0,0.25]

        self._out = Output()
        with self._out:
            self._fig, self._ax = plt.subplots(figsize = (14, 14*60/89), constrained_layout = True)
        self._ax.set_ylim(-260, 340)
        self._ax.set_xlim(-445,445)
        self._ax.set_yticklabels([])
        self._ax.set_xticklabels([])

        self._crystall_collection = PatchCollection(self._ecal.patches, color = self._crystall_colors)
        self._crystall_artist = self._ax.add_collection(self._crystall_collection)

        self._selection_collection = LineCollection([])
        self._selection_artist = self._ax.add_collection(self._selection_collection)
        
        self._blit_manager = BlitManager(self._fig.canvas, self._crystall_artist, self._selection_artist)

        self._lasso = LassoSelector(self._ax, onselect = self.on_select)
        
        self._energy_labels = []
        box_list = []
        for i in range(self._particles_manager.n_particles):
            self._energy_labels.append(Text(description = "Gesamte Energie der ausgewählten Kristalle in GeV:", value = "0", disabled=True))
            box_list.append(HBox([self._energy_labels[i]]))
        
        self._particle_selector = Accordion(children=box_list,  titles = [f"Teilchen {str(i)}" for i in list(range(1,self._particles_manager.n_particles+1))], )
        self._particle_selector.observe(self.change_particle, names="selected_index")
        self._final_box = VBox(children=[self._particle_selector, self._out])
    
    def show(self) -> None:
        with self._out:
            plt.show()
        display(self._final_box)
        self.on_select()

    def change_particle(self, change) -> None:
        if self._particle_selector.selected_index is not None:
            self._sel_particle = self._particle_selector.selected_index
        else:
            self._sel_particle = None
        self.on_select()
    
    def on_select(self, verts = None) -> None:
        if (verts is not None) and (self._sel_particle is not None):
            print("verts")
            path = Path(verts)
            self._selected_crystalls[self._sel_particle] = np.nonzero(path.contains_points(self._ecal._patch_coordinates.T))[0]
            energy = np.sum(self._crystall_content[abs(self._crystall_content[self._selected_crystalls[self._sel_particle]])])
            self._particles_manager.energy_measurement(self._sel_particle, energy)
            self._energy_labels[self._sel_particle].value = str(round(self.energies[self._sel_particle], 5))
        edge_colors = np.clip(self._crystall_colors,0,1)
        for n in range(self._particles_manager.n_particles):
            edge_colors[self._selected_crystalls[n]] = [0,0,1,0.5] if n==self._sel_particle else [1,1,0,0.3]
        self._crystall_artist.set_edgecolors(edge_colors)
        self._blit_manager.update()