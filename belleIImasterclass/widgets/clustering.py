from matplotlib.path import Path
from matplotlib.widgets import LassoSelector
from belleIImasterclass.particlesmanager import ParticlesManager
from belleIImasterclass.elements.ecal import ECal
from belleIImasterclass.widgets.blitmanager import BlitManager
from ipywidgets import Output, Accordion, Text, HBox, VBox, Button
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
from matplotlib.transforms import Affine2D
from matplotlib.patches import Rectangle

class ECLWidget:
    def __init__(self,particles_manager: ParticlesManager, noise_ratio = 0.2, true_particles=False) -> None:
        self._particles_manager = particles_manager
        self._noise_ratio = noise_ratio
        self._true_particles=true_particles
        self._ecal = ECal()        
        self._selected_crystalls = [np.array([], dtype = int)]*self._particles_manager.n_particles

        self._sel_particle = None
        self._crystall_colors = np.zeros((self._ecal._n_patches,4))
        self._crystall_colors[:,0] = 1
        for n_particle in range(len(self._particles_manager._df)):
            self._crystall_colors[:,3] = self._crystall_colors[:,3] + self._particles_manager.get_crystall_content(n_particle)

        noise_std=0.008
        noise_mean=0.004
        self.noise = np.clip(np.random.normal(noise_mean,noise_std, self._ecal._n_patches),0,0.2)*np.random.choice(a=[1, 0], size=(self._ecal._n_patches), p=[self._noise_ratio, 1-self._noise_ratio])
        self._crystall_colors[:,3] += self.noise
        
        self._crystall_content = deepcopy(self._crystall_colors[:,3])
        self._crystall_colors[:,3] = np.clip(1.5*np.sqrt(self._crystall_colors[:,3]),0,1)
        self._crystall_colors[np.where(self._crystall_colors[:,3]<=0.1),:]=[0,0,0,0.1]

        self._out = Output()
        
        self._energy_labels = []
        box_list = []
        self.update_button = Button(description='Update!',disabled=False,tooltip='Update',icon='rotate-right')
        for i in range(self._particles_manager.n_particles):
            self._energy_labels.append(Text(description = "Gesamte Energie der ausgewählten Kristalle in GeV:", value = "0", disabled=True))
            box_list.append(HBox([self._energy_labels[i],self.update_button]))
        
        self._particle_selector = Accordion(children=box_list,  titles = [f"Teilchen {str(i)}" for i in list(range(self._particles_manager.n_particles))], )
        self._particle_selector.observe(self.update, names="selected_index")
        self.update_button.on_click(self.update)
        self._final_box = VBox(children=[self._particle_selector, self._out])
    
    def show(self) -> None:
        with self._out:
            self._fig, self._ax = plt.subplots(figsize = (14, 14*60/89), constrained_layout = True)
        self._ax.set_xlim(-450, 450)
        self._ax.set_ylim(-280, 340)
        self._ax.set_yticklabels([])
        self._ax.set_xticklabels([])

        self._crystall_collection = PatchCollection(self._ecal.patches, color = self._crystall_colors)
        self._crystall_artist = self._ax.add_collection(self._crystall_collection)
        self._crystall_artist.set_aa(True)
        self._crystall_artist.set_linewidth(1.4)
        self._crystall_artist.set_animated(True)

        self.circle_collection = PatchCollection([])
        self._circle_artist = self._ax.add_collection(self.circle_collection)
        self._circle_artist.set_aa(True)
        self._circle_artist.set_linewidth(2)
        self._circle_artist.set_animated(True)

        self._blit_manager = BlitManager(self._fig.canvas, self._crystall_artist, self._circle_artist)

        self._lasso = LassoSelector(self._ax, onselect = self.on_select)
        
        with self._out:
            plt.show()
        display(self._final_box)
        self.update()


    def update(self, change=0) -> None:
        if self._particle_selector.selected_index is not None:
            self._sel_particle = self._particle_selector.selected_index
        else:
            self._sel_particle = None

        if self._true_particles:
            phi = self._particles_manager._df["phi"].to_numpy()[self._particles_manager.index]
            theta = self._particles_manager._df["theta"].to_numpy()[self._particles_manager.index] 
        else:
            phi = self._particles_manager._df["tracker_phi"].to_numpy()[self._particles_manager.index]
            theta = np.arctan2(self._particles_manager._df["tracker_pt"].to_numpy(),self._particles_manager._df["pz"].to_numpy())[self._particles_manager.index]

        self.circle_coordinates=self._ecal.get_cell_coordinates(theta,phi)

        self.on_select()
    
    def on_select(self, verts = None) -> None:
        if (verts is not None) and (self._sel_particle is not None):
            path = Path(verts)
            self._selected_crystalls[self._sel_particle] = np.nonzero(path.contains_points(self._ecal._patch_coordinates.T))[0]
            energy = np.sum(abs(self._crystall_content[self._selected_crystalls[self._sel_particle]]))
            self._particles_manager.energy_measurement(self._sel_particle, energy)
            self._energy_labels[self._sel_particle].value = str(round(energy, 5))

        if (self._sel_particle is not None) and ((self._true_particles) or (self._particles_manager._df.loc[self._sel_particle,"tracker_pt"] != 0)):
            circle=Circle((self.circle_coordinates[0,self._sel_particle],self.circle_coordinates[1,self._sel_particle]),radius=50)
            self._circle_artist.set_paths([circle])
            self._circle_artist.set_facecolor([[0,0,0,0]])
            self._circle_artist.set_edgecolor([[0,0,1,1]])

        edge_colors = np.clip(self._crystall_colors,0,1)
        for n in range(self._particles_manager.n_particles):
            edge_colors[self._selected_crystalls[n]] = [0,0,1,0.5] if n==self._sel_particle else [1,1,0,0.3]
        self._crystall_artist.set_edgecolors(edge_colors)
        self._blit_manager.update()
        if self._sel_particle is not None:
            self.make_patches(self._sel_particle)

    def make_patches(self,index,rand=20):
        selection_corners=np.zeros((2,5))
        midpoint=np.zeros((2))
        selection_size=1
        ptchcrds=self._ecal._patch_coordinates[:,self._selected_crystalls[index]]
        if len(ptchcrds[0])>0:
                midpoint[0]=(np.amax(ptchcrds[0])-np.amax(-ptchcrds[0]))/2
                midpoint[1]=(np.amax(ptchcrds[1])-np.amax(-ptchcrds[1]))/2           
                selection_size=np.amax(ptchcrds[0])+np.amax(-ptchcrds[0]) 
                if selection_size < (np.amax(ptchcrds[1])+np.amax(-ptchcrds[1])):  
                    selection_size = np.amax(ptchcrds[1])+np.amax(-ptchcrds[1])
                selection_corners[0,0]=midpoint[0]-selection_size/2-rand
                selection_corners[1,0]=midpoint[1]-selection_size/2-rand
                selection_corners[0,1]=midpoint[0]+selection_size/2+rand
                selection_corners[1,1]=midpoint[1]-selection_size/2-rand
                selection_corners[0,2]=midpoint[0]+selection_size/2+rand
                selection_corners[1,2]=midpoint[1]+selection_size/2+rand
                selection_corners[0,3]=midpoint[0]-selection_size/2-rand
                selection_corners[1,3]=midpoint[1]+selection_size/2+rand
                selection_corners[0,4]=selection_corners[0,0]
                selection_corners[1,4]=selection_corners[1,0]    

        path = Path(selection_corners.T)
        patchindices=np.nonzero(path.contains_points(self._ecal._patch_coordinates.T))[0]
        patches=[]
        facecolors=np.clip(self._crystall_colors[patchindices],0,1)
        edgecolors=np.clip(self._crystall_colors,0,1)
        edgecolors[self._selected_crystalls[index]]=[0,0,1,0.5]
        edgecolors=edgecolors[patchindices]
        for l in range(len(patchindices)):
            patches.append(Rectangle((self._ecal._patch_coordinates[:,patchindices[l]]+self._ecal._patch_offsets[:,patchindices[l]]-midpoint)*(50/(selection_size+2*rand)),
                                     width=self._ecal._crystal_size*(50/(selection_size+2*rand)),height=self._ecal._crystal_size*(50/(selection_size+2*rand)),
                                     angle=self._ecal._patch_angles[patchindices[l]]*180/np.pi,linewidth=20,))
        self._particles_manager.ecal_patches(index,patches,facecolors,edgecolors)


