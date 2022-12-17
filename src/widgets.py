from faulthandler import disable
import ipywidgets as widgets
from matplotlib import pyplot as plt
import numpy as np
from src.particle import Particle
from src.ecal import ECal
from src.tracker import Tracker
import pandas as pd

from copy import deepcopy

from matplotlib.patches import FancyArrowPatch
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
from matplotlib.colors import to_rgba_array

def make_box_layout():
     return widgets.Layout(
        border='solid 1px black',
        margin='0px 10px 10px 0px',
        padding='5px 5px 5px 5px',
        height = "1000px",
        width = "1000px"
     )

def get_arrow(posx,posy,phi,size=1,granularity=100): #returns (granularity,2) x,y array of an arrow at position x,y pointing at angle phi 
    x=(np.array([0,0,-2,2,0,0])*np.cos(phi)-np.array([-3,3,0,0,3,-3])*np.sin(phi))*size
    y=(np.array([0,0,-2,2,0,0])*np.sin(phi)+np.array([-3,3,0,0,3,-3])*np.cos(phi))*size
    arrow = np.append(np.array([x,y]).T,np.zeros((granularity-6,2)),axis=0)
    arrow = arrow + [posx,posy]
    return arrow

def get_sphere(x, y, z):
    phi = np.arctan2(y, x)
    theta = np.arctan2(np.sqrt(x**2+y**2), z)
    return phi, theta

class BlitManager: #manages the blitting for tracker and ecal widget
    def __init__(self, canvas, artist, artist2=None):
        """copy from matplotlib website (blitting)"""
        self.canvas = canvas
        self._bg = None
        self.artist = artist
        if artist2 is not None:
            self.artist2=artist2
            self.twoartists=True
        else:
            self.twoartists=False
        
        # grab the background on every draw
        self.cid = canvas.mpl_connect("draw_event", self.on_draw)

    def on_draw(self, event):
        cv = self.canvas
        if event is not None:
            if event.canvas != cv:
                raise RuntimeError
        self._bg = cv.copy_from_bbox(cv.figure.bbox)
        self._draw_animated()


    def _draw_animated(self):
        fig = self.canvas.figure
        fig.draw_artist(self.artist)
        if self.twoartists:
            fig.draw_artist(self.artist2)

    def update(self):
        cv = self.canvas
        fig = cv.figure
        if self._bg is None:
            self.on_draw(None)
        else:
            cv.restore_region(self._bg)
            self._draw_animated()
            cv.blit(fig.bbox)
        cv.flush_events()

class TrackingWidget:
    def __init__(self, data_path, B = 0.2, layers = 15, n_segments = 8, ecl_segments = 30, k = 3, dist = 0.1, noise = 0.05, linewidth = 5, show_truthbutton = False, continuous_update=True, truthvalues=False, ignore_noise=False, trackercolor="gray"):
        self.continuous_update=continuous_update
        self.truthvalues=truthvalues
        self.granularity=100      #must be multiple of 4 otherwise everyone is gonna die
        self.show_truthbutton = show_truthbutton
        self.particles_df = pd.read_hdf(data_path)
        self.particles_df.loc[:,'charge'] = self.particles_df.loc[:,'pdg']/abs(self.particles_df.loc[:,'pdg'])
        self.particles_df.reset_index(inplace = True, drop = True)
        self.tracker = Tracker(layers = layers, n_segments = n_segments, ecl_segments=ecl_segments, k=k,dist=dist, noise = noise, linewidth = linewidth, ignore_noise = ignore_noise,granularity=self.granularity,trackercolor=trackercolor)
        self.n_particles = len(self.particles_df)
        self.B = B*15/layers
        self.particles_df.loc[:, "radius"] = self.particles_df.loc[:,"pt"]/self.B
        self.particles = []
        self.select_particles = []
        self.truth_particles = []
        self.index = 0
        self.arrows = []
        for i in range(self.n_particles):
            # build actual particles
            p_df = self.particles_df.iloc[i]
            p = Particle(p_df["radius"], p_df["phi"], self.B, p_df["charge"], granularity=self.granularity)
            self.truth_particles.append(p)

            # make an arrow to corresponding ecl crystal
            arrow_phi = self.tracker.get_arrowangle(p)
            arrow_x = np.cos(arrow_phi)*(self.tracker.layers+2.8)
            arrow_y = np.sin(arrow_phi)*(self.tracker.layers+2.8)
            self.arrows.append(get_arrow(posx=arrow_x,posy=arrow_y,phi = arrow_phi + np.pi/2,size=0.18,granularity=self.granularity))

            # build select particles used for simulation
            p = Particle(0.00001, 0, self.B, np.random.randint(0,1)*2-1, granularity=self.granularity)
            self.select_particles.append(p)

        self.tracker.make_tracker_mask(self.truth_particles) #color the hits of the true particles red 
            
    def change_particle(self,change):
        self.index = self.tabs.selected_index
        self.update(1)

    def update(self,change): #update function for drawing

        if self.index is None: #if no particle is selected, no particle is drawn
            drawtrace = False  
        else:  
            self.select_particles[self.index].phi = self.phi[self.index].value+self.phi_fine[self.index].value
            self.select_particles[self.index].charge = -1 if self.charge[self.index].value == "negative el. Ladung" else 1
            self.select_particles[self.index].radius = (self.r[self.index].value+self.r_fine[self.index].value)/self.B

            drawtrace = True
            if self.show_truthbutton:      
                if self.truthbutton.value:
                    trace=self.truth_particles[self.index].trace_array()
                    self.r_label[self.index].value = str(round(self.truth_particles[self.index].radius*self.B,6))+"(wahrer Wert)"
                    self.phi_label[self.index].value = str(round(self.truth_particles[self.index].phi,6))+"(wahrer Wert)"
                else:
                    trace=self.select_particles[self.index].trace_array()
                    self.r_label[self.index].value = str(round(self.select_particles[self.index].radius*self.B,6))
                    self.phi_label[self.index].value = str(round(self.select_particles[self.index].phi,6))
            else:
                trace=self.select_particles[self.index].trace_array()
                self.r_label[self.index].value = str(round(self.select_particles[self.index].radius*self.B,6))
                self.phi_label[self.index].value = str(round(self.select_particles[self.index].phi,6))

            hits,misses=self.tracker.get_hits_and_misses(self.select_particles[self.index],self.index)
            self.hit_n_misses[self.index].value = str(hits) + " hits & " + str(misses) + " misses"
        
        segments,colors=self.tracker.get_hit_lines(self.select_particles,self.index) #segments and colors of segments hit by the simulated particles 
        if drawtrace == True: #add the particle trace and an arrow to the things to be drawn
            segments=np.append(segments,[trace.T],axis=0)
            segments=np.append(segments,[self.arrows[self.index]],axis=0)
            colors=np.append(colors,["blue","green"])
        #to make drawing faster, the trackersegments, the trace and the arrow are all drawn by the same LineCollection artist

        self.artist.set_segments(segments)
        self.artist.set_color(colors)
        self.bm.update()

            
    def show(self):
        self.out = widgets.Output()
        with self.out:
            self.fig, self.ax = plt.subplots(figsize=(7,7),constrained_layout=True)

        limit = self.tracker.layers +3
        self.ax.set_xlim(-limit,limit)
        self.ax.set_ylim(-limit,limit)
        self.ax.set_axis_off()
        artist = self.ax.add_collection(LineCollection([]))
        self.artist=artist
        self.artist.set_animated(True)
        self.fig.canvas.toolbar_position = 'left'
        self.ax.add_collection(self.tracker.get_tracker_collection()) #draw tracker (gray and red segments)
        self.bm = BlitManager(self.fig.canvas , self.artist)

        self.hit_n_misses = []
        self.r_label = []
        self.r = []
        self.r_fine = []
        self.phi_label = []
        self.phi = []
        self.phi_fine = []
        self.charge = []
        self.box_list = []
        if self.show_truthbutton: #make all sliders and buttons
            self.truthbutton = widgets.ToggleButton(value = False, description = "Zeige wahres Teilchen")
            self.truthbutton.observe(self.update, names = "value")
        for i in range(self.n_particles):
            self.hit_n_misses.append(widgets.Text(description = "", value = "0 hits & 0 misses", disabled=True))
            self.r_label.append(widgets.Text(description = "$p_T$:", value = "0", disabled=True))
            self.r.append(widgets.FloatSlider(self.particles_df.loc[i,"pt"] if self.truthvalues == True else 0,min = 0, max = 5, step = 0.01, description = "$p_T$",continuous_update=self.continuous_update))
            self.r[i].observe(self.update, names = "value")
            self.r_fine.append(widgets.FloatSlider(0 ,min = 0, max = 0.2, step = 0.001, description = "$p_T fine$",continuous_update=self.continuous_update))
            self.r_fine[i].observe(self.update, names = "value")
            self.phi_label.append(widgets.Text(description = "$\phi$:", value = "0", disabled=True))
            self.phi.append(widgets.FloatSlider(self.particles_df.loc[i,"phi"] if self.truthvalues == True else 0,min = -np.pi, max = np.pi, step = 0.01, description = "$\phi$",continuous_update=self.continuous_update))
            self.phi[i].observe(self.update, names = "value")
            self.phi_fine.append(widgets.FloatSlider(0 ,min = -0.15, max = 0.15, step = 0.001, description = "$\phi fine$",continuous_update=self.continuous_update))
            self.phi_fine[i].observe(self.update, names = "value")
            self.charge.append(widgets.RadioButtons(options=['positive el. Ladung', 'negative el. Ladung'],  description=''))
            self.charge[i].observe(self.update, names = "value")
            if self.show_truthbutton:
                p_box = widgets.VBox([self.hit_n_misses[i],self.r_label[i],self.r[i], self.r_fine[i],self.phi_label[i], self.phi[i], self.phi_fine[i], self.charge[i], self.truthbutton])
            else:
                p_box = widgets.VBox([self.hit_n_misses[i],self.r_label[i],self.r[i], self.r_fine[i],self.phi_label[i], self.phi[i], self.phi_fine[i], self.charge[i]])
            self.box_list.append(p_box)
        
        self.tabs = widgets.Accordion()
        self.tabs.children = self.box_list
        for i in range(self.n_particles):
            self.tabs.set_title(i,f"particle {i}")
        self.tabs.observe(self.change_particle, names = "selected_index")
        self.tabs_box = widgets.HBox([self.tabs])
        self.tabs_box.layout =widgets.Layout(
                                border='solid 1px black',
                                margin='0px 10px 10px 0px',
                                padding='5px 5px 5px 5px',
                                height = "750px ",
                                width = "500px"
                            )
        self.plot_box = widgets.HBox([self.out])
        self.plot_box.layout = widgets.Layout(
                                border='solid 1px black',
                                margin='0px 10px 10px 0px',
                                padding='5px 5px 5px 5px',
                                height = "750px ",
                                width = "750px"
                            )
        self.final_box = widgets.HBox([  self.tabs_box,self.plot_box])
        with self.out:
            plt.show()
            #plt.pause(.1)
        self.fig.canvas.draw() 
        display(self.final_box)  
        self.update(1)   

    @property
    def get_fitted_particles(self):
        df = pd.DataFrame(columns = ["pt", "phi", "el. Ladung", "radius"])
        for i in range(self.n_particles):
            df.loc[i,:] = [self.select_particles[i].momentum(), self.select_particles[i].phi, self.select_particles[i].charge, self.select_particles[i].radius]
        df.loc[:,"px"] = np.cos(df.loc[:,"phi"].astype("float"))*df.loc[:,"pt"]
        df.loc[:,"py"] = np.sin(df.loc[:,"phi"].astype("float"))*df.loc[:, "pt"]
        df.loc[:, "pz"] = self.particles_df.loc[:,"pz"]
        return df

    @property
    def get_ecl_position(self):
        phi=[]
        theta=[]
        r=0
        for i in range(self.n_particles):
            phi_0=self.select_particles[i].charge*np.pi/2-self.select_particles[i].phi+np.arccos(r/(2*self.select_particles[i].radius))*self.select_particles[i].charge-np.pi/2
            if phi_0 < 0:
                phi_0 += 2*np.pi  
            if phi_0 < 0:
                phi_0 += 2*np.pi  
            if phi_0 > np.pi:
                phi_0 -= 2*np.pi  
            if phi_0 > np.pi:
                phi_0 -= 2*np.pi                
            phi.append(phi_0)
            theta.append(np.arctan2(self.select_particles[i].momentum(),self.particles_df.loc[i,"pz"])) 
        return np.array(phi)*180/np.pi,np.array(theta)*180/np.pi

class TestDetektor:
    def __init__(self, B=0.1, layers=8, n_segments=3,ecl_segments=10, k=2):
        if layers > 20:
            print("Es sind maximal 20 Schichten möglich!")
            layers = 20
        self.tracker = Tracker(layers = layers, n_segments = n_segments,k=k ,ecl_segments=ecl_segments,noise = False, linewidth = 2)
        self.B = B
        self.particle = Particle(1, 0, B,-1)
        self.pt = 10

    def update(self,change):
        self.particle.charge = -1 if self.charge_widget.value == "negative el. Ladung" else 1
        self.B = self.b_widget.value/5
        self.particle.B = self.B    
        self.particle.radius = self.pt_widget.value/self.B if self.B != 0 else 100000
        trace=self.particle.trace_array()
        segments,colors=self.tracker.get_hit_lines([self.particle],0)
        segments=np.append(segments,[trace.T],axis=0)
        colors=np.append(colors,["blue"])
        self.artist.set_segments(segments)
        self.artist.set_color(colors)
        self.bm.update()

    def show(self):
        self.fig, self.ax = plt.subplots(figsize=(7,7))

        lim = self.tracker.layers+3
        self.ax.set_ylim([-lim,lim])
        self.ax.set_xlim([-lim,lim])
        artist = self.ax.add_collection(LineCollection([]))
        self.artist=artist
        self.artist.set_animated(True)
        self.fig.canvas.toolbar_position = 'left'
        self.ax.add_collection(self.tracker.get_tracker_collection())
        self.bm = BlitManager(self.fig.canvas , self.artist)        
        
        self.pt_widget= widgets.FloatSlider(1 ,min = 0.1, max = 4, step = 0.1, description = f"$p_t$")
        self.pt_widget.observe(self.update, names = "value")
        self.b_widget= widgets.Checkbox((False), description = "B-Feld")
        self.b_widget.observe(self.update, names = "value")
        self.charge_widget= widgets.RadioButtons(options=['positive el. Ladung', 'negative el. Ladung'],  description='')
        self.charge_widget.observe(self.update, names = "value")
        self.out = widgets.Output()
        p_box = widgets.VBox([self.pt_widget, self.b_widget,self.charge_widget])
        plt.show()
        display(self.out, p_box)
        self.update(1)    
        self.update(1)   

  
class ECLWidget:

    def __init__(self, data_path, noise_rate = 0.05, tw=None, idx=None):
        data = pd.read_hdf(data_path)
        if tw is not None:
            self.line_position=-tw.get_ecl_position[0]*2+890
            mask=self.line_position>720
            self.line_position[mask]-=720
            mask=self.line_position<0
            self.line_position[mask]+=720
            self.drawline=True
        else:
            self.drawline=False
        coords = [f'{i}' for i in np.arange(0, 6624)]
        if idx is None:
            hits = data[coords]
        else:
            hits = data[coords].iloc[idx:(idx+1)]
        hits = hits.reset_index(drop=True)
        self.edge_size = 5
        self.ecal = ECal(144,46,hits, crystal_edge = self.edge_size, noise_rate = noise_rate)   
        content = np.array(self.ecal.crystals_df["content"].to_numpy(),dtype=float)
        content = 2*np.sqrt(abs(content))
        self.alphas = np.clip(content,0.2,1)
        
        self.out = widgets.Output()
        with self.out:
            fig, ax = plt.subplots(figsize=(15,6),constrained_layout=True)
        ax.set_ylim(-10,46*5+10)
        ax.set_xlim(-10,144*5+10)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        self.artist = ax.add_collection(self.ecal.collection)
        self.lineartist=ax.add_collection(LineCollection([]))
        self.xys = np.array([self.ecal.crystals_df["x"] + self.edge_size/2,self.ecal.crystals_df["y"] + self.edge_size/2],dtype = "float64").T
        self.Npts = len(self.xys)
        self.lasso = LassoSelector(ax, onselect=self.onselect)
        self.artist.set_animated(True)
        self.lineartist.set_animated(True)
        self.bm_ecal = BlitManager(fig.canvas , self.artist, self.lineartist)
        self.allhidden=False
        self.ind = []
        self.vertices = [[(0,0)]]*len(hits)
        self.particle_index = 0
        
    def onselect(self, verts):
        self.vertices[self.particle_index]=verts
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.ecal.select_particles.loc[self.particle_index, :] = 0
        self.ecal.select_particles.loc[self.particle_index, self.ind.astype(str)] = 1
        self.ecal.set_colors(self.particle_index,self.allhidden,self.drawline)
        facecolors = to_rgba_array(self.ecal.crystals_df.loc[:,"facecolor"].to_numpy())
        content_mask = (self.ecal.crystals_df["content"]>0).to_numpy()
        facecolors.T[-1] = 0.5
        facecolors[content_mask,-1] = self.alphas[content_mask]
        edgecolors = to_rgba_array(self.ecal.crystals_df.loc[:,"edgecolor"].to_numpy())
        self.ecal.collection.set_edgecolors(edgecolors)
        self.ecal.collection.set_facecolors(facecolors)
        if self.drawline==True and self.allhidden==False:
            self.lineartist.set_segments([np.array([np.ones(100)*self.line_position[self.particle_index],np.linspace(-5,46*5+5,100)]).T])
            self.lineartist.set_color(["blue"])
        self.bm_ecal.update()
        particle_mask = self.ecal.select_particles.loc[self.particle_index, :].to_numpy()>0
        energy = self.ecal.crystals_df.loc[particle_mask, "content"].sum()
        self.energy_labels[self.particle_index].value = str(round(energy,4))
        
    def change_particle(self,change):
        if self.particle.selected_index is not None:
            self.particle_index = self.particle.selected_index
            self.allhidden=False
        else:
            self.allhidden=True
        self.onselect(self.vertices[self.particle_index])
        
    def show(self):
        self.particle = widgets.Accordion()
        self.particle.observe(self.change_particle, names = "selected_index")
        self.energy_labels = []
        self.box_list = []
        for i in range(self.ecal.n_particles):
            self.particle.set_title(i, f"Teilchen {i}")
            self.energy_labels.append(widgets.Text(description = "Gesamte Energie der ausgewählten Kristalle in GeV:", value = "0", disabled=True))
            self.box_list.append(widgets.HBox([self.energy_labels[i]]))
        self.particle.children = self.box_list
        self.final_box = widgets.VBox(children=[self.particle, self.out])
        with self.out:
            plt.show()
        display(self.final_box)
        self.onselect([(0,0)])
    
    @property
    def get_particles_energy(self):
        energys = []
        for i in range(self.ecal.n_particles):
            particle_mask = self.ecal.select_particles.loc[i, :].to_numpy()>0
            energys.append(self.ecal.crystals_df.loc[particle_mask, "content"].sum())
        return pd.DataFrame(energys, columns = ["Energie"])
    
    @property
    def get_particles_radius(self):
        radius = []
        for i in range(self.ecal.n_particles):
            particle_mask = self.ecal.select_particles.loc[i, :].to_numpy()>0
            selected_crystals = self.ecal.crystals_df.loc[particle_mask]
            selected_hits = selected_crystals.query("content>0").index
            xdiff = abs(int(self.ecal.center[i]/144.)%46 - (selected_hits/144).astype('int')%46)
            ydiff = abs(self.ecal.center[i]%144 - selected_hits%144)
            radius.append(max(xdiff.max(), ydiff.max()))
        return pd.DataFrame(radius, columns= ["Radius"])

    @property
    def get_selections(self):
        selectionvertices=[]
        points=[]
        xypoints=[]
        patches=[]
        minxs=[1000]*len(self.vertices)
        maxxs=[0]*len(self.vertices)
        minys=[1000]*len(self.vertices)
        maxys=[0]*len(self.vertices)
        colors=[]
        for i in range(len(self.vertices)):
            for l in range(len(self.vertices[i])):
                if self.vertices[i][l][0]<minxs[i]:
                    minxs[i]=self.vertices[i][l][0]
                if self.vertices[i][l][0]>maxxs[i]:
                    maxxs[i]=self.vertices[i][l][0]
                if self.vertices[i][l][1]<minys[i]:
                    minys[i]=self.vertices[i][l][1]
                if self.vertices[i][l][1]>maxys[i]:
                    maxys[i]=self.vertices[i][l][1]         
            midx=(maxxs[i]+minxs[i])/2
            midy=(maxys[i]+minys[i])/2
            dummypatches=[]
            selectionvertices.append([(midx-45,midy-45),(midx+45,midy-45),(midx+45,midy+45),(midx-45,midy+45)])
            selpath = Path(selectionvertices[i])
            points.append(np.nonzero(selpath.contains_points(self.xys))[0])
            xypoints.append(self.xys[points[i]])
            midx = (np.amax(xypoints[i][:,0])-np.amax(-xypoints[i][:,0]))/2
            midy = (np.amax(xypoints[i][:,1])-np.amax(-xypoints[i][:,1]))/2
            for l in range(len(xypoints[i])):
                dummypatches.append(Rectangle((xypoints[i][l][0]-midx-2, xypoints[i][l][1]-midy-2), 4, 4, edgecolor = "black", facecolor = "gray", linewidth = 3))
            patches.append(dummypatches)
            hitmask=(self.ecal.crystals_df["content"] > 0)[points[i]]
            dummycolors=np.array(["gray"]*len(xypoints[i]))
            dummycolors[hitmask]= "red"
            dummycolors = to_rgba_array(dummycolors)
            dummycolors.T[-1] = 0.5
            dummycolors[hitmask,-1] = self.alphas[points[i]][hitmask]
            colors.append(dummycolors)
        return patches,colors


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
truth_particles = pd.DataFrame(columns = ["Masse", "el. Ladung","Image","K_L0","E_p"], data=true_particle_data, index=true_particle_names)
truth_particles.loc[:, "Masse"] = truth_particles["Masse"]*10**(-3)

class MatchingWidget:
    def __init__(self, ew, tw, kl=None, cheat_mode=True, cheating_threshhold = 1e-2) -> None:
        self.energies = ew.get_particles_energy
        self.radius = ew.get_particles_radius
        self.momenta = tw.get_fitted_particles
        if kl is not None:
            self.KLM_hits = np.array(kl.KLM_hit)
        else:
            self.KLM_hits = np.array([0]*len(self.momenta))
        self.patches=[]
        for i in range(len(ew.get_selections[0])):
            self.patches.append(PatchCollection(np.array(ew.get_selections[0][i]),color=ew.get_selections[1][i]))
        columns = ["Ladung", "Energie", "impuls", "Masse", "E_p","KLM"]
        self.true_df = tw.particles_df
        self.res_df = pd.DataFrame(data = np.zeros((len(self.energies), len(columns))), columns = columns)
        #self.res_df.loc[:,"pt"] = np.sqrt(( self.momenta.loc[:,["px", "py"]]**2).sum())
        if(cheat_mode):
            self.momenta_cheat_mask = ((self.true_df["pt"]-self.momenta["pt"])<cheating_threshhold).to_numpy()
            self.energies_cheat_mask = ((self.true_df["energy"]-self.energies["Energie"])<cheating_threshhold).to_numpy()
            self.momenta.loc[self.momenta_cheat_mask, ["px", "py", "pz"]] =  self.true_df.loc[self.momenta_cheat_mask, ["px", "py", "pz"]]
            self.energies.loc[self.energies_cheat_mask, "Energie"] = self.true_df.loc[self.energies_cheat_mask, "energy"] 

            
    def update(self, change = 0):
        sele_index = self.tabs.selected_index
        self.res_df.loc[sele_index, "Energie"] = self.energies.loc[sele_index, "Energie"]
        self.res_df.loc[sele_index, "el. Ladung"] = self.momenta.loc[sele_index, "el. Ladung"]
        self.res_df.loc[sele_index, "impuls"] = np.sqrt((self.momenta.loc[sele_index, ["px", "py", "pz"]]**2).sum().astype("float"))
        self.res_df.loc[sele_index, "E_p"] = self.res_df.loc[sele_index, "Energie"]/self.res_df.loc[sele_index, "impuls"]
        self.res_df.loc[sele_index, "KLM"] = self.KLM_hits[sele_index]
        # if self.res_df.loc[:, "Energie"] > self.res_df.loc[:, "impuls"]:
        self.res_df.loc[:, "Masse"] = np.sqrt(abs(self.res_df.loc[:, "Energie"]**2 - self.res_df.loc[:, "impuls"]**2))
        self.res_df.loc[:, "Masse"] = self.res_df.loc[:, "Masse"].fillna(0)
        self.sel_charge[sele_index].value = str(truth_particles.loc[self.part_ids[sele_index].value, "el. Ladung"])
        self.sel_mass[sele_index].value = str(truth_particles.loc[self.part_ids[sele_index].value, "Masse"])
        self.sel_image[sele_index].value = truth_particles.loc[self.part_ids[sele_index].value, "Image"]
        self.sel_label[sele_index].value = "So sieht ein typisches "+self.part_ids[sele_index].value + " Teilchen im Ecal aus:"
        self.sel_KL0[sele_index].value = truth_particles.loc[self.part_ids[sele_index].value, "K_L0"] + " im KLM Detektor"
        self.sel_E_p[sele_index].value = truth_particles.loc[self.part_ids[sele_index].value, "E_p"]

        #for i in range(len(self.res_df)):
        self.KL0_txt[sele_index].value = "kein Hit im KLM Detektor" if self.res_df.loc[sele_index, "KLM"] == 0 else "Hit im KLM Detektor"
        self.energy_txt[sele_index].value = str(self.res_df.loc[sele_index, "Energie"])
        self.charge_txt[sele_index].value = str(self.res_df.loc[sele_index, "el. Ladung"])
        self.moment_txt[sele_index].value = str(self.res_df.loc[sele_index, "impuls"])
        self.invmas_txt[sele_index].value = str(self.res_df.loc[sele_index, "Masse"])
        self.E_p_txt[sele_index].value = str(self.res_df.loc[sele_index, "E_p"])
        self.px_txt[sele_index].value = str(self.momenta.loc[sele_index, "px"])
        self.py_txt[sele_index].value = str(self.momenta.loc[sele_index, "py"])
        self.pz_txt[sele_index].value = str(self.momenta.loc[sele_index, "pz"])
        self.ax1.cla()
        self.ax1.set_ylim(-50,50)
        self.ax1.set_xlim(-50,50)
        self.ax1.set_yticklabels([])
        self.ax1.set_xticklabels([])
        self.ax1.add_collection(self.patches[sele_index])

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
        self.out1 = widgets.Output()
        with self.out1:
            self.fig1, self.ax1 = plt.subplots(figsize=(3.4,3.4),constrained_layout=True)
        self.ax1.set_ylim(-3,103)
        self.ax1.set_xlim(-3,103)
        self.label1=widgets.Text(value = "Resultate", disabled = True)
        self.label2=widgets.Text(value = "Bekannte Teilchen zum Vergleichen", disabled = True)

        for i in range(len(self.res_df)):
            self.px_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "$p_x$", disabled = True))
            self.py_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "$p_y$", disabled = True))
            self.pz_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "$p_z$", disabled = True))
            self.energy_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "Energie", disabled = True))
            self.charge_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "el. Ladung", disabled = True))
            self.moment_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "Impuls", disabled = True))
            self.invmas_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "Masse", disabled = True))
            self.E_p_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "$E/p$", disabled = True))
            self.KL0_txt.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "$K_L^0$", disabled = True))   
            self.res_box = widgets.VBox(children=[self.label1, self.energy_txt[i], self.charge_txt[i], self.moment_txt[i], self.invmas_txt[i], self.E_p_txt[i],self.px_txt[i],self.py_txt[i],self.pz_txt[i],self.KL0_txt[i],self.out1])
            self.res_box.layout = widgets.Layout(border='solid 1px black',margin='0px 10px 10px 0px',padding='5px 5px 5px 5px',height = "700px ",width = "370px")            

            self.part_ids.append(widgets.Select(options = truth_particles.index, value = "e+", description = "Teilchen"))
            self.part_ids[i].observe(self.update, "value")
            self.sel_mass.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "Masse", disabled = True))
            self.sel_charge.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "el. Ladung", disabled = True))
            self.sel_E_p.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "$E/p$", disabled = True))
            self.sel_KL0.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", description = "$K_L^0$", disabled = True))
            self.sel_label.append(widgets.Text(placeholder = "kein Teilchen ausgewählt", disabled = True))
            self.sel_image.append(widgets.Image(value=truth_particles.loc["e+", "Image"],format='png',width=320,height=320))
            self.sel_box = widgets.VBox(children=[self.label2, self.part_ids[i], self.sel_mass[i], self.sel_charge[i], self.sel_KL0[i],self.sel_E_p[i],self.sel_label[i],self.sel_image[i]])
            self.sel_box.layout = widgets.Layout(border='solid 1px black',margin='0px 10px 10px 0px',padding='5px 5px 5px 5px',height = "700px ",width = "370px")  

            box = widgets.HBox(children=[self.res_box, self.sel_box])
            boxes.append(box)
        self.tabs = widgets.Tab(children=boxes)
        self.tabs.observe(self.update, "selected_index")
        for i in range(len(self.res_df)):
            self.tabs.set_title(i,f"Teilchen {i}")
        self.update()
        with self.out1:
            plt.show()
        display(self.tabs)

class MissingWidget():

    def calc_missing_part(self, dummy):
        fourvecs = np.zeros((4,4))
        for i in range(4):
            for j in range(4):
                fourvecs[i,j] = self.boxes[i].children[j+1].value
        missing_part = np.array([15.580,0,0,0])-fourvecs.sum(0)
        print(missing_part)
        for j in range(4):
            self.boxes[4].children[j+1].value = missing_part[j]
        self.boxes[4].children[5].value = np.sqrt(((missing_part*np.array([1,-1,-1,-1]))**2).sum())


    def show(self):
        self.boxes = []
        for i in range(4):
            w_label1 = widgets.Label(value=f"Teilchen {i}")
            w_energy1 = widgets.FloatText(description="Energie", placeholder = "kein Eintrag", value = 0., layout = widgets.Layout(width="200px"))
            w_px1 = widgets.FloatText(description="$p_x$", placeholder = "kein Eintrag", value = 0., layout = widgets.Layout(width="200px"))
            w_py1 = widgets.FloatText(description="$p_y$", placeholder = "kein Eintrag", value = 0., layout = widgets.Layout(width="200px"))
            w_pz1 = widgets.FloatText(description="$p_z$", placeholder = "kein Eintrag", value = 0., layout = widgets.Layout(width="200px"))
            t1_box = widgets.VBox(children=[w_label1, w_energy1, w_px1, w_py1, w_pz1])
            self.boxes.append(t1_box)
        w_label1 = widgets.Label(value=f"fehlendes Teilchen")
        w_energy1 = widgets.FloatText(description="Energie", placeholder = "kein Eintrag", disabled = True, layout = widgets.Layout(width="200px"))
        w_px1 = widgets.FloatText(description="$p_x$", placeholder = "kein Eintrag", disabled = True, layout = widgets.Layout(width="200px"))
        w_py1 = widgets.FloatText(description="$p_y$", placeholder = "kein Eintrag", disabled = True, layout = widgets.Layout(width="200px"))
        w_pz1 = widgets.FloatText(description="$p_z$", placeholder = "kein Eintrag", disabled = True, layout = widgets.Layout(width="200px"))
        w_mass = widgets.FloatText(description="Masse",  placeholder = "kein Eintrag", disabled = True)
        t1_box = widgets.VBox(children=[w_label1, w_energy1, w_px1, w_py1, w_pz1, w_mass])
        self.boxes.append(t1_box)
        button = widgets.Button(description="Berechne fehlendes Teilchen")
        button.on_click(self.calc_missing_part)
        self.boxes.append(button)
        self.final_box = widgets.HBox(children=self.boxes)
        display(self.final_box)


class KLMWidget:
    def __init__(self,data_path,tw=None,always_hit=False,a=35,b=1.2):
        self.rows=64
        self.columns=200
        edge_gap_ratio=2/4
        self.blocksize=4
        self.edgesize=self.blocksize*edge_gap_ratio
        self.gapsize=self.blocksize*(1-edge_gap_ratio)        
        data = pd.read_hdf(data_path)
        self.data=[]
        self.detec=[]
        for i in range(len(data)):
            self.data.append([int((np.argmax(data.iloc[i,np.arange(3,6627)].to_numpy())%144)*(self.columns/144)),int(int(np.argmax(data.iloc[i,np.arange(3,6627)].to_numpy())/144)*(self.rows/46))])
            if (data.iloc[i,6628]==13 or data.iloc[i,6628]==-13 or always_hit):
                self.detec.append(True)
            else:
                self.detec.append(False)

        if tw is None:
            phi_position=[]
            theta_position=[]
            for i in range(len(self.data)):
                phi_position.append(self.data[i][0])
                theta_position.append(self.data[i][1])
        else:
            phi_position=-tw.get_ecl_position[0]*2+890
            mask=phi_position>720
            phi_position[mask]-=720
            mask=phi_position>720
            phi_position[mask]-=720
            mask=phi_position<0
            phi_position[mask]+=720
            phi_position = phi_position/5/144*self.columns
            theta_position=(tw.get_ecl_position[1]-33)/1.45*(self.rows/64)

        self.x_position = np.array(phi_position)*self.blocksize
        self.y_position = np.array(theta_position)*self.blocksize
        self.patches=[]
        self.patchcolors=np.zeros((self.columns*self.rows,4))
        self.patchcolors[:,2]=1

        for c in range(self.columns):
            for r in range(self.rows):
                self.patches.append(Rectangle((c*self.blocksize,r*self.blocksize), self.edgesize, self.edgesize, edgecolor = "blue", facecolor = "gray", linewidth = 1))     
                intensity=0
                for i in range(len(self.data)):
                    if self.detec[i]:
                        intensity=intensity+np.exp((-(c-self.data[i][0])**2-(r-self.data[i][1])**2)*0.04)
                self.patchcolors[c*self.rows+r,3]=np.random.normal(loc=intensity, scale=0.13)

        self.patchcolors[:,3] = np.clip(self.patchcolors[:,3],0.1,1)
        
        self.out = widgets.Output()
        with self.out:
            fig, ax = plt.subplots(figsize=(15,6),constrained_layout=True)

        ax.set_yticklabels([])
        ax.set_xticklabels([])

        ax.set_ylim(-10,self.rows*self.blocksize+10-self.edgesize)
        ax.set_xlim(-10,self.columns*self.blocksize+10-self.edgesize)
        self.patchartist = ax.add_collection(PatchCollection(self.patches,color=self.patchcolors))
        self.lineartist = ax.add_collection(LineCollection([]))
        self.patchartist.set_animated(True)
        self.lineartist.set_animated(True)
        self.bm = BlitManager(fig.canvas , self.patchartist, self.lineartist)
        
    def update(self, change):
        if self.tabs.selected_index is None:
            segments=[]
        else:
            segments=[[[self.x_position[self.tabs.selected_index]-10*self.blocksize,self.y_position[self.tabs.selected_index]-10*self.blocksize],
                       [self.x_position[self.tabs.selected_index]-10*self.blocksize,self.y_position[self.tabs.selected_index]+10*self.blocksize],
                       [self.x_position[self.tabs.selected_index]+10*self.blocksize,self.y_position[self.tabs.selected_index]+10*self.blocksize],
                       [self.x_position[self.tabs.selected_index]+10*self.blocksize,self.y_position[self.tabs.selected_index]-10*self.blocksize],
                       [self.x_position[self.tabs.selected_index]-10*self.blocksize,self.y_position[self.tabs.selected_index]-10*self.blocksize]]]
        
        self.lineartist.set_segments(segments)
        self.lineartist.set_color("red")
        self.bm.update()
        
    def show(self):

        self.tabs = widgets.Accordion()
        self.tabs.observe(self.update, names = "selected_index")
        self.tickbox = []
        self.box_list = []
        self.boxtext=widgets.Text(value = "Wurde hier ein Teilchen erkannt?", disabled = True)
        for i in range(len(self.data)):
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