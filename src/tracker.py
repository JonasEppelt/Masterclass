from copy import deepcopy
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd

def get_track_line(radius, begin, end, granularity = 100):
    t = np.linspace(begin, end, granularity)
    return radius*np.cos(t), radius*np.sin(t)

def get_ecl_lines(radius, begin, end, width, granularity = 100):
        t = np.linspace(begin, end, int(granularity/4))   
        ti = np.linspace(end, begin, int(granularity/4))
        inner_x = radius*np.cos(t)
        inner_y = radius*np.sin(t)
        outer_x = (radius+width)*np.cos(ti)
        outer_y = (radius+width)*np.sin(ti)
        
        right_x = np.linspace(inner_x[-1], outer_x[0], int(granularity/4))
        right_y = np.linspace(inner_y[-1], outer_y[0], int(granularity/4))
        left_x = np.linspace(outer_x[-1], inner_x[0], int(granularity/4))
        left_y = np.linspace(outer_y[-1], inner_y[0], int(granularity/4))

        retx=np.append(inner_x,[right_x,outer_x,left_x])
        rety=np.append(inner_y,[right_y,outer_y,left_y])

        return retx,rety

class Tracker:
    def __init__(self, layers, n_segments, ecl_segments, k = 2, dist = 0.2, noise = False, linewidth = 8, ignore_noise = False, granularity=100, trackercolor="gray"):
        self.layers = layers                                                    #number of layers in the tracker (including ecl segments)
        self.noise = noise                                                      #noisefloor constant
        self.trackercolor=trackercolor
        self.granularity = granularity                                          #number of points in the segment lines                   
        self.ignore_noise = ignore_noise                                        #only important for hits&misses , if true noise hits will not be counted
        self.total_segments=0                                                   #will contain number segments aufter initialisation
        self.n_segments = n_segments - k                                        #number of segments in the first layer
        self.ecl_segments = ecl_segments                                        #number of ecl segments
        self.segments = pd.DataFrame(columns = ["begin","end","radius"])        #dataframe of each segments radius, starting angle and ending angle
        self.all_lines = []                                                     #this will contain the x,y arrays for all segments in one single array in the shape (total_segments,2,granularity)
        self.linewidths = []                                                    #contains the linewidths for all segments
        self.noisemask = []                                                     #mask for noise hits
        for l in range(1,layers+1):                                             #making the normal segments
            len_segment = 2*np.pi/(n_segments+k*l)
            for i in range(n_segments+k*l):
                self.total_segments+=1
                radius = l
                begin = len_segment*i+dist/(l+1)#+0.1
                end = len_segment*(i+1)-dist/(l+1)#+0.1
                size = linewidth
                self.noisemask.append(True if np.random.rand()<self.noise else False) 
                self.linewidths.append(size*2) 
                self.all_lines.append(np.array(get_track_line(radius, begin, end,self.granularity)))
                self.segments.loc[self.total_segments,:] = [begin,end,radius]
        l = layers+1
        len_segment = 2*np.pi/(self.ecl_segments)#+k*l)
        for i in range(self.ecl_segments):#+k*l):                               #making the ecl segments
            self.total_segments+=1
            radius = l
            begin = len_segment*i+(dist/(l+1))
            end = len_segment*(i+1)-(dist)/(l+1)
            size = 3*linewidth/5
            self.linewidths.append(size*1.3)                      
            self.noisemask.append(True if np.random.rand()<self.noise else False) 
            self.all_lines.append(np.array(get_ecl_lines(radius, begin, end, 1,self.granularity)))
            self.segments.loc[self.total_segments,:] = [begin,end,radius]

        self.all_lines=np.array(self.all_lines).transpose((0,2,1))              
        self.particle_masks=[]
        self.tracker_mask=np.full(self.total_segments,False)
        self.lines = [] 
        self.segments.loc[:,"radius"] = self.segments.loc[:,"radius"].astype("float")

    def make_tracker_mask(self,truth_particles=[]):                                #this only has to be called once, makes the mask for all segemnts that have been hit by a true particle or noise. these will be the red segments
        for i in range(len(truth_particles)):
            self.particle_masks.append(self.check_hit(truth_particles[i]))
            self.tracker_mask=np.logical_or(self.tracker_mask,self.particle_masks[i])
        self.tracker_mask=np.logical_or(self.tracker_mask,self.noisemask)

    def get_tracker_collection(self):                                           #returns Linecollection for the Tracker (gray) and the segments that have been hit by the true particles or noise (red)
        #just the tracker
        tracker = self.all_lines
        colors = [self.trackercolor]*tracker[:,0,0].size
        linewidth = (np.array(self.linewidths)*15/self.layers).tolist()
        #hits in the tracker
        hits=self.all_lines[self.tracker_mask,:,:]
        tracker=np.append(tracker,hits,axis=0)
        colors.extend(["red"]*hits[:,0,0].size)
        linewidth.extend([2.5]*hits[:,0,0].size)
        return LineCollection(tracker, color = colors, linewidths = linewidth)

    def get_hits_and_misses(self,particle,particle_index):                      #compares hits of the true particles and hits of the simulated particles
        if(self.ignore_noise == True):
            hits=np.logical_and(self.check_hit(particle),self.particle_masks[particle_index])
            misses=np.logical_and(self.check_hit(particle),np.logical_not(self.particle_masks[particle_index]))
        else:
            hits=np.logical_and(self.check_hit(particle),self.tracker_mask)
            misses=np.logical_and(self.check_hit(particle),np.logical_not(self.tracker_mask))
        return [hits.sum(),misses.sum()]

    def get_arrowangle(self,particle):                                          #returns angle for the arrow (angle of the outer most segment)
        last_hit_segment = (self.all_lines[self.check_hit(particle),:,:])[-1,:,:]
        x,y = last_hit_segment.T
        phi = np.arctan2(x,y)
        return -np.mean(phi)+np.pi/2
    
    def get_hit_lines(self,particles,particle_index):                           #returns x,y arrays for segments that have been hit by the given particles and color list for these segments
        segments=np.empty((0,self.granularity,2))
        colors=[]
        for i in range(len(particles)):
            check_hit=self.check_hit(particles[i])
            segments=np.append(segments,self.all_lines[check_hit,:,:],axis=0)
            if(i == particle_index):
                colors.extend(["blue"]*check_hit.sum())
            else:
                colors.extend(["yellow"]*check_hit.sum())
        return segments,colors

    def check_hit(self, particle):                                              #returns mask for segments that have benn hit by a given particle
        theta=particle.charge*np.pi/2-particle.phi+np.arccos(self.segments.radius/(2*particle.radius))*particle.charge-np.pi/2

        mask = theta < 0
        theta[mask] += 2*np.pi
        mask = theta < 0
        theta[mask] += 2*np.pi
        mask = theta > 2*np.pi
        theta[mask] -= 2*np.pi


        return (theta > self.segments.begin) & (theta < self.segments.end) & (self.segments.radius <= abs(2*particle.radius))
