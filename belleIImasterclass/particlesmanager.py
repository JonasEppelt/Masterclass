import pandas as pd
import numpy as np
from copy import deepcopy

class ParticlesManager:
    '''
    class to manage particle measurements across multiple widgets
    '''
    def __init__(self, path: str,Dark_matter_particles=0) -> None:
        #alle Darkmatterteilchen sind im darkmatter df
        #alle sichtbaren sind im _df
        #immer die ersten die ersten Teicleh aus dem h5 file sind die dark_matter teichen        
        self._path = path
        self._df = pd.read_csv(path)
        self.total_n_particles=len(self._df)
        self._df["charge"] = 0# * (-1)
        for i in range(len(self._df)):
            if self._df.loc[i,"pdg"]!=22 and self._df.loc[i,"pdg"]!=111: #gamma's und pi0's haben Ladung 0
                self._df.at[i,"charge"]=np.sign(self._df.loc[i,"pdg"])
            else:
                self._df.at[i,"charge"]=0
        self._df["tracker_pt"] = 0
        self._df["tracker_phi"] = 0
        self._df["tracker_charge"] = 0
        self._df["ecl_energy"] = 0
        self._df["ID_mass"] = 0
        self._df["klm_detect"] = False
        self._df["patches"] = [[]]*len(self._df)
        self._df["patch_edgecolors"] = [np.array([])]*len(self._df)
        self._df["patch_facecolors"] = [np.array([])]*len(self._df)

        self.dark_matter_df=self._df.copy(deep=True)   
        self._df=self._df.drop(labels=np.arange(0,Dark_matter_particles),axis=0)
        self.dark_matter_df=self.dark_matter_df.drop(labels=np.arange(Dark_matter_particles,self.total_n_particles),axis=0)

        self.dark_matter_df=self.dark_matter_df.drop(labels=["tracker_pt","tracker_phi","tracker_charge","ecl_energy","ID_mass",
                                                             "klm_detect","patches","patch_edgecolors","patch_facecolors"],axis=1)
        self.dark_matter_df["ew_charge"]=0
        self.dark_matter_df["ew_px"]=0
        self.dark_matter_df["ew_py"]=0
        self.dark_matter_df["ew_pz"]=0
        self.dark_matter_df["ew_energy"]=0
        self.dark_matter_df["ew_mass"]=0

        self._df=self._df.set_index(np.arange(len(self._df)))
        self.dark_matter_df=self.dark_matter_df.set_index(np.arange(len(self.dark_matter_df)))

    @property
    def index(self) -> pd.Index:
        return self._df.index
    @property
    def n_particles(self) -> int:
        return len(self._df)
    @property
    def crystal_column_names(self):
        return [str(i) for i in range(0,8736)]

    def save_measurements_as_csv(self):
        
        measurements=pd.DataFrame(index=self.index,columns=["pt","phi","pz","klm_detect","energy","charge","missing","mass"])

        measurements.loc[self.index,"missing"]=False
        measurements.loc[self.index,"pt"]=self._df["tracker_pt"]
        measurements.loc[self.index,"phi"]=self._df["tracker_phi"]
        measurements.loc[self.index,"pz"]=self._df["pz"]
        measurements.loc[self.index,"klm_detect"]=self._df["klm_detect"]
        measurements.loc[self.index,"mass"]=0
        measurements.loc[self.index,"energy"]=self._df["ecl_energy"]
        measurements.loc[self.index,"charge"]=self._df["tracker_charge"]

        measurements.at[self.n_particles,"missing"]=True
        measurements.at[self.n_particles,"pt"]=np.sqrt(self.dark_matter_df.loc[0,"ew_px"]**2+self.dark_matter_df.loc[0,"ew_py"]**2)
        measurements.at[self.n_particles,"phi"]=np.arctan2(self.dark_matter_df.loc[0,"ew_py"],self.dark_matter_df.loc[0,"ew_px"])
        measurements.at[self.n_particles,"pz"]=self.dark_matter_df.loc[0,"ew_pz"]
        measurements.at[self.n_particles,"klm_detect"]=False
        measurements.at[self.n_particles,"mass"]=self.dark_matter_df.loc[0,"ew_mass"]
        measurements.at[self.n_particles,"energy"]=self.dark_matter_df.loc[0,"ew_energy"]
        measurements.at[self.n_particles,"charge"]=self.dark_matter_df.loc[0,"ew_charge"]

        measurements.to_csv(path_or_buf="Ergebnisse.csv")

    def load_measurements_from_csv(self,path="Ergebnisse.csv"):
        
        measurements= pd.read_csv(path,index_col=0)

        self._df["tracker_pt"]     = measurements.loc[self.index,"pt"]
        self._df["tracker_phi"]    = measurements.loc[self.index,"phi"]
        self._df["tracker_charge"] = measurements.loc[self.index,"charge"]
        self._df["ecl_energy"]     = measurements.loc[self.index,"energy"]
        self._df["klm_detect"]     = measurements.loc[self.index,"klm_detect"]

        self.dark_matter_df.at[0,"ew_px"]     =measurements.loc[self.n_particles,"pt"]*np.cos(measurements.loc[self.n_particles,"phi"])
        self.dark_matter_df.at[0,"ew_py"]     =measurements.loc[self.n_particles,"pt"]*np.sin(measurements.loc[self.n_particles,"phi"])
        self.dark_matter_df.at[0,"ew_pz"]     =measurements.loc[self.n_particles,"pz"]
        self.dark_matter_df.at[0,"ew_energy"] =measurements.loc[self.n_particles,"energy"]
        self.dark_matter_df.at[0,"ew_mass"]   =measurements.loc[self.n_particles,"mass"]
        self.dark_matter_df.at[0,"ew_charge"] =measurements.loc[self.n_particles,"charge"]        

    def __getitem__(self,i) -> pd.Series:
        return self._df.loc[i,:]
    def __len__(self) -> int:
        return len(self._df)
    def tracker_measurement(self, index, pt, phi, charge) -> None:
        self._df.at[index, "tracker_pt"] = pt
        self._df.at[index, "tracker_phi"] = phi
        self._df.at[index, "tracker_charge"] = charge  

    def energy_measurement(self, index, energy) -> None:
        self._df.at[index, "ecl_energy"] = energy

    def Klm_measurement(self, index, hit) -> None:
        self._df.at[index, "klm_detect"] = hit

    def ecal_patches(self, index, patches, facecolors, edgecolors) -> None:
        self._df.at[index,"patches"]=patches
        self._df.at[index,"patch_edgecolors"] = edgecolors
        self._df.at[index,"patch_facecolors"] = facecolors
        

    def missing_particle_measurement(self,index,px,py,pz,energy,mass,charge) -> None:
        if len(self.dark_matter_df) <= index:
            self.dark_matter_df.loc[index]=["0"]*len(self.dark_matter_df.columns)
        self.dark_matter_df.at[index,"ew_px"]=px
        self.dark_matter_df.at[index,"ew_py"]=py
        self.dark_matter_df.at[index,"ew_pz"]=pz
        self.dark_matter_df.at[index,"ew_energy"]=energy
        self.dark_matter_df.at[index,"ew_mass"]=mass
        self.dark_matter_df.at[index,"ew_charge"]=charge

    def get_crystall_content(self, n_particle):
        return np.clip(self._df.loc[n_particle,self.crystal_column_names].to_numpy(),0,100000)
    