import pandas as pd
import numpy as np
from copy import deepcopy

class ParticlesManager:
    '''
    class to manage particle measurements across multiple widgets
    '''
    def __init__(self, path: str) -> None:
        self._path = path
        self._df = pd.read_hdf(path)
        self._df["charge"] = np.sign(self._df["pdg"]) * (-1)
        self._df["tracker_pt"] = 0
        self._df["tracker_phi"] = 0
        self._df["tracker_charge"] = 0
        self._df["ecl_energy"] = 0
        self._df["klm_detect"] = False
        self._df["patches"] = [[]]*len(self._df)
        self._df["patch_edgecolors"] = [np.array([])]*len(self._df)
        self._df["patch_facecolors"] = [np.array([])]*len(self._df)
        self.missing_df=pd.DataFrame(data=[[0,0,0,0,0,0,0]],columns=["px","py","pz","p","energy","mass","charge"])
    
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
        measurements.at[self.n_particles,"pt"]=np.sqrt(self.missing_df.loc[0,"px"]**2+self.missing_df.loc[0,"py"]**2)
        measurements.at[self.n_particles,"phi"]=np.arctan2(self.missing_df.loc[0,"py"],self.missing_df.loc[0,"px"])
        measurements.at[self.n_particles,"pz"]=self.missing_df.loc[0,"pz"]
        measurements.at[self.n_particles,"klm_detect"]=False
        measurements.at[self.n_particles,"mass"]=self.missing_df.loc[0,"mass"]
        measurements.at[self.n_particles,"energy"]=self.missing_df.loc[0,"energy"]
        measurements.at[self.n_particles,"charge"]=self.missing_df.loc[0,"charge"]

        measurements.to_csv(path_or_buf="Ergebnisse.csv")

    def load_measurements_from_csv(self,path="Ergebnisse.csv"):
        
        measurements= pd.read_csv(path)

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
        

    def missing_particle_measurement(self,px,py,pz,energy,mass,charge) -> None:
        self.missing_df.at[0,"px"]=px
        self.missing_df.at[0,"py"]=py
        self.missing_df.at[0,"pz"]=pz
        self.missing_df.at[0,"p"]=np.sqrt(px**2+py**2+pz**2)
        self.missing_df.at[0,"energy"]=energy
        self.missing_df.at[0,"mass"]=mass
        self.missing_df.at[0,"charge"]=charge

    def get_crystall_content(self, n_particle):
        return np.clip(self._df.loc[n_particle,self.crystal_column_names].to_numpy(),0,100000)
    