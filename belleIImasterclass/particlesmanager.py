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
        self._df["colors"] = [np.array([])]*len(self._df)
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

    def get_measurements_csv(self):
        self._df.to_csv(path_or_buf="Ergebnisse.csv",columns=["tracker_pt","tracker_phi","tracker_charge","ecl_energy","klm_detect"])
        self.missing_df.to_csv(path_or_buf="fehlendes_Teilchen.csv",columns=["px","py","pz","p","energy","mass","charge"])

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

    def ecal_patches(self, index, patches, colors) -> None:
        self._df.at[index,"patches"]=patches
        self._df.at[index,"colors"] = colors

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
    