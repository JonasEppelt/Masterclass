import pandas as pd
import numpy as np


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
    
    @property
    def index(self) -> pd.Index:
        return self._df.index
    @property
    def n_particles(self) -> int:
        return len(self._df)
    @property
    def crystal_column_names(self) -> list[str]:
        return [str(i) for i in range(0,8736)]

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

    def get_crystall_content(self, n_particle):
        return np.clip(self._df.iloc[n_particle][self.crystal_column_names].to_numpy(),0,100000)
    