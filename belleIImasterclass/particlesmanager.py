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
        self._df["ecl_E"] = 0
        self._df["tracker_charge"] = 0
        self._df["ecl_energy"] = 0
    
    @property
    def index(self) -> pd.Index:
        return self._df.index
    @property
    def n_particles(self) -> int:
        return len(self._df)

    def __getitem__(self,i) -> pd.Series:
        return self._df.loc[i,:]
    def __len__(self) -> int:
        return len(self._df)
    def tracker_measurement(self, index, pt, phi, charge) -> None:
        self._df.loc[index, "tracker_pt"] = pt
        self._df.loc[index, "tracker_phi"] = phi
        self._df.loc[index, "tracker_charge"] = charge  