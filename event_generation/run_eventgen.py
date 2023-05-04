try:
    from simulation import add_simulation
    from reconstruction import add_reconstruction
except ImportError:
    # Output expected ImportErrors.
    print("Import Error")
import basf2
from ROOT import Belle2
import pandas as pd
import numpy as np

class ECLInfoExtractor(basf2.Module):
    def __init__(self):
        super().__init__()

    def initialize(self):
        self.mc_particles = Belle2.PyStoreArray('MCParticles')
        self.ecl_cal_digits = Belle2.PyStoreArray('ECLCalDigits')
        self.b2_ecl_clusters = Belle2.PyStoreArray('ECLClusters')
        self.ecllocalmaxima = Belle2.PyStoreArray('ECLLocalMaximums')
        self.eventinfo = Belle2.PyStoreObj('EventMetaData')

    def event(self) -> None:
        uni_event_id = f"{self.eventinfo.getExperiment()}_{self.eventinfo.getRun()}_{self.eventinfo.getEvent()}"
        df = pd.DataFrame()

        # get mc particle information
        for mc_particle_id, mc_particle in enumerate(self.mc_particles):
            if mc_particle.getPDG() != 300553:
                df[mc_particle_id] = pd.Series()
                vector = mc_particle.get4Vector()
                df.loc[mc_particle_id, "px"] = vector.Px()
                df.loc[mc_particle_id, "py"] = vector.Py()
                df.loc[mc_particle_id, "pz"] = vector.Pz()
                df.loc[mc_particle_id, "energy"] = vector.P()
                df.loc[mc_particle_id, "p"] = vector.E()
                df.loc[mc_particle_id, "pt"] = vector.Pt()
                df.loc[mc_particle_id, "theta"] = vector.Theta()
                df.loc[mc_particle_id, "phi"] = vector.Phi()
                df.loc[mc_particle_id, "pdg"] = mc_particle.getPDG()
                df.loc[mc_particle_id, "mass"] = mc_particle.getMass()
                df.loc[mc_particle_id, "charge"] = mc_particle.getCharge()
                
                df.loc[mc_particle_id, [str(i)  for i in range(8735+1)]] = np.zeros(8735+1)
            # get mc depositions
                mc_relations = mc_particle.getRelationsWith("ECLCalDigits")
                for mc_relations_id in range(mc_relations.size()):
                    mc_energy = mc_relations.weight(mc_relations_id)
                    cell_id = mc_relations.object(mc_relations_id).getCellId()
                    df.loc[mc_particle_id, str(cell_id)] = mc_energy

        df.to_csv(f"{uni_event_id}.csv")

    def terminate(self):
        pass

if __name__ == "__main__":
    main = basf2.create_path()
    main.add_module('EventInfoSetter',
            evtNumList=100)
    
    decfile = basf2.find_file('./decfile')
    main.add_module('EvtGenInput', userDECFile=decfile)
    

    main.add_module('Gearbox')
    main.add_module('Geometry', useDB=True)

    add_simulation(path=main)

    add_reconstruction(path=main)

    ecl_output = ECLInfoExtractor()

    main.add_module(ecl_output)
    main.add_module('Progress')
    basf2.process(main)