from belleIImasterclass.particlesmanager import ParticlesManager
from ipywidgets import Output, Accordion, Text, HBox, VBox


class ECLWidget:
    def __init__(self,particles_manager: ParticlesManager, noise_rate = 0.05) -> None:
        self._particles_manager = particles_manager
        self._noise_rate = noise_rate


        self._out = Output