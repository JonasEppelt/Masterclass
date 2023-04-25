# Belle II Masterclass

by Jonas Eppelt, Filipp Gostner, Isabel Haide, Torben Ferber,

Basierend auf Jupyter-Notebooks werden die grundlegende Prinzipien im Arbeiten mit Teilchendetektoren wie Belle~II vermittelt.
Mit interaktiven Widgets lernt man wie die Ladung, der Impuls und die Energie von Teilchen bestimmt wird.
Diese Informationen werden kombiniert um verschiedene Teilchen zu identifizieren.
Dadurch können Ereignisse interpretiert werden und nach neuen Teilchen wie Dunkler Materie gesucht werden.

## Die Masterclass nutzen
Clone dieses Repository mit
```bash
git clone https://git.scc.kit.edu/ukjyf/belleii_masterclass.git
```
Diese Masterclass nutzt python3 und benötigt zusätzliche packete die folgendermasen installiert werden können:
```bash
pip install jupyterlab numpy pandas tables matplolib ipympl jupyterlab-git jupyter-lsp ipywidgets=8.0.4
```
Um die Notebooks zu starten, starte jupyter lab
```
jupyter lab
```
und wähle eines der Notebooks `Introduction.ipynb`, `Event_Interpretation.ipynb` oder `Fehlen_hier_Teilchen.ipynb` aus.


## Generieren von Ereignissen
Um neue Ereignisse zu simulieren, wird die Belle~II software `basf2` und b2luigi benötigt.
---

## development install
Im Verzeichnis ausführen:
```bash
python3 -m pip install -e .
```

## Tracking widget

TrackingWidget(particles_manager,...)

Stellt ein Teilchen Tracker dar. Der Tracker ließt aus dem particles_manager die Daten über mehrere Teilchen aus. Im Tracker werden alle Segmente die die Teilchen überflogen haben rot gekennzeichnet, zusätzlich kommt noch Rauschen dazu. 
Für jedes dieser Teilchen(truthparticles) kann man ein Teilchen simulieren und an die wahren Teilchen anpassen indem man den Impuls und Winkel des Teilchens variiert. An welches Teilchen man anpassen soll zeigt der Pfeil an.

### Funktionen:

TrackingWidget.show():
Zeigt das tracking Widget an.

### Parameter zum Tracker:

particles_manager: 
particles manager klasse, welche nötig ist um die events zu laden und die Ergebnisse fest zu halten.<br>

B(=None,float): 
B-feld Einstellung für den Tracker, bei None wird das B-feld automatisch angepasst zu der Snzahl von layers gewählt.<br>

layer(=50,int): 
Anzahl der Layers im Tracker. Umso mehr Layer man hat, umso besser kann man die Teilchen anpassen. Mehr als 60 wird unübersichtlich.<br>

n_segments(=8,int): 
Anzahl der Segmente im ersten layer.<br>

add_segments(=3,int): 
Gibt an, wieviel Segmente bei jedem Layer dazu kommen (layer i hat (n_segments+k*i) Segmente).<br>

noise_ratio(=0.1,float): 
Konstante für Noisefloor im Tracker, 0.1 bedeutet 10% Prozent aller Segmente werden fälschlicherweise Rot gekennzeichnet.<br>

linewidth(=2.5,float): 
Linienbreite, am besten bei 2.5 lassen.<br>

continuous_update(=True,bool): 
Einstellung ob die Slider im Widget kontinuierlich ausgelesen werden oder nur zu gewissen Punkten. Nur auf False stellen, wenn es Performance Probleme gibt.<br>

truthvalues(=False,bool): 
Wenn auf True werden direkt die richtigen Werte für Phi und impuls der simulierten Teilchen eingestellt.<br>

granularity(=100,int):
Sollte einfach auf 100 bleiben.

show_hitcounter(=False,bool):
Zeigt ein Hitcounter für die simulierten Teilchen an.

## ECL widget:
Stellt ein Kaloriemeter dar. Ziel hier ist die Energie der einzelnen Teilchen herauszufinden. der ECL ließt aus dem selben particles manager wie der Tracker die Information darüber, in welchem Segment wieviel Energie ist. In dem Widget
kann man dann jeweils die angezeigten Teilchen umkreisen(lasso selector) und deren Energie auslesen. 

### Funktionen:

ECLWidget.show():
Zeigt das ECL Widget an.

### Parameter zum Widget:

particles_manager: 
particles manager klasse, welche nötig ist um die events zu laden und die Ergebnisse fest zu halten.<br>

noise_rate(=0.2,float): 
Konstante für den Noisefloor. <br>

true_particles:(=False,bool)
... <br>

## Matching Widget
Hier werden jetzt die gesammelten Ergebnisse aus dem Tracking_widget und ECL_widget dargestellt. Man kann vergleichen und nachvollziehen was für Teilchen man gefunden hat.

### Funktionen:

MatchingWidget.show():
Zeigt das matching Widget an.

### Parameter zum Widget:

ew: 
ECL_widget object <br>

tw: 
Tracking_widget object <br>

cheat_mode(=True,bool): 
Da bei Elektronen es oft kaum möglich ist sie richtig zu fitten kann das Matchingwidget die Werte anpassen, falls die gefundenen Werte nahe genug an den echten sind(wenn auf True). <br>

cheating_threshhold(=1e-2,float): 
Threshhold für differenz ab dem die Werte angepasst werden. (nur relevant für cheat_mode=True) <br>