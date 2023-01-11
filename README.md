# Belle II Masterclass

by Torben Ferber, Isabel Haide, Jonas Eppelt, Filipp Gostner

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
pip install jupyterlab numpy pandas tables matplolib ipympl jupyterlab-git jupyter-lsp ipywidgets=7.7.2
```
Um die Notebooks zu starten, starte jupyter lab
```
jupyter lab
```
und wähle eines der Notebooks `Introduction.ipynb`, `Event_Interpretation.ipynb` oder `Fehlen_hier_Teilchen.ipynb` aus.


## Generieren von Ereignissen
Um neue Ereignisse zu simulieren, wird die Belle~II software `basf2` und b2luigi benötigt.
---

## Tracking widget

Stellt ein Teilchen Tracker dar. Der Tracker ließt aus einem .h5 File die Daten über mehrere Teilchen aus. Im Tracker werden alle Segmente die die Teilchen überflogen haben rot gekennzeichnet, zusätzlich kommt noch Rauschen dazu. 
Für jedes dieser Teilchen(truthparticles) kann man ein Teilchen simulieren und an die wahren Teilchen anpassen indem man den Impuls und Winkel des Teilchens variiert. An welches Teilchen man anpassen soll zeigt der Pfeil an. Wie gut die 
Anpassung ist kann man aus der Anzeige "x hits & y misses" auslesen. Diese gibt an wie oft das simulierte Teilchen, die selben Segmente getroffen hat wie das wahre Teilchen, bzw. wie oft nicht. 

### Funktionen:

TrackingWidget.show():
Zeigt das tracking Widget an.

### Parameter zum Tracker:

data_path(.h5): 
Pfad für ein .h5 File die Teilchen und Information über deren Impuls, Masse, Energie und ... enthält. Siehe 2part_events für den Aufbau solch einer Datei.<br>

B(=0.2,float): 
B-feld Einstellung für den Tracker am besten zwischen 0.1 und 0.4<br>

layer(=15,int): 
Anzahl der Layers im Tracker. Umso mehr Layer man hat, umso besser kann man die Teilchen anpassen. Mehr als 60 wird unübersichtlich.<br>

n_segments(=8,int): 
Anzahl der Segmente im ersten layer.<br>

ecl_segments(=30,int): 
Anzahl der größeren ECL_segmente im letzten Layer.<br>

k(=3,int): 
Gibt bei den normalen Segmenten an, wieviel Segmente bei jedem Layer dazu kommen (layer i hat (n_segments+k*i) Segmente).<br>

dist(=0.1,float): 
Gibt die Distanz zwischen den einzelnen Segmenten an.<br>

noise(=0.1,float): 
Konstante für Noisefloor im Tracker, 0.1 bedeutet 10% Prozent aller Segmente werden fälschlicherweise Rot gekennzeichnet.<br>

linewidth(=5,float): 
Linienbreite, am besten bei 5 lassen.<br>

show_truthbutton(=False,bool): 
Mit dieser Funktion auf True wird im Widget ein Button angezeigt mit dem man die Flugbahn des wahren Teilchen Anzeigen kann.<br>

continuous_update(=True,bool): 
Einstellung ob die Slider im Widget kontinuierlich ausgelesen werden oder nur zu gewissen Punkten. Nur auf False stellen, wenn es Performance Probleme gibt.<br>

truthvalues(=False,bool): 
Wenn auf True werden direkt die richtigen Werte für Phi und impuls der simulierten Teilchen eingestellt.<br>

ignore_noise(=False,bool): 
Mit dieser Option auf True wird für die Anzeige x hit & y misses, das simulierte Teilchen NUR mit dem wahren Teilchen verglichen. Auf False werden auch das Rauschen und die anderen Teilchen beachtet.<br>

trackercolor(="gray",string(color)): 
Farbe des Trackers.<br>

### propertys:

get_fitted_particles:
Returnt ein Dataframe mit Gesamtimpuls,phi,Ladung,Radius und der Impulskomponenten.

## ECL widget:
Stellt ein Kaloriemeter dar. Ziel hier ist die Energie der einzelnen Teilchen herauszufinden. der ECL ließt aus dem selben .h5 file wie der Tracker die Information darüber, in welchem Segment wieviel Energie ist. In dem Widget
kann man dann jeweils die angezeigten Teilchen umkreisen(lasso selector) und deren Energie auslesen. 

### Funktionen:

ECLWidget.show():
Zeigt das ECL Widget an.

### Parameter zum Widget:

data_path(.h5): 
Pfad für ein .h5 File mit den Teilchen. <br>

noise_rate(=0.05,float): 
Konstante für den Noisefloor. <br>

idx(=None,int or None): 
hiermit kann man nur ein einzelnen Teilchen zeigen lassen, sollte normal auf None bleiben. <br>

### propertys:

get_particles_energy:
Returnt ein Dataframe mit den Energien die man harausgefunden hat. <br>

get_particles__radius:
Returnt ein Dataframe mit den Radien der Teilchen. <br>

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