Tracking widget:

Stellt ein Teilchen Tracker dar. Der Tracker ließt aus einem .h5 File die Daten über mehrere Teilchen aus. Im Tracker werden alle Segmente die die Teilchen überflogen haben rot gekennzeichnet, zusätzlich kommt noch Rauschen dazu. 
Für jedes dieser Teilchen(truthparticles) kann man ein Teilchen simulieren und an die wahren Teilchen anpassen indem man den Impuls und Winkel des Teilchens anpasst. An welches Teilchen man anpassen soll zeigt der Pfeil an. Wie gut die 
Anpassung ist kann man aus der Anzeige "x hits & y misses" auslesen. Diese gibt an wie oft das simulierte Teilchen, die selben Segmente getroffen hat wie das wahre Teilchen, bzw. wie oft nicht. 

Tracking_widget(data_path,B=0.2,layers=15,n_segments = 5, ecl_segments = 30, k = 3, dist = 0.1, noise = 0.05, linewidth = 5, show_truthbutton = False, continuous_update=True, truthvalues=False, ignore_noise=False, trackercolor="gray")

data_path: Pfad für ein .h5 File die Teilchen und Information über deren Impuls, Masse, Energie und ... enthält. Siehe 2part_events für den Aufbau solch einer Datei.
B: B-feld Einstellung für den Tracker am besten zwischen 0.1 und 0.4
layer: Anzahl der Layers im Tracker. Umso mehr Layer man hat, umso besser kann man die Teilchen anpassen. Mehr als 60 wird unübersichtlich.
n_segments: Anzahl in der Segmente im ersten layer
ecl_segments: Anzahl der größeren ECL_segmente im letzten Layer
k: Gibt bei den normalen Segmenten an, wieviel Segmente bei jedem Layer dazu kommen
dist: Gibt distanz zwischen den einzelnen Segmenten an
noise: konstante für Noisefloor im Tracker, 0.1 bedeutet 10% Prozent aller Segmente werden fälschlicherweise Rot gekennzeichnet
linewidth: Linienbreite, am besten bei 5 lassen
show_truthbutton: Mit dieser Funktion auf True wird im Widget ein Button angezeigt mit dem man die Flugbahn des wahren Teilchen Anzeigen kann.
continuous_update: Einstellung ob die Slider im Widget kontinuierlich ausgelesen werden oder nur zu gewissen Punkten. Nur auf False stellen, wenn es Performance Probleme gibt.
truthvalues: Wenn auf True werden direkt die richtigen Werte für Phi und impuls der simulierten Teilchen eingestellt.
ignore_noise: Mit dieser Option auf True wird für die Anzeige x hit & y misses, das simulierte Teilchen NUR mit dem wahren Teilchen verglichen. Auf false werden auch das Rauschen und die anderen Teilchen beachtet.
trackercolor: Farbe des Trackers


ECL widget:
Stellt ein kaloriemeter dar. Ziel hier ist die Energie der einzelnen Teilchen herauszufinden. der ECL ließt aus dem selben .h5 file wie der Tracker die Information darüber, in welchem segment wieviel Energie ist. Mit dem widget
kann man dann jeweils die angezeigten Teilchen umkreisen(lasso selector) und deren Energie auslesen. 

ECL_widget(data_path, noise_rate = 0.05, idx=None)

data_path: Pfad für ein .h5 File die Teilchen.
noise_rate: konstante für den Noisefloor
idx: hiermit kann man nur ein einzelnen Teilchen zeigen lassen, sollte normal auf None bleiben.

Matching widget:
Hier werden jetzt die gesammelten Ergebnisse aus dem Tracking_widget und ECL_widget dargestellt. Man kann vergleichen und nachvollziehen was für Teilchen man gefunden hat.

MatchingWidget(ew, tw, cheat_mode=True, cheating_threshhold = 1e-2)

ew: ECL_widget object
tw: Tracking_widget object
cheat_mode: Da bei Elektronen es oft kaum möglich ist sie richtig zu fitten kann das Matching widget die Werte anpassen, falls die gefundenen Werte nahe genug an den echten sind(wenn auf True)
cheating_threshhold: Threshhold für differenz ab dem die Werte angepasst werden. (nur relevant für cheat_mode=True)