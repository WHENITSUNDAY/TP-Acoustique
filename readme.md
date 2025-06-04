# Simulation de la propagation d'ondes acoustiques dans un milieu 2D

Ce code a pour but de simuler la propagation d'ondes acoustiques dans un milieu bidimensionnel entre deux dispositifs piézoélectriques, l'un émetteur et l'autre récepteur.
On considère un premier cas où le milieu est seulement composé d'où, puis un second cas où on ajoute une plaque d'aluminium entre les deux dispositifs piézoélectriques.
On peut donc visualiser la propagation complexe de l'onde acoustique dans l'eau, puis dans l'eau et l'aluminium (ce que l'on ne peut pas faire pendant le TP ! Car on a seulement accès au signal reçu)
En parallèle, on modélise notamment le signal électrique reçu par le récepteur piézoélectrique, qui est proportionnel à la pression acoustique moyenne dans la zone du récepteur.
Cela permet donc de comparer (visuellement) avec les résultats expérimentaux obtenus pendant le TP au Fablab.
Le but est donc de montrer aux élèves que la simulation numérique est un outil puissant pour comprendre et visualiser des phénomènes physiques complexes, 
parfois difficilement accessibles par l'expérimentation directe, et ce avec un simple code Python d'une centaine de lignes.

## Résultats de la simulation

### Simulation eau seule
<video src="onde_acoustique_eau_alu.mp4" width="80%" controls></video>

### Simulation eau + aluminium
<video src="onde_acoustique_eau_alu.mp4" width="80%" controls></video>