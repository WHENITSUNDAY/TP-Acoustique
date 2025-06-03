#!/bin/bash

echo "Suppression des fichiers précédents..."

find ./data -type f -delete
find ./images -type f -delete

echo "Compilation en cours..."
make
if [ $? -ne 0 ]; then
    echo "Erreur lors de la compilation"
    exit 1
fi

echo "Exécution du programme..."
./run
if [ $? -ne 0 ]; then
    echo "Erreur lors de l'exécution du programme"
    exit 1
fi

echo "Génération des graphiques..."
gnuplot plot/ondes_2D.gplot
if [ $? -ne 0 ]; then
    echo "Erreur lors de la génération des graphiques"
    exit 1
fi

echo "Conversion des images en GIF..."

convert -delay 20 -loop 0 images/P_*.png simulation.gif
if [ $? -eq 0 ]; then
    echo "GIF créé avec succès: simulation.gif"
else
    echo "Erreur lors de la création du GIF"
    exit 1
fi