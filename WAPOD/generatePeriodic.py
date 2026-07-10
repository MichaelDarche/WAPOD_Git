#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 15:31:08 2025

@author: michael
"""

import sys
import re

def lire_fichier_milieu(nom_fichier):
    with open(nom_fichier, "r") as f:
        lignes = [l.strip() for l in f.readlines() if l.strip() != ""]

    nbcell = int(lignes[0].split(":")[1].strip())
    nb_milieux = int(lignes[1].split(":")[1].strip())

    blocs = []
    bloc = []
    
    # On commence à la ligne 2
    for ligne in lignes[2:]:
        if ligne.startswith("milieu"):   # debut d'un nouveau bloc
            if bloc:  
                blocs.append(bloc)
            bloc = [ligne]
        else:
            bloc.append(ligne)

    # Ajouter le dernier bloc
    if bloc:
        blocs.append(bloc)

    return nbcell, nb_milieux, blocs


def ecrire_milieu(nom_fichier, nbcell, nb_milieux, blocs):
    compteur = 0

    with open(nom_fichier, "w") as f:
        f.write(f"nombre de milieux : {nbcell*nb_milieux}\n\n")
        for _ in range(nbcell):
            for bloc in blocs:
                for idx, ligne in enumerate(bloc):
                    if idx == 0:  # ligne "milieu X"
                        f.write(f"milieu {compteur}\n")
                        compteur += 1
                    else:
                        f.write(ligne + "\n")
                f.write("\n")
                
                
   
def lire_fichier_frontieres(nom_fichier):
    with open(nom_fichier, "r") as f:
        lignes = [l.strip() for l in f.readlines() if l.strip() != ""]

    nbcell = int(lignes[0].split(":")[1].strip())
    nb_frontieres = int(lignes[1].split(":")[1].strip())

    blocs = []
    bloc = []

    for ligne in lignes[2:]:
        if ligne.startswith("frontiere"):
            if bloc:
                blocs.append(bloc)
            bloc = [ligne]
        else:
            bloc.append(ligne)

    if bloc:
        blocs.append(bloc)

    return nbcell, nb_frontieres, blocs


def extraire_alpha(ligne):
    """Retourne la valeur de Alpha dans une ligne du type 'Alpha : xxxxx' """
    m = re.search(r"Alpha\s*:\s*([0-9.+Ee-]+)", ligne)
    return float(m.group(1)) if m else None


def ecrire_frontieres_sortie(nom_fichier, nbcell, nb_frontieres, blocs):
    compteur = 1
    
    # Alpha du dernier bloc => taille de la cellule
    tailleCellule = None
    for ligne in blocs[-1]:
        a = extraire_alpha(ligne)
        if a is not None:
            tailleCellule = a
            break

    if tailleCellule is None:
        raise RuntimeError("Impossible de trouver Alpha dans le dernier bloc.")

    with open(nom_fichier, "w") as f:
        f.write(f"nombre de frontieres : {nbcell*nb_frontieres-1}\n\n")
        for n in range(nbcell):
            for bloc in blocs:
                alpha_original = None

                # On recupere l'Alpha original du bloc courant
                for ligne in bloc:
                    a = extraire_alpha(ligne)
                    if a is not None:
                        alpha_original = a
                        break

                for idx, ligne in enumerate(bloc):

                    if idx == 0:
                        f.write(f"frontiere {compteur}\n")
                        compteur += 1

                    elif ligne.startswith("Alpha"):
                        # Nouveau Alpha = Alpha original + n × tailleCellule
                        alpha_nouveau = alpha_original + n * tailleCellule
                        f.write(f"Alpha   :  {alpha_nouveau:.14E}\n")

                    else:
                        f.write(ligne + "\n")

                f.write("\n")             
            
