Finite Element Method (MECA1120)
================================

Tous les documents et codes sources relatifs aux problèmes et projets du cours d'éléments finis.

Homework1: Intégration sur un quadrilatère
------------------------------------------

Programme simple permettant d'intégrer une fonction 2D quelconque sur un quadrilatère quelconque utilisant la règle d'intégration de Gauss-Legendre sur 4 points.

Homework2: Volume de la mer Méditerranée
--------------------------------------------------

Programme permettant l'intégration de la bathymétrie de la mer Méditerranée sur un maillage non-structuré de quadrilatères quelconque en utilisant une règle d'intégration quelconque (à la base, Gauss-Legendre). Le mesh provient du NOAA_.

Homework3: Calcul de la longueur du littoral de la mer
-------------------------------------------------------

Programme permettant le gestion, le tri, et la sélection des segments dans un meshGrid. Une fois les segments extérieurs trouvés et sélectionnés, la détermination de la longueur du littoral est possible. 

Homework4: Premier problème d'éléments fini
-------------------------------------------

Mise en pratique de la méthode des éléments finis pour l'équation de poisson pour différents mesh.

Homework5: Les solveurs linéaires
---------------------------------

Améliorer les performances de notre de notre premier programme en rénumérotant les noeuds et en appliquant des soldeurs bandes

Homework5: Augmenter l'ordre
----------------------------

Jusque là, nous avons utiliser des éléments (bi)linéaires. Nous allons augmenter la précisions des méthodes en augmentant l'ordre. Il s'agit donc de transformer nos éléments (bi)linéaires en éléments (bi)quadratique ou (bi)cubique en leur rajoutant des noeuds, et en calculant les fonctions de formes correspondantes.


.. _NOAA: http://www.ngdc.noaa.gov/mgg/gdas/gd_designagrid.html