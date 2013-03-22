Finite Element Method (MECA1120)
================================

Tous les documents et codes sources relatifs aux probl�mes et projets du cours d'�l�ments finis.

Homework1: Int�gration sur un quadrilat�re
------------------------------------------

Programme simple permettant d'int�grer une fonction 2D quelconque sur un quadrilat�re quelconque utilisant la r�gle d'int�gration de Gauss-Legendre sur 4 points.

Homework2: Volume de la mer M�diterran�e
--------------------------------------------------

Programme permettant l'int�gration de la bathym�trie de la mer M�diterran�e sur un maillage non-structur� de quadrilat�res quelconque en utilisant une r�gle d'int�gration quelconque (� la base, Gauss-Legendre). Le mesh provient du NOAA_.

Homework3: Calcul de la longueur du littoral de la mer
-------------------------------------------------------

Programme permettant le gestion, le tri, et la s�lection des segments dans un meshGrid. Une fois les segments ext�rieurs trouv�s et s�lectionn�s, la d�termination de la longueur du littoral est possible. 

Homework4: Premier probl�me d'�l�ments fini
-------------------------------------------

Mise en pratique de la m�thode des �l�ments finis pour l'�quation de poisson pour diff�rents mesh.

Homework5: Les solveurs lin�aires
---------------------------------

Am�liorer les performances de notre de notre premier programme en r�num�rotant les noeuds et en appliquant des soldeurs bandes

Homework5: Augmenter l'ordre
----------------------------

Jusque l�, nous avons utiliser des �l�ments (bi)lin�aires. Nous allons augmenter la pr�cisions des m�thodes en augmentant l'ordre. Il s'agit donc de transformer nos �l�ments (bi)lin�aires en �l�ments (bi)quadratique ou (bi)cubique en leur rajoutant des noeuds, et en calculant les fonctions de formes correspondantes.


.. _NOAA: http://www.ngdc.noaa.gov/mgg/gdas/gd_designagrid.html