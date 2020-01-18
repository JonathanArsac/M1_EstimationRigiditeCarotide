Ce répertoire contient les outils informatiques du groupe US pour l'estimation de mouvement programmé en langage MatLab..
ll contient trois repertoires :
	
	 -Method

		- contient des méthodes d'estimation du mouvement (qui utilisent éventuellement des estimateurs se trouvant dans le répertoire LocalEstimation)
		- des images tests sont disponibles pour chaque methode
		- pour estimer il faut lancer les fichiers lancer.m, disponibles dans le répertoire de chaque méthode

	-LocalEstimation 

		- contient des estimateurs locaux des décalages (utilisables avec des méthodes de type block matching)
		- les fonctions contiennent des exemples d'execution
		
	
	-Tools
 
		- procédures d'exploitation des cartes des déplacements estimés :
				- dérivation pour calculer des cartes de déformations.
				- recalage des images en tenant compte du déplacement estimés (pour évaluer l'erreur d'estimation par exemple)
		- procédures pour la visualisation des champs de mouvement estimés




Le suivi de ces outils est assuré pas François Duboeuf  (6208) et Adrian Basarab (8781)

Document mis a jour le 7 Décembre 2007