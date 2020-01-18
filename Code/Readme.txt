Ce r�pertoire contient les outils informatiques du groupe US pour l'estimation de mouvement programm� en langage MatLab..
ll contient trois repertoires :
	
	 -Method

		- contient des m�thodes d'estimation du mouvement (qui utilisent �ventuellement des estimateurs se trouvant dans le r�pertoire LocalEstimation)
		- des images tests sont disponibles pour chaque methode
		- pour estimer il faut lancer les fichiers lancer.m, disponibles dans le r�pertoire de chaque m�thode

	-LocalEstimation 

		- contient des estimateurs locaux des d�calages (utilisables avec des m�thodes de type block matching)
		- les fonctions contiennent des exemples d'execution
		
	
	-Tools
 
		- proc�dures d'exploitation des cartes des d�placements estim�s :
				- d�rivation pour calculer des cartes de d�formations.
				- recalage des images en tenant compte du d�placement estim�s (pour �valuer l'erreur d'estimation par exemple)
		- proc�dures pour la visualisation des champs de mouvement estim�s




Le suivi de ces outils est assur� pas Fran�ois Duboeuf  (6208) et Adrian Basarab (8781)

Document mis a jour le 7 D�cembre 2007