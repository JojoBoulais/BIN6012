Salut mon beau Marco!

J'utilise le fichier project.sh pour lancer des scripts avec slurm
Sans doute changer ou retirer la valeur de --account

`sbatch project.sh`

De manière générale j'utilise le script `preprocess.py` pour preprocess les fichiers directement download de CellxGene
J'input ensuite ce fichier dans `grn.py`

Les autres fichier .py c'est plus des tests pour le moment.
