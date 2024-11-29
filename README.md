# Pipeline automatisé pour l'appel de variants

## Description
Ce pipeline est un outil automatisé qui permet de passer des données brutes de séquençage (FASTQ) à un fichier de variants annoté (VCF). Il inclut des étapes pour le contrôle de qualité, le trimming, l'alignement, le traitement des fichiers alignés, l'appel de variants et, optionnellement, l'annotation fonctionnelle des variants.

## Fonctionnalités
- Démultiplexage des données avec Sabre.
- Analyse de la qualité des lectures avant et après trimming avec FastQC.
- Suppression des adaptateurs et des bases de faible qualité avec Cutadapt.
- Alignement des lectures sur un génome de référence avec BWA.
- Conversion des fichiers SAM en BAM triés et indexés avec SAMtools.
- Appel de SNPs et indels avec SAMtools et BCFtools.
- Optionnel : Annotation fonctionnelle des variants avec SnpEff.

## Fichiers d'entrée
- **Fichier FASTQ** : Données brutes de séquençage (`data/FC20150701_1.fq.gz`).
- **Fichier de barcodes** : Pour le démultiplexage (`data/FC20150701_1.txt`).
- **Génome de référence** : Fichier FASTA du génome de référence (`refgenome/Gmax_275_v2.0.fa`).

## Fichiers de sortie
- **Fichiers FastQC** : Rapports de qualité avant et après trimming.
- **Fichiers BAM** : Alignements triés et indexés.
- **Fichiers VCF** : Variants détectés.
- **Fichier annoté (optionnel)** : Variants avec annotations fonctionnelles.

## Dépendances
Les outils suivants doivent être installés et disponibles dans votre `$PATH` :
- [Sabre](https://github.com/najoshi/sabre)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Cutadapt](https://cutadapt.readthedocs.io/)
- [BWA](http://bio-bwa.sourceforge.net/)
- [SAMtools](http://www.htslib.org/)
- [BCFtools](http://www.htslib.org/)
- (Optionnel) [SnpEff](http://snpeff.sourceforge.net/)

## Installation
1. Clonez ce dépôt GitHub :
   ```bash
   git clone https://github.com/votre-repo/variant-calling-pipeline.git
   cd variant-calling-pipeline
2. Assurez-vous que toutes les dépendances sont installées et configurées correctement.
3. Placez vos fichiers d'entrée dans le dossier data/.

## Utilisation
Pour exécuter le pipeline, utilisez la commande suivante :
```bash
bash pipeline.sh
```

## Explication :
- La ligne de début de la boîte de code est marquée par trois accents graves suivis d'une langue (par exemple `bash` pour la coloration syntaxique).
- La ligne de fin de la boîte de code est simplement trois accents graves (\`\`\`) sans rien après.




## Structure du dépôt 
variant-calling-pipeline/
│
├── data/                   # Fichiers d'entrée (FASTQ, barcodes, génome de référence)
├── results/                # Résultats générés par le pipeline
├── scripts/                # Scripts auxiliaires (e.g., sabre.sh)
├── logs/                   # Fichiers journaux pour chaque étape
├── README.md               # Documentation du pipeline
└── pipeline.sh             # Script principal du pipeline

## Dépannage 
1. Vérifiez les fichiers journaux dans le dossier logs/ pour identifier les erreurs.
2. Assurez-vous que vos fichiers d'entrée sont correctement formatés.
3. Testez chaque étape individuellement si nécessaire.

## Bonus : Étape d'annotation
Si SnpEff est installé, le script ajoute une étape d'annotation fonctionnelle pour enrichir les fichiers VCF. Si cette étape échoue, vérifiez que le fichier de configuration SnpEff est correctement configuré.
