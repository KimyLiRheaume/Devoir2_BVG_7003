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

## Fichiers d'entrée
- **Fichier FASTQ** : Données brutes de séquençage (`Data/FC20150701_1.fq.gz`).
- **Fichier de barcodes** : Pour le démultiplexage (`Data/FC20150701_1.txt`).
- **Génome de référence** : Fichier FASTA du génome de référence (`Data/Gmax_275_v2.0.fa`).

Si vous souhaitez utiliser le pipeline avec d'autres fichiers de données, le pipeline a été conçu pour reconnaitre des fichiers avec des extensions .fq.gz, .fa.gz et .txt. 

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



## Installation
1. Clonez ce dépôt GitHub :```git clone https://github.com/KimyLiRheaume/Devoir2_BVG_7003.git```
2. Assurez-vous que toutes les dépendances sont installées et configurées correctement.
3. Placez vos fichiers d'entrée dans le dossier data/.

## Utilisation
Pour exécuter le pipeline, utilisez la commande suivante :
```
bash ScriptFinal.sh
```




## Structure du dépôt 

variant-calling-pipeline/
```
│
├── Data/                   # Fichiers d'entrée (FASTQ, barcodes, génome de référence)
├── Résultats/              # Résultats générés par le pipeline
├── Script/                # Script principal du pipeline
├── Logs/                   # Fichiers journaux pour chaque étape
└──README.md                # Documentation du pipeline

```

## Dépannage 
1. Vérifiez les fichiers journaux dans le dossier logs/ pour identifier les erreurs.
2. Assurez-vous que vos fichiers d'entrée sont correctement formatés.
3. Testez chaque étape individuellement si nécessaire.
4. Consultez la section Wiki pour plus d'informations et pour collaborer. 


