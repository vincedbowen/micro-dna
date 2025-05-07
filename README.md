# Micro DNA

Add data files from the public google drive [directory](https://drive.google.com/drive/folders/1WYPwCzSv__28iQlHwvZUZpNlmjTM7Fea) into the data folder.
```
  $ cd pairwise-aligner-comparer
  $ python compare_pwa.py
```
![image](https://github.com/user-attachments/assets/e149035d-4c0c-41e1-b397-25a744adc24b)

```
$ cd micro-dna-finder/cli
$ python3 main.py \
-bf ../data/SRR413984.sorted.NC_000001.10.bam \
-ff ../data/GCF_000001405.13_GRCh37_genomic.NC_000001.10.fna \
-of ../results/results.txt
```
![image](https://github.com/user-attachments/assets/a57e6111-bf06-4e5a-956d-038dee66786a)

```
$ cd micro-dna-finder/cli
$ python3 main.py \
-bf ../data/SRR413984.sorted.NC_000001.10.bam \
-ff ../data/GCF_000001405.13_GRCh37_genomic.NC_000001.10.fna \
-of ../results/results.csv
```
![image](https://github.com/user-attachments/assets/45d7517d-edce-4adb-854e-766a29a30f06)
