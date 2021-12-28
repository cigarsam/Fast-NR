# Fast-NR
|           |           |
| :-------- |  --------:| 
|  Name     |  Fast-NR  |
|  Author   |  Na He    |
|  Version  |  1.0      | 

Fast-NR is a software which applied on STARR-seq data to identify the negative regulatory elements, like silencer or insulator.

### Workflow
The workflow of Fast-NR as following:

![image](https://user-images.githubusercontent.com/81353915/112604452-1d4cce00-8e51-11eb-985c-9b5830a31fe6.png)


### Installation
You can install Fast-NR either with using pip or git repository.

  
1. Using git repository

  `git clone https://github.com/Na-He/Fast-NR.git`
  
  `unzip Fast-NR-main.zip`
  
  `cd Fast-NR-main`
  
  `python setup.py install`

2. Using git repository

   `pip install git+https://github.com/Na-He/Fast-NR.git`
  
### Dependencies

  Fast-NR requires
  
*  numpy 
*  pysam
*  scipy 


### Command

`Fast-NR -t <controlFile> -c <treatmentFile> -g <genomeFile> -o <outputName> [options]`

### Parameters

#### requires:

| parameter | type      | Description|
| :-------- |  :--------| :----------|
| `-t`      |  *file*     | Control file path and name, usually the STARR-seq plasmid library. Only .bam and .bed file accpected. |
| `-c`      |  *file*     | Treatment file path and name, usually the STARR-seq cDNA library. Only .bam and 3 column or 6 column .bed file accpected. But this program will not consider the strand information. |
| `-g`      |  *file*     | Genome chromosome size file, which can download from the UCSC. File include two column, chromosome, length of chromosome. |
| `-o`      |  *file*     | Output path and prefix name of the final result peaks. |

#### options:

| parameter | type      | Description|
| :-------- |  :--------| :----------|
| `-p`      |  *value*    | The cut-off of p value. Default 10<sup>-5</sup>.|
| `-cp`     |  *string*    | The correct method of p value. Include "BH", "Bonferroni". Default "Bonferroni". |
| `-cm`     | *string*    | The method used to calculate the curve similarity. Include 'Cosine', 'Pearson', 'Euclidean', and 'Gradiente'. Default 'Cosine'.|
| `-ct`     | *value*    | The percent cut-off of similarity distance. From 0 to 1. Default 0.9.|
| `-ws`     |  *value*    | The size of window, used to find differential coverage region. Default 600 bp.|
| `-l`      | *value*     | Extend the fragment length to this fixed bp length. Default 0.|
| `--remove`     | *logical*   | Remove the slop or not. Set to "TRUE" if want to remove the results in slop regions. Default "FALSE".|

### Output Format

The final output file is a .bed file, include 9 columns.

*  column 1 | *chromosome*. The chromosome name.
*  column 2 | *start*. The start position of peak region.
*  column 3 | *end*. The end position of peak region.
*  column 4 | *summit*. Summit of peak region, also the mid position of peak.
*  column 5 | *tCount*. Number of fragments in treatment library in peak region.
*  column 6 | *cCount*. Number of fragments in control library in peak region.
*  column 7 | *FC*. The fold change of fragment counts in peak region compared with treatment count and control count. 
*  column 8 | *P*. P value of peak region.
*  column 9 | *cP*. Corrected p value.

***Example:***

`chr2L	103500	104100	103800	120	224	0.535714	2.198757e-13	3.493824e-10`

`chr2L	114264	114864	114564	8	112	0.071429	1.877090e-36	2.982696e-33`

`chr2L	143773	144373	144073	260	440	0.590909	7.688970e-19	1.221777e-15`

`chr2L	148969	149569	149269	139	329	0.422492	1.780870e-30	2.829802e-27`


***
### Citation:
Once you used this program, please cite the following paper:

(https://github.com/Na-He/Fast-NR)

