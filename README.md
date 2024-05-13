# Sandrine_internship
allele frequency, here you can find all code and data

## Text files

### All_bam_list

Here you will file the bam file name of every population used in this study. It is used to create a sync file and to know the position of each population while separating the sync file.

### Create_sync_file

In this text file, there is a bash code used to create data used in the simulation code: Create sync file taking all populations together, cleaning data for every chromosome number and concatenation, filtering where one generation in one replicate is taking with his respective ancestor together.


## Data

### Reads

Here is an Excel file representing every population with its corresponding number of reads before trimming, after trimming and also after mapping respectively.

### Fst-value

Here is a two-collum table, the first collum gives the population name from where you can access his type, replicate and generation. The second collum gives the obtained Fst value between the given population and his respective ancestor.


## Codes

### Save_Polymorphic_sites

Here is the R code used for cleaning the data process on bash in Create_sync_file.

### Bash.ssh

In this text file, there is a bash code used to create data used in the simulation code: Create sync file taking all populations together, cleaning data for every chromosome number and concatenation, filtering where one generation in one replicate is taking with his respective ancestor together.


### fst

Here is the R code simulating Fst between each population with its ancestor. here it gives the average and multilocus FST. This code allows plotting Fst value over generation and along the genome.

### Ne

This code corresponds to Ne's estimation for each population considering his effective size and his ancestor's effective sizes.

### CA_cmh_chisq

This code performs the Cochran-Mantel-Haenszel test and chi-square test of Dioecy


### CA_cmh_chisq
This code performs the Cochran-Mantel-Haenszel test and chi-square test of Androdioecy

