# blasperin

Tool to BLAST Fasta files.


## Getting Started

**blasperin** is a crucial part of a major set of packages: 

* keggreader
* keggimporter
* anendbtofasta
* clusperin
* clusteringloader

## Requirements

To run **blasperin** you only need the **BLAST** software.


## Installing

* Download (and extract) the zip file from this repository or clone it using **git** command.
* Go to the opened directory.
* Run:

```
python setup.py install
```

* Installation is done.


## Run blasperin

* Simply type **blasperin** in your console.

**Example:**

```
blasperin --label testing \
          --author 'your name' \
          --fasta-files-directory /var/clustering \
          --log-file /var/log/blaspering.log
```

## Explaining the parameters

* --label

A name to identify the process.

* --author

A name to identify who is executing the process.

* --fasta-files-directory

Full path to where's the Fasta files to be blasted.
 
* --log-file

Full path the a log file. That's necessary because the process can take several hours/days depending on the amount of files.



