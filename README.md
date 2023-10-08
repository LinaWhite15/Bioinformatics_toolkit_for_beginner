# Bioinformatics toolkit for beginner
The utility is designed for processing protein sequences, as well as working with DNA sequences in the fastQ format (new functions for RNA and DNA processing will be added in future versions).

## Table of content 
+ [Overview](https://github.com/LinaWhite15/Bioinformatics_toolkit_for_beginner/edit/Development/README.md#overview)
+ [Installation](https://github.com/LinaWhite15/Bioinformatics_toolkit_for_beginner/edit/Development/README.md#installation)
+ [Usage](https://github.com/LinaWhite15/Bioinformatics_toolkit_for_beginner/edit/Development/README.md#usage)
+ [Credits](https://github.com/LinaWhite15/Bioinformatics_toolkit_for_beginner/edit/Development/README.md#credits)

## Overview
The toolkit contains a number of functions that allow you to filter data in the format of fastQ аnd also to analyze protein sequences according to a number of important physical and chemical properties. The programm is suitable for a wide range of users, including biologists with minimal knowledge of Python.

## Installation
To install the program, download files main_script.py and contents of the folder "modules".
OR

You can simple clone this repository using
```
git clone git@github.com:LinaWhite15/Bioinformatics_toolkit_for_beginner.git
```
(for Linux and WSL users)

**Python3 is required.**

## Usage
Before running the script, you must import additional modules
```
from modules.fastq_toolkit import *
from modules.protein_toolkit import *
```
### Input
For running **protein_toolkit** you must enter the protein sequence in one-letter or three-letter format and select one of the operations.
### List of operation in **protein_toolkit**
* ```content_check``` - Analyzes the amino acid composition of the protein. The output gives the percentage content of each molecule in the peptide. 

* ```seq_length``` - Measures the length of the peptide and gives the number of amino acids.

* ```protein_formula``` - Gives the atomic composition of the polymer.

* ```protein_mass``` - Сalculates the molecular mass of a protein in g/mol

* ```charge``` - Determines the charge of a protein when pH = 7.
### fastq_toolkit 
**fastq_toolkit** takes 4 arguments as input: ```seqs```, ```gc_bounds```, ```length_bounds```, ```quality_threshold```:
* ```seqs``` - a dictionary consisting of fastq sequences.  Key: string, containing name of sequence. The value is a tuple of two strings: sequence and quality.
* ```gc_bounds``` - composition GC interval (in percent) for filtering. The default is (0 :100). If you pass one number as an argument, other will be considered as the upper limit.
* ```length_bounds``` - length interval for filtering, by default it is equal to (0, 2**32).
* ```quality_threshold``` - threshold value of average read quality for filtering, default is 0 (phred33 scale).

### Output
* **protein_toolkit** - a string with result of performed operation.
* **fastq_toolkit** - a dictionary containing sequences corresponding to user-specified conditions.

## Credits
**Team:** 
Belikova Angelina - kiit@gmail.com Implemented: ```protein_formula```, ```protein_mass```, ```seq_length```, **fastq_toolkit**.

Aryuna Ayusheeva - aryuna.ayusheeva.1998@mail.ru Implemented: ```aa_content_check```, ```aa_chain_charge```.

Bredov Denis - d2707bredov@gmail.com Implemented: ```Mann_Whitney_U```, ```decomposition```, ```seq_transform```, ```check_and_procees_seq```, ```print_result```, ```run_protein_analyzer_tool```.
