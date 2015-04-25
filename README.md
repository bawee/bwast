# bwast
Performs blast comparisons using genbank files

Description
--------------
Python script to run blast on Genbank/EMBL files without having to first convert to fasta. 

Usage
-----------

**Blast 2 genbank files**

``bwast.py <genbank1> <genbank2>``

**Blast 2 or more files**

``bwast.py <genbank1> <genbank2> <fasta1> <genbank3>``

**Blast 2 or more genbank files and load in Artemis Comparison Tool**

``bwast.py -a <genbank1> <genbank2>``

**Blast subregions of genbank files**

``bwast.py <genbank1> 200..2000 <genbank2> 7000..9000``


Cool features: 
------------------
* Allows subregion of genbank or fasta file to be specified
* Allows loading of files into Artemis Comparison Tool (ACT) - act must be on your PATH

Requirements
-----------------

blast+ (blastn and tblastx need to be on your path)
act (act needs to be on your path)
biopython

Cut and converted files will be written into the same directory as the file specified in the arguments, NOT the current working directory.




F.A.Q
----------------

**Can I use relative/absolute paths to point to input files?**
The path name (e.g. ../../seq/) is currently read in when parsing the given arguments. It would be better to copy or link your sequence files into your current working directory.

**How do I get ACT on my PATH?**

On the Mac, add this line to your .profile file (~/.profile): 
``export PATH="$PATH:/Applications/Artemis.app/Contents``

Then reload your .profile:
``source ~/.profile``



Written by Bryan Wee.


Version 0.0.1 

"""
