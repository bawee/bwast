# bwast
Performs blast comparisons using genbank files

Description
--------------
Python script to run blast on Genbank/EMBL files without having to first convert to fasta. 

Usage
-----------

**Blast 2 genbank files**

``bwast.py <genbank1> <genbank2>``

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

**How do I get ACT on my PATH?**

On the Mac, add this line to your .profile: 
``export PATH="$PATH:/Applications/Artemis.app/Contents``



Written by Bryan Wee.


Version 0.0.1 

"""
