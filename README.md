# bwast
Performs blast comparisons using genbank files

Description
--------------
Python script to run blast on Genbank/EMBL files without having to first convert to fasta. 

Allows:

* A subregion of genbank or fasta file to be specified after filename (start..stop)
* Automated loading of the sequences and comparison into the Artemis Comparison Tool (ACT) - act must be on your PATH (-a option)
* An infinite number of input sequences, depending on memory available

Requires:

* Files to have recognisable suffixes (i.e. gb, gbk, fa, fas, fna, fasta, emb, embl)
* Files to be present or sym-linked in current directory (i.e. NO relative/absolute paths to files)
* Properly formatted genbank/embl files (According to biopython requirements)

Produces:

* Fasta version of genbank/embl file supplied
* Blast tab delimited (-outfmt 6) output with blast options used, in the filename
* A new genbank/embl/fasta file of a specified subregion, if coordinates were given


Examples
---------------

**Blast two genbank files**

``bwast.py <genbank1> <genbank2>``

**Blast two or more files**

``bwast.py <genbank1> <genbank2> <fasta1> <genbank3>``

**Blast two or more genbank files and load in Artemis Comparison Tool**

``bwast.py -a <genbank1> <genbank2>``

**Blast subregions of genbank files**

``bwast.py <genbank1> 200..2000 <genbank2> 7000..9000``

**Blast with blastn instead of megablast and use a lower e-value cutoff**

``bwast.py <genbank1> <genbank2> -f '-task blastn -evalue 0.0001'``

**See optional arguments and usage**

``bwast.py -h``


Requirements
-----------------

* **NCBI's blast+** (``blastn`` and ``tblastx`` need to be on your path)
* **Artemis Comparison Tool (ACT)** (``act`` needs to be on your path)
* **Biopython**


F.A.Q
----------------

**Why do I get the error: ``BLAST engine error: Empty CBlastQueryVector`` or ``Command line argument error: Query is Empty!``?**

Genbank files output by Artemis can sometimes cause this due to the absence of a valid header. Please contact me if you need a script to add a dummy header.


**Can I use relative/absolute paths to point to input files?**

The path name (e.g. ../../seq/) is currently read in when parsing the given arguments. It would be better to copy or link your sequence files into your current working directory.

E.g. 

``cd directory/you/want/to/work/in``

``ln -s path/to/sequence/files .``

``path/to/bwast.py genbank1.gb genbank2.gb``


**How do I get ACT on my PATH?**

On the Mac, run this command: 

``export PATH="$PATH:/Applications/Artemis.app/Contents"``

or put this line into your ~/.profile to have it set permanently.


**I get error messages when using a genbank or embl file. Is there something wrong with my file format?**

Biopython (used by bwast) is very fussy about the exact genbank/embl format used. Try manually adding the expected lines/text in the header. Genbank/embl files output by Artemis and RAST are known to be incompatible with biopython.


**My BLAST hits were flipped and are matching the wrong end of the sequence**

If your query and reference are very similar in length, ACT can sometimes flip the BLAST hits incorrectly.


**What does the name *bwast* mean?**

Blast Wrapper And Sequence Truncator 


Info
-----------
Written by Bryan Wee.

Inspired by EasyfigCL (MJ Sullivan) and Seqhandler (NF Alikhan)

Merge subroutine adapted from Seqhandler (NF Alikhan) github.com/happykhan/seqhandler

Version 0.0.1 

