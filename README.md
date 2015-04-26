# bwast
Performs blast comparisons using genbank files

Description
--------------
Python script to run blast on Genbank/EMBL files without having to first convert to fasta. 

* Allows subregion of genbank or fasta file to be specified after filename (start..stop)
* Allows loading of files into Artemis Comparison Tool (ACT) - act must be on your PATH (-a option)
* Can take in infinite number of input sequences, depending on memory available

* Requires files to have recognisable suffixes (i.e. gb, gbk, fa, fas, fna, fasta, emb, embl)
* Requires files to be present or sym-linked in current directory

* Produces fasta version of genbank/embl file
* Produces blast tab delimited (-outfmt 6) output with blast options specified
* Produces genbank/embl/fasta of subregion, if coordinates given



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

**Can I use relative/absolute paths to point to input files?**

The path name (e.g. ../../seq/) is currently read in when parsing the given arguments. It would be better to copy or link your sequence files into your current working directory.

E.g. 

``cd directory/you/want/to/work/in``

``ln -s path/to/sequence/files .``

``path/to/bwast.py genbank1.gb genbank2.gb``

**How do I get ACT on my PATH?**

On the Mac, add this line to your .profile file (~/.profile): 

``export PATH="$PATH:/Applications/Artemis.app/Contents``

Then reload your .profile:

``source ~/.profile``

**My BLAST hits were flipped and are matching the wrong end of the sequence**

If your query and reference are very similar in length, ACT can sometimes flip the BLAST hits incorrectly.

**What does the name *bwast* mean?**

Blast Wrapper And Sequence Truncator 


Info
-----------
Written by Bryan Wee.

Inspired by EasyfigCL (MJ Sullivan) and Seqhandler (NF Alikhan)


Version 0.0.1 


