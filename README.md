# bwast
Performs command-line BLAST using Fasta/Genbank/EMBL files and loads the comparison up in ACT, *in just one step*.


Why use **bwast**?
-------------------

Allows:

* A subregion of genbank or fasta file to be specified after filename (start..stop)
* Automated loading of the sequences and comparison into the Artemis Comparison Tool (ACT) - act must be on your PATH (Please use ``-a`` flag)
* An infinite number of input sequences, depending on memory available

Produces (Output files):

* Fasta version of genbank/embl file supplied
* Blast tab delimited (-outfmt 6) output with blast options used, in the filename
* A new genbank/embl/fasta file of a specified subregion, if coordinates were given


Requirements
-----------------

**Dependencies**

* [NCBI's blast+] []
[NCBI's blast+]: http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

    (``blastn`` and ``tblastx`` need to be on your path)

* [Artemis Comparison Tool (ACT)] []
[Artemis Comparison Tool (ACT)]: https://www.sanger.ac.uk/resources/software/act/

    (``act`` needs to be on your path)

* [Biopython] []
[Biopython]: http://biopython.org/wiki/Main_Page

**Input format requirements**

* Input files should have a recognisable suffix (e.g. gb, gbk, emb, embl, fas, fasta, fna, fa)
* Input files should be present or sym-linked in the working directory. Please do not use relative/absolute path with filenames.
* Genbank and EMBL files must be formatted according to Biopython requirements. Valid headers must be present. Genbank output from Artemis and RAST are known to cause issues. 


Quick Start Instructions
--------------------------

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


F.A.Q
----------------

**I think I found a bug in the script. How do I let you know?**

Thanks for taking the time to report it! Please submit an issue on GitHub and I will be in touch shortly. You can also contact me on twitter [@bawee] [].

[@bawee]: https://twitter.com/bawee

**Why do I get the error: ``BLAST engine error: Empty CBlastQueryVector`` or ``Command line argument error: Query is Empty!``?**

Genbank files output by Artemis can sometimes cause this due to the absence of a valid header. Please contact me on twitter if you need a script to add a dummy header.


**Can I use relative/absolute paths to point to input files?**

The path name (e.g. ../../seq/) is currently read in when parsing the given arguments. It would be better to copy or link your sequence files into your current working directory.

E.g. 

```
cd directory/you/want/to/work/in

ln -s path/to/sequence/files .

path/to/bwast.py genbank1.gb genbank2.gb
```

**How do I get ACT on my PATH?**

On the Mac, run this command: 

```
export PATH="$PATH:/Applications/Artemis.app/Contents"
```

or put this line into your ~/.profile to have it set permanently.


**I get error messages when using a genbank or embl file. Is there something wrong with my file format?**

Biopython (used by bwast) is very fussy about the exact genbank/embl format used. Try manually adding the expected lines/text in the header. Genbank/embl files output by Artemis and RAST are known to be incompatible with biopython.


**My BLAST hits were flipped and are matching the wrong end of the sequence**

If your query and reference are very similar in length, ACT can sometimes flip the BLAST hits incorrectly. I do not know how to prevent this from happening. If you are specifying a subregion, try changing the sequence length.


**What does the name *bwast* mean?**

Blast Wrapper And Sequence Truncator 


Info
-----------
Written by Bryan Wee.

Inspired by EasyfigCL (MJ Sullivan) and Seqhandler (NF Alikhan)

Merge subroutine adapted from Seqhandler (NF Alikhan) github.com/happykhan/seqhandler

Version 0.0.1 - 1st version
Version 0.0.2 - added ability to detect and merge multifasta inputs - Thanks to B Forde for the suggestion.
