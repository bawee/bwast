# bwast

**B**last **W**rapper **A**nd **S**equence **T**rimmer 

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

* [NCBI's blast+](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

    (``blastn`` and ``tblastx`` need to be on your path)

* [Artemis Comparison Tool (ACT)](http://www.sanger.ac.uk/science/tools/artemis-comparison-tool-act)
 
    (``act`` needs to be on your path)

* [Biopython](http://biopython.org/wiki/Main_Page)
    
    (I use pip to install biopython)

**Input format requirements**

* Input files should have a recognisable suffix (e.g. gb, gbk, gbff, emb, embl, fas, fasta, fna, fa)
* Genbank and EMBL files must be formatted according to Biopython requirements. Valid headers must be present. Genbank output from Artemis and RAST are known to cause issues. 


Quick Start Instructions
--------------------------

**Blast two genbank files**

``bwast.py <genbank1> <genbank2>``

**Blast more than two files**

``bwast.py <genbank1> <genbank2> <fasta1> <genbank3>``

**Blast two or more genbank files and load in Artemis Comparison Tool**

``bwast.py -a <genbank1> <genbank2>``

**Blast subregions of genbank files**

``bwast.py <genbank1> 200..2000 <genbank2> 7000..9000``

**Blast with blastn instead of megablast and use a custom e-value cutoff**

``bwast.py <genbank1> <genbank2> -f '-task blastn -evalue 0.0001'``

**See optional arguments and usage**

``bwast.py -h``


Detailed Instructions
----------------------------

**Getting bwast**

To use, simply use ``git clone https://github.com/bawee/bwast.git``, if you have Git installed (recommended).

*--or--*

Just download the ``bwast.py`` Python script on it's own.


**Running bwast**

``path/to/bwast.py sequence1.gbk sequence2.fna``


**Optional arguments**

* ``-h, --help`` show this help message and exit

* ``-f FLAGS, --flags FLAGS`` Custom BLAST options, enclosed in quotes. E.g. ``-f '-task blastn -evalue 0.001'``

* ``-v, --verbose`` Verbose mode

* ``-a, --act`` Run ACT after performing BLAST

* ``-b {blastn,tblastx}, --blast {blastn,tblastx}`` Blast program to use. Either ``tblastn`` or ``blastn``. Default is ``blastn``
* ``-o, --outfmt`` Change the output format. Default is ``--outfmt 6``. 


F.A.Q
----------------

1. **I think I found a bug in the script. How do I let you know?**

    Thanks for taking the time to report it! Please submit an [issue](https://github.com/bawee/bwast/issues) on GitHub and I will be in touch shortly. You can also contact me on twitter [@bawee](https://twitter.com/bawee).

2. **Why do I get the error: ``BLAST engine error: Empty CBlastQueryVector`` or ``Command line argument error: Query is Empty!``?**

    Genbank files output by Artemis can sometimes cause this due to the absence of a valid header. Please contact me on twitter if you need a script to add a dummy header.


3. **Can I use relative/absolute paths to point to input files?**

    Yes you can!

4. **How do I get ACT on my PATH?**

    On a Mac OSX, run this command: 

    ```
    export PATH="$PATH:/Applications/Artemis.app/Contents"
    ```

    or put this line into your ~/.profile to have it set permanently.

    On Windows 10, you can a combination of Windows Subsystem for Linux shell (e.g. ubuntu) and a window X-server such as VcXsrv or Xming. You will still need to add the ACT executable (for unix) to your path, for example:
    ```
    export PATH="$PATH:/path/to/artemis/download"
    ```

5. **I get error messages when using a genbank or embl file. Is there something wrong with my file format?**

    Biopython (used by bwast) is very fussy about the exact genbank/embl format used. Try manually adding the expected lines/text in the header. Genbank/embl files output by Artemis and RAST are known to be incompatible with biopython.


6. **My BLAST hits were flipped and are matching the wrong end of the sequence**

    If your query and reference are very similar in length, ACT can sometimes flip the BLAST hits incorrectly. I do not know how to prevent this from happening. If you are specifying a subregion, try changing the sequence length.


7. **What does the name *bwast* mean?**

    Blast Wrapper And Sequence Trimmer 


Info
-----------
Written by Bryan Wee.

Inspired by EasyfigCL (MJ Sullivan) and Seqhandler (NF Alikhan)

Merge subroutine adapted from Seqhandler (NF Alikhan) github.com/happykhan/seqhandler

1. Version 0.0.1 - 1st version
2. Version 0.0.2 - added ability to detect and merge multifasta inputs - Thanks to B Forde for the suggestion.
