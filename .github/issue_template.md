<!--

If you are reporting a bug, please provide:

* version of sambamba

* version of samtools if the issue is with pileup tool

* the command line or bash script used (feel free to shorten paths and filenames)

* whenever possible, a set of BAM/BED files to reproduce the issue
  * bonus points if you try to minimize the test case yourself, as issues are often localized:
    - try to use sambamba slice to first extract the reference where the error occurs
    - if that succeeds (the error is still reproducible), continue to crop the file in binary-search fashion
  * suggestions regarding privacy concerns for human data:
    - try to minimize the test case (few variants => much harder to tell anything about the patient)
    - write a script that changes the reference and/or applies a fixed shift to all record & record mate positions
    - apply a random base permutation, e.g. A->C, C->T, T->A, G->G
    
You can also suggest a new feature, but please understand that neither of maintainers uses sambamba in daily jobs anymore.
As such, while it has some chances to get picked up by another person, if you possess any programming skills whatsoever, we encourage you to try your hand at it! (Feel free to ask us questions about the codebase or D quirks over email)
      
-->
