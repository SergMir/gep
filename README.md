Research of Gene Expression Programming (GEP) algorithm
---
This project is yet another implementation of basic GEP, invented and described by Candida Ferreira, and modifications by different authors (including [me](https://github.com/SergMir)) related to:
* Coding and representations of syntactic trees;
* Selection methods;
* Numerical constants;
* Fitness fuctions;
* Applying new features: incrementation of chromosome's length, differential models etc.

GEP, as an evolutionary computation techique, is a highly randomized algorithm, so needs to be executed multiple times to generate several models and best one is picked.

This framework allows to setup configuration of different features and gather statistical characteristics of chosen config: probability of finding good/appropriate solution, numerical value of quality of solutions (usually mean squared error), etc.

File structure
---
* papers - published, drafts and draft of PhD thesis (to be finished in 2015), LaTex, mainly in russian;
* stats - gathered statistics of different algorithm's configurations;
* src
  * test_files_gen - generates samples files (set of coordinate points) of known functions, like sin, Gabor, Rosenbrock, that are next to be passed to GEP;
  * gep_src - GEP library (written in C), API described at include/gep.h;
  * test.c - example of usage and testing procedure;
  * samples_comparator.c - tool to compare samples files;
  * scripts - generates statistics.

Links
---
[Website of GEP's author](http://www.gene-expression-programming.com/)

[List of used papers](papers/draft_thesis_half/related_work.bib)
