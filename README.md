# MACHINA - Metastatic And Clonal History INtegrative Analysis

MACHINA is a computational framework for inferring migration patterns between a primary tumor and metastases using DNA sequencing data.
![Overview of MACHINA](doc/overview.png)

## Contents

  1. [Installation](#installation)
     * [Bioconda](#bioconda)
     * [Manual compilation](#compilation)
  2. [Usage instructions](#usage)
     * [I/O formats](#io)
       - [Clone tree](#clonetree)
       - [Leaf labeling](#leaflabeling)
       - [Vertex labeling](#vertexlabeling)
       - [Frequencies](#frequencies)
     * [Parsimonious Migration History](#pmh)
     * [Parsimonious Migration History with Tree Resolution](#pmh_tr)
     * [Parsimonious Migration History with Tree Inference](#pmh_ti)

<a name="installation"></a>
## Installation 

<a name="bioconda"></a>

### bioconda
1. Install [Anaconda](https://docs.anaconda.com/anaconda/install/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) if you do not already have one installed.
2. (recommended) Create a new conda environment for `machina` and activate it:
```
conda create -n machina
conda activate machina
```
3. Set up `conda` channels for `bioconda` (once per Anaconda/Miniconda installation):
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
4. Install `machina` from bioconda:
```
conda install machina
```

<a name="compilation"></a>

### Manual compilation
*Note that binaries for macOS and linux are available [here](https://github.com/raphael-group/machina/releases). These binaries require a valid Gurobi installation and license key. License key location can be specified via the environment variable GRB_LICENSE_KEY. In addition, installation of Gurobi in a non-standard location will require updating LD_LIBRARY_PATH (linux) and DYLD_LIBRARY_PATH (macOS).*

Also note that to run the below examples, you must either provide the full path to the executable (e.g., `/path/to/machina/build/pmh_sankoff`) or add the `build` directory to your PATH.

<a name="dep"></a>

#### Dependencies

MACHINA is written in C++11 and thus requires a modern C++ compiler (GCC >= 4.8.1, or Clang). In addition, MACHINA has the following dependencies.

* [CMake](http://www.cmake.org/) (>= 3.0)
* [Boost](http://www.boost.org) (>= 1.38)
* [LEMON](http://lemon.cs.elte.hu/trac/lemon) graph library (>= 1.3)
* [Gurobi](http://www.gurobi.com) (>= 6.0)

[Graphviz](http://www.graphviz.org) is required to visualize the resulting DOT files, but is not required for compilation.

[Gurobi](http://www.gurobi.com) is a commercial ILP solver with two licensing options: (1) a single-host license where the license is tied to a single computer and (2) a network license for use in a compute cluster. Both options are freely available for users in academia.

In case [doxygen](http://www.stack.nl/~dimitri/doxygen/) is available, extended source code documentation will be generated.

<a name="comp"></a>

#### Compilation

To compile MACHINA, execute the following commands from the root of the repository:

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make

In case CMake fails to detect LEMON or Gurobi, run the following command with adjusted paths:

    $ cmake -DLIBLEMON_ROOT=~/lemon \
    -DGUROBI_HOME=/path/to/gurobiXXX

where `XXX` is the 3-digit version of gurobi.

The compilation results in the following files in the `build` directory:

COMMAND | DESCRIPTION
-----------|-------------
`cluster` | Cluster mutations using a combinatorial algorithm that models variant read counts using a binomial distribution.
`generatemigrationtrees` | Generates all migration trees given anatomical site labels. These migration trees can be used to constrain the search space of the `pmh`, `pmh_pr` and `pmh_cti` algorithms.
`generatemutationtrees` | Generates all mutation trees given a frequency matrix.
`pmh_sankoff`  | Enumerates all minimum-migration vertex labelings given a clone tree.
`pmh` | Solves the Parsimonious Migration History (PMH) problem given a migration pattern restriction and a clone tree.
`pmh_tr`    | Solves the Parsimonious Migration History with Polytomy Resolution (PMH-PR) problem given a migration pattern restriction and a clone tree.
`pmh_ti`     | Solves the Parsimonious Migration History and Tree Inference (PMH-TI) given a migration pattern restriction, a mutation tree and a frequency matrix.
`simulate` | Simulates a metastatic tumor.
`visualizeclonetree` | Visualizes a clone tree and optional vertex labeling.
`visualizemigrationgraph` | Visualizes the migration graph given a clone tree and vertex labeling.

<a name="usage"></a>
## Usage instructions

<a name="io"></a>
### I/O formats

Below we describe the various formats used by the algorithms of the `MACHINA` framework.

<a name="clonetree"></a>
#### Clone tree

A clone tree is provided as an edge list. Each line specifies an edge by listing the labels of the incident vertices separated by a space or tab character. For example:

    A A1
    A A2
    A A3
    A A4
    A A5
    A A6
    ...

See [patient1.tree](data/mcpherson_2016/patient1.tree) for the complete clone tree.

<a name="leaflabeling"></a>
#### Leaf labeling

A leaf labeling assigns an anatomical site label to each leaf of a clone tree. Each line contains two values, the leaf label and the anatomical site label separated by a space or tab character. For example:

    A1 Om
    A2 SBwl
    A3 LFTB
    A4 LOv
    A5 ApC
    A6 RFTA
    ...

See [patient1.labeling](data/mcpherson_2016/patient1.labeling) for the complete leaf labeling.

<a name="vertexlabeling"></a>
#### Vertex labeling

A vertex labeling assigns an anatomical site label to each vertex of a clone tree (including the leaves). Each line contains two values, the vertex label and the anatomical site label separated by a space or tab character. For example:

    A ROv
    B SBwl
    D ROv
    F ROv
    H ROv
    A1 Om
    A2 SBwl
    ...

See [patient1.reported.labeling](data/mcpherson_2016/patient1.reported.labeling) for the complete vertex labeling.

#### Frequencies

A frequency file encodes the frequency of every mutation (cluster) in an anatomical site (sample). It is a tab separated file. The first line lists the number of anatomical sites followed by the number of samples and then the number of mutations, each on separate lines. The fourth line is ignored but describes the format of the rest of the file. Each subsequent line encodes the cell frequency of a mutation in a sample: first the sample 0-based index is given, followed by the label of the sample, the 0-based index of the anatomical site, the anatomical site label, the 0-based index of the mutation, the label of the mutation, the frequency lower bound and upper bound.

    6 #anatomical sites							
    6 #samples							
    10 #mutation clusters							
    #sample_index	sample_label	anatomical_site_index	anatomical_site_label	character_index	character_label	f_lb	f_ub
    0	breast	0	breast	0	1	0.503628522	0.545237495
    0	breast	0	breast	1	2	0	0.01213794
    ...

See [F.tsv](data/hoadley_2016/A7/A7_MACHINA_0.95.tsv) for the complete frequency file. For an example on how obtain this file from read data, please see: [process_A7_new.ipynb](data/hoadley_2016/A7/raw/process_A7_new.ipynb). Specifically, you will need to process the bulk DNA sequencing data by first calling single-nucleotide variants and copy-number aberrations. Then, SNVs that occur in copy-neutral regions need to be clustered (e.g., using SciClone or PyClone). Confidence intervals can then be obtained by first pooling for each sample the read counts of the mutations that belong to the same cluster followed by using a beta distribution. Please see the supplement of the MACHINA paper for more details.

<a name="pmh"></a>
### Parsimonious Migration History (`pmh_sankoff` and `pmh`)

In the parsimonious migration history we are given a clone tree `T` whose leaves are labeled by anatomical sites. The task is to label the inner vertices of `T` such that the resulting migration graph `G` has minimum migration number and comigration number. Additionally, it is possible to specify constraints on the topology of the migration graph.

PATTERN | DESCRIPTION
--------|------------
PS (parallel single-source seeding) | Each metastatic site is seeded directly from the primary tumor, i.e. `G` is a multi-tree such that the primary `P` is the only vertex with out-degree greater than 1.
S (single-source seeding) | Each metastatic site is seeded from only one other anatomical site, i.e. `G` is a multi-tree.
M (multi-source seeding) | A metastatic site may be seeded from multiple anatomical sites, but no directed cycles are introduced. That is, `G` is multi-DAG.
R (reseeding) | Directed cycles in `G` are allowed.

In our algorithms we allow for the following restrictions on the migration pattern:

1. Unrestricted: PS, S, M and R
2. No reseeding: PS, S, M
3. No reseeding and no multi-source seeding: PS and S
4. No reseeding, no multi-source seeding and no single-source seeding: PS

The unconstrained PMH problem can be solved by running `pmh_sankoff`, which is an adaptation of the Sankoff algorithm and enumerates all migration histories:

    Usage:
    pmh_sankoff [--help|-h|-help] [-c str] [-o str] [-p str] T leaf_labeling
    Where:
    T
        Clone tree
    leaf_labeling
        Leaf labeling
    --help|-h|-help
        Print a short help message
    -c str
        Color map file
    -o str
        Output prefix
    -p str
        Primary anatomical sites separated by commas (if omitted, every
        anatomical site will be considered iteratively as the primary)

An example execution of the `pmh_sankoff` algorithm (executed from the root directory of the MACHINA repository):

    $ mkdir patient1
    $ pmh_sankoff -p LOv,ROv -c data/mcpherson_2016/coloring.txt data/mcpherson_2016/patient1.tree \
    data/mcpherson_2016/patient1.labeling -o patient1/ 2> patient1/result.txt
    
    $ cat patient1/result.txt
    Clone tree has 7 anatomical sites
    Found 4 maximum parsimony labelings with primary 'LOv'
    Found 2 labelings with 7 comigrations, 2 seeding sites and pR
    Found 2 labelings with 11 comigrations, 3 seeding sites and pR
    Labeling 0: 13 migrations, 11 comigrations, 3 seeding sites and pR
    Labeling 1: 13 migrations, 11 comigrations, 3 seeding sites and pR
    Labeling 2: 13 migrations, 7 comigrations, 2 seeding sites and pR
    Labeling 3: 13 migrations, 7 comigrations, 2 seeding sites and pR
    Found 1 maximum parsimony labelings with primary 'ROv'
    Found 1 labelings with 10 comigrations, 2 seeding sites and pM
    Labeling 0: 13 migrations, 10 comigrations, 2 seeding sites and pM

The above command considers the left ovary (LOv) and right ovary (ROv) as the primary tumor site and enumerates all minimum migration vertex labelings of the given clone tree and leaf labeling. The output is stored in the directory patient1. In the DOT files the given color map is used for coloring the anatomical sites.

To further constrain the migration graph, we can use `pmh`:


    Usage:
      pmh [--help|-h|-help] [-G str] [-OLD] [-UB_gamma int] [-UB_mu int]
         [-UB_sigma int] -c str [-e] [-g] [-l int] [-log] [-m str] [-o str]
         -p str [-t int] T leaf_labeling
    Where:
      T
         Clone tree
      leaf_labeling
         Leaf labeling
      --help|-h|-help
         Print a short help message
      -G str
         Optional file with migration graphs
      -OLD
         Use old ILP (typically much slower)
      -UB_gamma int
         Upper bound on the comigration number (default: -1, disabled)
      -UB_mu int
         Upper bound on the migration number (default: -1, disabled)
      -UB_sigma int
         Upper bound on the seeding site number (default: -1, disabled)
      -c str
         Color map file
      -e
         Export ILP
      -g
         Output search graph
      -l int
         Time limit in seconds (default: -1, no time limit)
      -log
         Gurobi logging
      -m str
         Allowed migration patterns:
           0 : PS
           1 : PS, S
           2 : PS, S, M
           3 : PS, S, M, R
         If no pattern is specified, all allowed patterns will be enumerated.
      -o str
         Output prefix
      -p str
         Primary anatomical site
      -t int
         Number of threads (default: -1, #cores)


An example execution of the `pmh` algorithm (executed from the root directory of the MACHINA repository):

    $ mkdir patient1_constrained
    $ pmh -p LOv -c data/mcpherson_2016/coloring.txt data/mcpherson_2016/patient1.tree \
    data/mcpherson_2016/patient1.labeling -o patient1_constrained/ > patient1_constrained/result.txt
    
    $ cat patient1_constrained/result.txt
    LOv-    (PS)    15      6       1       pPS     15.125  15.125  0.106424
    LOv-    (PS, S) 15      6       1       pPS     15.125  15.125  0.057286
    LOv-    (PS, S, M)      15      6       1       pPS     15.125  15.125  0.060853
    LOv-    (PS, S, M, R)   13      7       2       pR      13.148  13.148  0.44595

Each line lists the solution found by MACHINA. First the primary anatomical site is given, then the provided migration pattern restriction set, followed by the migration number, comigration number and seeding site number. Finally, the identified migration pattern is given, followed by a lower bound (LB) on the optimal solution and then an upper bound (UB), ending with the total running time in seconds. In case LB == UB, the identified solution is optimal.

<a name="pmh_tr"></a>
### Parsimonious Migration History with Tree Resolution (`pmh_tr`)

In the parsimonious migration history with polytomy resolution we are given a clone tree `T` whose leaves are labeled by anatomical sites. The task is to find a refinement `T'` of `T` and label its inner vertices such that the resulting migration graph `G` has minimum migration number, comigration number and seeding site number. It is possible to specify constraints on the topology of the migration graph. 

**(Update 10/4/21)** - Now it can enumerate through multiple optimal/sub-optimal solutions. the flag `-N` can be used to specify the number of solutions to fetch. One problem is that `pmh_tr` outputs multiple identical trees when their internal representation inside `pmh_tr` is different. To remove such duplicate trees, `-P` flag can be used.

**(Update 11/9/21)** - Symmetry breaking constraints added, so no need to do post-processing. 

**(Update 11/15/21)** - A new flag, '-C,' has been added to count the number of solutions without printing the trees. The flag outputs the precise size of the solution space if '-N' is not used (assuming the number of solution is less than 2,000,000,000). Using **pmh tr** without '-N' may be slow. The count will be bounded by the value of '-N' if '-N' is used. Using `-C` with `-N` may be useful to know if the number of solution is greater than a specific threshold.

    Usage:
      pmh_tr [--help|-h|-help] [-G str] [-OLD] [-UB_gamma int] [-UB_mu int]
         [-UB_sigma int] -c str [-e] [-g] [-l int] [-log] [-m str] [-o str]
         -p str [-t int] T leaf_labeling [-N int] [-P]
    Where:
      T
         Clone tree
      leaf_labeling
         Leaf labeling
      --help|-h|-help
         Print a short help message
      -C
         Only output the number of solutions (default: False)
      -G str
         Optional file with migration graphs
      -OLD
         Use old ILP (typically much slower)
      -UB_gamma int
         Upper bound on the comigration number (default: -1, disabled)
      -UB_mu int
         Upper bound on the migration number (default: -1, disabled)
      -UB_sigma int
         Upper bound on the seeding site number (default: -1, disabled)
      -c str
         Color map file
      -e
         Export ILP
      -g
         Output search graph
      -l int
         Time limit in seconds (default: -1, no time limit)
      -log
         Gurobi logging
      -m str
         Allowed migration patterns:
           0 : PS
           1 : PS, S
           2 : PS, S, M
           3 : PS, S, M, R
         If no pattern is specified, all allowed patterns will be
         enumerated (default: '0,1,2,3')
      -o str
         Output prefix
      -p str
         Primary anatomical site
      -t int
         Number of threads (default: -1, #cores)
      -N int
         Number of solutions (default: 1)
      -P 
         Enable post-processing (default: False)

An example execution (executed from the root directory of the MACHINA repository):


    $ mkdir patient1_tr
    $ pmh_tr -p LOv -c data/mcpherson_2016/coloring.txt data/mcpherson_2016/patient1.tree \
    data/mcpherson_2016/patient1.labeling -o patient1_tr/ -N 5 > patient1_tr/result.txt
    
    $ cat patient1_tr/result.txt
    LOv-    (PS)    0       12      6       1       pPS     12125   12125   0.302463
    LOv-    (PS)    1       12      6       1       pPS     12125   12125   0.303952
    LOv-    (PS)    2       12      6       1       pPS     12125   12125   0.305179
    LOv-    (PS)    3       12      6       1       pPS     12125   12125   0.3064
    LOv-    (PS)    4       12      6       1       pPS     12125   12125   0.307605
    LOv-    (PS, S) 0       12      6       1       pPS     12125   12125   0.50031
    LOv-    (PS, S) 1       12      6       1       pPS     12125   12125   0.501627
    LOv-    (PS, S) 2       12      6       1       pPS     12125   12125   0.502769
    LOv-    (PS, S) 3       12      6       1       pPS     12125   12125   0.503921
    LOv-    (PS, S) 4       12      6       1       pPS     12125   12125   0.505132
    LOv-    (PS, S, M)      0       12      6       1       pPS     12125   12125   0.509809
    LOv-    (PS, S, M)      1       12      6       1       pPS     12125   12125   0.511246
    LOv-    (PS, S, M)      2       12      6       1       pPS     12125   12125   0.512531
    LOv-    (PS, S, M)      3       12      6       1       pPS     12125   12125   0.513785
    LOv-    (PS, S, M)      4       12      6       1       pPS     12125   12125   0.515015
    LOv-    (PS, S, M, R)   0       11      7       2       pR      11148   11148   3.53926
    LOv-    (PS, S, M, R)   1       11      7       2       pR      11148   11148   3.54072
    LOv-    (PS, S, M, R)   2       11      7       2       pR      11148   11148   3.54199
    LOv-    (PS, S, M, R)   3       11      7       2       pR      11148   11148   3.54325
    LOv-    (PS, S, M, R)   4       11      7       2       pR      11148   11148   3.54452

Each line lists the solution found by `pmh_tr`. First the primary anatomical site is given, then the provided migration pattern restriction set, followed by the solution number,  migration number, comigration number and seeding site number. Finally, the identified migration pattern is given, followed by a lower bound (LB) on the optimal solution and then an upper bound (UB), ending with the total running time in seconds. In case LB == UB, the identified solution is optimal. Here, the input number of solutions is 5 and there's no post-processing, so the top 5 solutions are being reported for each pattern. In case `-P` flag is activated, `pmh_tr` will remove the duplicate trees, so the number of actual solutions will often be less than the input number of solutions. 

```
$ mkdir patient1_tr
$ pmh_tr -p LOv -c data/mcpherson_2016/coloring.txt data/mcpherson_2016/patient1.tree \
data/mcpherson_2016/patient1.labeling -o patient1_tr/ -N 5 -P > patient1_tr/result.txt

$ cat patient1_tr/result.txt
LOv-    (PS)    0       12      6       1       pPS     12125   12125   0.439512
LOv-    (PS, S) 0       12      6       1       pPS     12125   12125   0.438603
LOv-    (PS, S, M)      0       12      6       1       pPS     12125   12125   0.53151
LOv-    (PS, S, M, R)   0       11      7       2       pR      11148   11148   4.28141
```

Here only one solution has been reported for each pattern as other four solutions were outputting the same tree. 

If '-C' is used, **pmh tr** prints the primary anatomical site, the specified migration pattern, the total number of feasible solutions or the size of solution space, and the number of optimal solutions in that order. If '-N' is used with `-C`, the count is bounded by the value of `-N`.

```
$ mkdir patient1_tr
$ pmh_tr -p LOv -c ../data/mcpherson_2016/coloring.txt ../data/mcpherson_2016/patient1.tree\ 
../data/mcpherson_2016/patient1.labeling -o patient1_tr -C > patient1_tr/result.txt

$ cat patient1_tr/result.txt
LOv-    (PS)    8       1
LOv-    (PS, S) 8       1
LOv-    (PS, S, M)      36      1
LOv-    (PS, S, M, R)   33248   1
```

Here, the number of solutions for **(PS, S, M, R)** pattern has been bounded by the value of $N=1000$.
```
$ pmh_tr -p LOv -c ../data/mcpherson_2016/coloring.txt ../data/mcpherson_2016/patient1.tree\ 
../data/mcpherson_2016/patient1.labeling -o patient1_tr -C -N 1000 > patient1_tr/result.txt

$ cat patient1_tr/result.txt
LOv-    (PS)    8       1
LOv-    (PS, S) 8       1
LOv-    (PS, S, M)      36      1
LOv-    (PS, S, M, R)   1000    1
```

<a name="pmh_ti"></a>
### Parsimonious Migration History with Tree Inference (`pmh_ti`)

Given a mutation tree `T` with mutation frequencies `F-` and `F+`, the task is to find a frequency assignment `F` yielding a refined clone tree `T'` that admits a vertex labeling `l` such that the resulting migration graph `G` has minimum number of migrations and co-migrations. It is possible to specify constraints on the topology of the migration graph.


    Usage:
      pmh_ti [--help|-h|-help] -F str [-G str] [-OLD] [-UB_gamma int]
         [-UB_mu int] [-UB_sigma int] -barT str -c str [-e] [-g] [-l int] [-log]
         [-m str] [-mutTreeIdx int] [-noPR] [-o str] -p str [-t int]
    Where:
      --help|-h|-help
         Print a short help message
      -F str
         Frequencies file
      -G str
         Optional file with migration graphs
      -OLD
         Use old ILP (typically much slower)
      -UB_gamma int
         Upper bound on the comigration number (default: -1, disabled)
      -UB_mu int
         Upper bound on the migration number (default: -1, disabled)
      -UB_sigma int
         Upper bound on the seeding site number (default: -1, disabled)
      -barT str
         Mutation trees
      -c str
         Color map file
      -e
         Export ILP
      -g
         Output search graph
      -l int
         Time limit in seconds for the ILP (default: -1, unlimited)
      -log
         Gurobi logging
      -m str
         Allowed migration patterns:
           0 : PS
           1 : PS, S
           2 : PS, S, M
           3 : PS, S, M, R
         If no pattern is specified, all allowed patterns will be
         enumerated (default: '0,1,2,3')
      -mutTreeIdx int
         Mutation tree index (default: -1)
      -noPR
         Disable polytomy resolution
      -o str
         Output prefix
      -p str
         Primary anatomical site
      -t int
         Number of threads (default: -1, #cores)

An example execution (executed from the root directory of the MACHINA repository):

    $ mkdir A7
    $ generatemutationtrees data/hoadley_2016/A7/A7_MACHINA_0.95.tsv > A7/mutation_trees.txt
    8/1/1/-1 (1)
    9/3/4/-1 (8)
    10/4/8/-1 (9)
    Found 2 mutation trees with 10 out of 10 mutations
    
    $ pmh_ti -p breast -c data/hoadley_2016/coloring.txt -m 1 -o A7 -F data/hoadley_2016/A7/A7_MACHINA_0.95.tsv \
    -barT A7/mutation_trees.txt > A7/result.txt
    Read 2 mutation trees.
    
    $ cat A7/result.txt
    0-      (PS, S) 5       5       2       mS      5146.83 5146.83 36.4647
    1-      (PS, S) 5       5       2       mS      5146.83 5146.83 37.8299

The program `generatemutationtrees` uses the SPRUCE algorithm to enumerate all mutation trees given a frequency matrix. The program `pmh_ti` considers solves the PMH-TI problem for each enumerated mutation tree. The `results.txt` file is formatted in exactly the same way as in `pmh`.

