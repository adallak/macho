# Introduction
The repository contains an algorithm to implement the Majorized  Cholesky (MaCho) method as described in <https://arxiv.org/abs/2404.05148>.
It contains three main files.

- `macho.r` - contains the main `macho` function to estimate a topological ordering and the `DAG_from_Ordering` function to estimate a DAG from
  ordering. The latter function is borrowed from  <https://github.com/WY-Chen/EqVarDAG/>
- `experiments.r`- contains part of the experiments described in the paper.
- `utility.r` - utility file for generating data and measuring the performance.
