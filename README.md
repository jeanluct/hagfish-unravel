# hagfish-unravel

This repository contains Matlab code used in the paper [_Unraveling
hagfish slime_][1] by Gaurav Chaudhary, Randy Ewoldt, and Jean-Luc
Thiffeault.

There are 5 different cases considered, and for each there is a Matlab
function `Lsolve_<case>.m` to compute the solution for the unraveled
thread length `L(t)`, and a corresponding script
`plot_Lsolve_<case>.m` that calls the function and plots the results.
The 5 cases are

* `Lsolve_pinned_thread.m`: Section II.A of the paper.

* `Lsolve_pinned_skein.m`: Section II.B of the paper.

* `Lsolve_freefree.m`: Section II.C of the paper.

* `Lsolve_two_skeins.m`: Section II.D of the paper.

* `Lsolve_pinned_skein_fishflow.m`: Section III of the Supplementary
  Information.

[1]: http://arxiv.org/abs/1809.xxxxx
