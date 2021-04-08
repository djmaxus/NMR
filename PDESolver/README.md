# PDE Solver
**Numerical solver of nonlinear diffusion for Ph.D. project**

## Author
|Contacts|Photo|
|-:|:-:|
|[Максим Елизарьев][djmaxus]<br>[elizarev@phystech.edu](mailto:elizarev@phystech.edu)|![](./djmaxus.jpeg)|
|||

## Contents
|                             File | Description|
|---------------------------------:|:-----------|
|     [solveDif](./solveDif.m)     | Numerical solver of nonlinear diffusion<br>* Fully-implicit two-layer scheme<br>* Newton – Gauss nonlinear solver<br>* ROM I-III functionality<br>* Calculation time logging|
|    [buildFlow](./buildFlow.m)    | Evaluate nonlinear flow derivatives|
|[buildJacobian](./buildJacobian.m)| Build Jacobian of spatial derivative out of nonlinear flow derivaties|
|[buildResidual](./buildResidual.m)| Build nonlinear residual|
|   [setProblem](./setProblem.m)   | Create data structure contating problem setup|
|      [setDiff](./setDiff.m)      | Create diffusion field|
|      [setInit](./setInit.m)      | Create initial state|
|      [plotSol](./plotSol.m)      | Plot results of `solveDif`|
|     [plotComp](./plotComp.m)     | Comparative plot of two `solveDif` results|
|  [definitions](./definitions.m)  | Some useful lambda-functions|
## References
* Carlberg, K., Bou-Mosleh, C., & Farhat, C. (2011). Efficient non-linear model reduction via a least-squares Petrov-Galerkin projection and compressive tensor approximations. International Journal for Numerical Methods in Engineering, 86(2), 155–181. https://doi.org/10.1002/nme.3050
* Everson, R., & Sirovich, L. (1995). Karhunen — Lòeve procedure for gappy data. J. Opt. Soc. Am. A, 12(8), 1657.

## Links
* https://djmaxus.github.io/djmaxus/
* http://databookuw.com/
* https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm

[//]: Links
[djmaxus]: https://djmaxus.github.io/djmaxus/
[email]:   mailto:elizariev.ma@cet-mipt.ru
[cet]:     https://www.cet-mipt.ru/
