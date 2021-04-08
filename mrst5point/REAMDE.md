# `mrst5point`
**[MRST][mrst]-based 5-point oil recovery model**

Requires [Matlab Reservoir Simulation Toolbox][mrst] installed on your computer.

## Author
|Contacts|Photo|
|-:|:-:|
|[Максим Елизарьев][djmaxus]<br>[elizarev@phystech.edu](mailto:elizarev@phystech.edu)|![](./djmaxus.jpeg)|
|||

## Contents
| File | Description |
|-:|:-|
|[startMRSTmodules.m](./startMRSTmodules.m)| Initializes MRST|
|[startup_user.m](./startup_user.m)| Custom MRST startup<br>_used in_ `startMRSTmodules.m`|
|[fivePoint.m](./fivePoint.m)|Creates data structure of problem parameters<br>Takes required permeability field as an input|
|[ResSimSPCA.m](./ResSimSPCA.m)|Takes problem parameters and performs simulation|
|[run5point.m](./run5point.m)|Combines two above functions<br>Takes permeability field and perform simulation|
|[FOM.m](./FOM.m)| **F**ull **O**rder **M**odel class<br>* performs simulations<br>* stores data<br>* calculates history matching objective function|


## Links
* https://djmaxus.github.io/djmaxus/
* https://www.sintef.no/projectweb/mrst/

[//]: Shorts
[djmaxus]: https://djmaxus.github.io/djmaxus/
[email]:   mailto:elizariev.ma@cet-mipt.ru
[cet]:     https://www.cet-mipt.ru/
[mrst]: https://www.sintef.no/projectweb/mrst/
