
# ASR_DET_DC

This repository contains the code of the following paper.
- Yujie Li, Benying Tan, Atsunori Kanemura, Shuxue Ding, and Wuhui Chen, "Analysis sparse representation for nonnegative signals based on determinant measure by DC programming", _Complexity_, vol. 2018, 12 pp., April 2018. [[link]](https://www.hindawi.com/journals/complexity/2018/2685745/)

This software is licensed under the Apache License 2.0.


## Information

This software has been tested on the following envirnment.

* Windows 10
* MATLAB R2017b
* 7-Zip

This software uses the following programs and data.

* M. Yaghoobi, S. Nam, R. Gribonval, M. E. Davies, "Constrained overcomplete analysis operator learning for cosparse signal modelling," _IEEE Transaction on Signal Processing_, vol. 61, no. 9, pp 2341–2355, 2013.
    * http://www.mehrdadya.com/
* The Extended Yale Face Database B
    * http://vision.ucsd.edu/~iskwak/ExtYaleDatabase/ExtYaleB.html
* KSVD Box v9
    * http://www.cs.technion.ac.il/~ronrubin/software.html


## How to run

1. Run `setup.ps1` with Windows PowerShell.
2. Expand `CAOLv1.tar.gz` with 7-Zip.
3. Copy the following file to current working directory from the expanded `CAOLv1.tar.gz`.
    * `ProjUNColBall.m`
4. Open `main.m` with MATLAB.
5. Run `main`.
