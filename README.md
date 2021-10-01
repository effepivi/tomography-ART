# tomography-ART
This repository contains the serial C++ code to reconstruct X-ray tomographic data using the algebraic reconstruction technique (ART). It is an iterative reconstruction algorithm.

Example of command line:
```bash
$ ./myART -s ../data/shepp_logan/sinogram.tif -r ../data/shepp_logan/reference_CT.tif  -o art_reconstruction.tif -i 200
```
