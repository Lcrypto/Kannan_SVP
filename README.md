# Kannan method for solving Short Vector Problem and Block Korkin-Zolotarev method for solve Short Basis Problem


![alt text](https://github.com/Lcrypto/Kannan_SVP/blob/master/Ideal%20Lattice%20%20Challenge%20(TU%20Darmstadt-U%20Wollongong)%202013%20result.jpg)


https://web.archive.org/web/20130921071646/https:/www.latticechallenge.org/ideallattice-challenge/index.php

The GitHub repository contains implementations of Kannan's methods for enumerating the shortest vector in a lattice with a norm less than some radius, also known as the Short(est) Vector Problem (SVP). The provided implementations include solutions in MATLAB and C++ using the NTL Library by Victor Shoup.

In addition to Kannan's method, the repository includes various lattice reduction algorithms such as LLL, KZ/HKZ, Seysen, and Brun. These algorithms can preprocess the lattice to speed up the enumeration process.

The base algorithm used for Sphere Decoder (MIMO, QR-MLD, Lattice Aided), ECC decoding (BDD), Hamming distance estimation (based on Kannan Embedding techniques), and Post Quantum Cryptography cryptoanalysis (before LWE-problems) can be found in the repository. However, it should be noted that these methods may not be directly applicable to non-binary or non-ternary cases.



One of important application of Number geometry it estimation of code distance using  Lattice-based method (Kannan embeding, SVP, SBP, Block Korkin-Zolotarev, BKZ for solution Shortest Basis Problem). According to our (Usatyuk Vasiliy) results in the code distance challenge at https://decodingchallenge.org/low-weight/, Lattice methods are superior:
![alt text](https://github.com/Lcrypto/Length-und-Rate-adaptive-code/blob/master/Code_distance_challenge.png)

Overall, this repository provides useful tools and implementations for solving various problems related to lattices, including SVP, short(est) basis problems (SBP), Hamming distance estimation, and cryptography-related analyses.
