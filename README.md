# Reed-Muller-Codes-Matlab
A MATLAB function library containing encoders, decoders and weight enumerators for Reed-Muller codes.

## Overview
Reed–Muller codes are a class of error-correcting codes that are used in wireless communications applications [1, 2].

Each Reed-Muller code is defined using two integer parameters, r <= m, and is notated RM(r,m). 
An RM(r,m) code has length 2^m, minimum distance 2^(m-r) and dimension \sum_{k=0}^{r} \nchoosek(m,k).

This Matlab library contains several encoders, unique decoders, list decoders and weight enumerators for Reed-Muller codes.

## User Information
### Dependencies
In order to use this project, you will need only Matlab. The library was tested on: Matlab R2020a.

The algorithmic code in this project mostly follows the notations in [6].

### User Manual
The main script is Tests.m. 
It contains testing procedures for all encoders and decoders implemented in the project. 
It can also be used as a reference for how you should call the different functions.

#### Encoders
1. rmenc_old.m, rmenc.m - encoding using a generator matrix.
2. rmenc_v2.m - encoding by multinomial evaluation.

#### Decoders
1. Reed's decoder - decodes any RM(r,m) using the majority-logic procedure [2]. implemented in rmdec_reed.m.
2. FHT decoder(s) - decodes RM(1,m) codes using Fast Hadamard Transform. 
   - The FHT unique decoder is implemented in rmdec_fht.m
   - The FHT list decoder is implemented in rmlistdec_fht.m.
3. Recursive decoder(s) - decodes any RM(r,m) by a recursive decomposition. , and t
   - The recursive unique decoder is implemented in rmdec_dumer.m [3].
   - The recursive list decoder is implemented in rmlistdec_dumer.m [4].
   - The recursive projection-aggregation unique decoder (RPA) is implemented in rmdec_rpa.m [5].
   - The recursive projection-aggregation list decoder (list RPA) is implemented in rmlistdec_rpa.m [5].

#### Demapping
rmdemap.m demaps an RM(r,m) codeword to its information bits using Reed's algorithm by means of polynomial evaluation.

#### Weight enumeration
reedmullerweights.m outputs the weights enumerator for any RM(1,m) and RM(2,m) [7] code as well as several specific RM(r,m) codes 
with r>2 and m<10 [8] (currently: RM(3,5), RM(3,6), RM(4,6), RM(3,7), RM(4,7), RM(3,8), RM(4,8), RM(4,9), RM(5,9)).

## References
[1] David E. Muller, "Application of Boolean algebra to switching circuit design and to error detection," Transactions of the IRE Professional Group on Electronic Computers, 1954.

[2] Irving S. Reed, "A class of multiple-error-correcting codes and the decoding scheme," Transactions of the IRE Professional Group on Information Theory, 1954.

[3] Ilya Dumer, "Soft-decision decoding of Reed-Muller codes: a simplified algorithm," IEEE transactions on Information Theory, 2006.

[4] Ilya Dumer and Kirill Shabunov, "Soft-decision decoding of Reed-Muller codes: recursive lists," IEEE Transactions on Information Theory, 2006.

[5] Min Ye and Emmanuel Abbe, "Recursive projection-aggregation decoding of Reed-Muller codes," IEEE Transactions on Information Theory, 2020.

[6] Emmanuel Abbe, Amir Shpilka and Min Ye, Reed-Muller Codes: Theory and Algorithms, IEEE Transactions on Information Theory, 2020.

[7] Neil J. A. Sloane and Elwyn R. Berlekamp, Weight Enumerator for Second-Order Reed-Muller Codes, IEEE Transactions on Information Theory, 1970.

[8] List of weight distributions, section 2.1: Reed-Muller codes, available: [OEIS list of weight distributions](http://oeis.org/wiki/List_of_weight_distributions#Reed-M.C3.BCller_codes).
