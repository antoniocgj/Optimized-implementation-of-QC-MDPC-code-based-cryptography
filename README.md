# Optimized-implementation-of-QC-MDPC-code-based-cryptography

This repository holds two optimized implementations of a QC-MDPC code-based KEM. Both versions use algorithms from [Chou's QcBits](https://www.win.tue.nl/~tchou/qcbits/). 

## Updated_QcBits 

Paper: [Guimarães, A, Aranha, DF, Borin, E. Optimized implementation of QC‐MDPC code‐based cryptography. Concurrency Computat Pract Exper. 2019; 31:e5089. ](https://doi.org/10.1002/cpe.5089)

**Abstract** This paper presents a new enhanced version of the QcBits key encapsulation mechanism, which is a constant‐time implementation of the Niederreiter cryptosystem using QC‐MDPC codes. In this version, we updated the implementation parameters to meet the 128‐bit quantum security level, replaced some of the core algorithms to avoid using slower instructions, vectorized the entire code using the AVX‐512 instruction set extension, and applied several other techniques to achieve a competitive performance level. Our implementation takes 928, 259, and 5008 thousand Skylake cycles to perform batch key generation (cost per key), encryption, and uniform decryption, respectively. Comparing with the current state‐of‐the‐art implementation for QC‐MDPC codes, BIKE, our code is 1.9 times faster when decrypting messages.

## Updated_QcBits_AFR

Paper: [Guimaraes, A., Borin, E., & Aranha, D. F. (2019, May). Introducing arithmetic failures to accelerate QC-MDPC code-based cryptography. In Code-Based Cryptography Workshop (pp. 44-68). Springer, Cham.](https://doi.org/10.1007/978-3-030-25922-8_3)

**Abstract** In this work, we optimize the performance of QC-MDPC code-based cryptosystems through the insertion of configurable failure rates in their arithmetic procedures. We present constant time algorithms with a configurable failure rate for multiplication and inversion over binary polynomials, the two most expensive subroutines used in QC-MDPC implementations. Using a failure rate negligible compared to the security level (2^{−128}), our multiplication is 2 times faster than NTL on sparse polynomials and 1.6 times faster than a naive constant-time sparse polynomial multiplication. Our inversion algorithm, based on Wu et al., is 2 times faster than the original algorithm and 12 times faster than Itoh-Tsujii using the same modulus polynomial (x^{32749} − 1). By inserting these algorithms in a version of QcBits at the 128-bit quantum security level, we were able to achieve a speedup of 1.9 on the key generation and up to 1.4 on the decryption time. Comparing with variant 2 of the BIKE suite, which also implements the Niederreiter Cryptosystem using QC-MDPC codes, our final version of QcBits performs the uniform decryption 2.7 times faster.

## Execution time

Execution time (cycles) of our updated versions of QcBits, with (Updated_QcBits_AFR) and without (Updated_QcBits) a negligible (< 2^{-128}) Arithmetic Failure Rate; and of the additional implementation of BIKE-2, for comparison.

|                          | Updated_QcBits | Updated_QcBits_AFR | BIKE-2         |
|--------------------------|:--------------:|:------------------:|----------------|
| Key Generation           |   40,265,904   |     21,332,058     | 12,944,920 (*) |
| Encryption               |     259,306    |       256,655      | 348,227        |
| Constant-time Decryption |    9,803,835   |      8,016,312     | (**)           |
| Uniform Decryption       |    5,008,429   |      3,505,079     | 9,572,412      |

(∗) BIKE’s polynomial inversion is not constant-time.
(**) BIKE does not present constant-time decryption

We benchmarked the implementations on an Intel i7-7820X processor with Hyper-Threading and TurboBoost disabled. The execution time of BIKE (Variant 2) is presented for comparison. More details in the paper.


## Build 

### Constant-time decryption

```console
cd updated_qcbits/ 
make
./qcbits
```

### Uniform Decryption 

```console
cd updated_qcbits/ 
make
./qcbits
```