# computingPicard
This repository has practical implementations in SAGE and PARI/GP which compute the L-polynomials of generic Picard curves using the algorithms in https://arxiv.org/pdf/2010.07247.pdf. 

We provide implementations in Sage and in Pari/GP. We recommend that you use the Pari/GP implementation if you are interesting in doing a large computation, as it is significantly faster.

To use any of the algorithms, a file with Cartier-Manin matrices needs to be inputted. These are computed using the algorithm in https://arxiv.org/abs/2004.10189. We provide some sample outputs of this algorithm to use in the data folder. When using these lists of Cartier-Manin matrices, for the Sage implementation you save them as a text file and input the path. For Pari/GP you use readstr() on the text file to input it. 

# Sage

For the Sage implementation, an implementation of Jacobian arithmetic is used (from https://github.com/emmacneil/c34-curves, written by Evan Macneil). Under the c34 curves folder, we have put Evan Macneil's Sage files. These should be put in the same folder as lpoly.sage in order to be loaded in.

As an example use of the Sage implementation, you might do the following in a Jupyter notebook:

import os
cwd = os.getcwd()
os.chdir(cwd)

cm_path = cwd + "/hwlpolys_1.txt"

findLpolys(cm_path, 1000, 10000, True)

This would print the L-polynomials for primes between 1000 and 10000 for y^3=x^4+x+1.

# Pari/GP

For the Pari/GP implementation, no additional files are needed. To use it, first save the Cartier Manin matrix text file and FrobTraces.txt file (if you intend to use it). Then, set CMList to be readstr() applied to your saved Cartier Manin matrix text file (e.g. hwlpolys_1.txt). The Frobenius traces in FrobTraces allow lpoly.gp to compute a cube root of unity modulo p deterministically, making the algorithm completely deterministic. This does speed up the algorithm, but the effect is not significant.

If you decide not to use the Frobenius traces, run liftLpolysAlt(CMList,C,i) on the intended range. Here, the list C = [C[1],C[2],C[3]] corresponds to the Picard curve y^3 = x^4+C[3]x^2+C[2]x+C[1]. The index i simply reads off the corresponding element in CMList, and computes the corresponding L-polynomial for C modulo p.

If you decide to use the Frobenius traces, you can directly input the array in FrobTraces (under data) into Pari/GP as TrList. Then liftLpolys(CMList,TrList,C,i) with this additional input again computes the L-polynomial, but now uses the corresponding trace of Frobenius to do it deterministically. Depending on the first prime in CMList, you will want to shift indices in TrList so that these match up.