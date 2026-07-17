# `multicharge/snippets` contains small codes to demonstrate an approac

## Index generation

The original code in `subroutine get_amat_0d` loops over a triangle
of matrix elements. The code is parallelised with OpenMP. The initial
code involves 2 loops, but OpenMP can only parallelise 1 of those
because the pair don't loop over a rectangle (or put otherwise the
loop limits of the inner loop depend on iteration number in the
outer loop). If we want to increase the amount of parallelism that
can be obtained we have to do something about the loop structure.

The 3 codes here are:
1. `prog_triangle.f90` generates the indeces analogously to the original
   code in `get_amat_0d` by looping over the upper triangle of a square
   matrix.
2. `prog_single.f90` replaces the 2 loops running over the upper triangle
   with 1 loop that runs over as many elements as there are in the upper
   triangle. The indeces into the matrix are computed from the loop
   counter. This means the loop can be arbitrarily broken while still
   generating the same index sequence.
3. `prog_mapping.f90` replaces the loops over a triangle with loops
   over a rectangle the size of the upper half of the square matrix. The
   part of the rectangle that overlaps with the upper triangle is kept
   as is. The part that overlaps with the lower triangle is mapped into
   the lower part of the upper triangle that the rectangle does not cover.

All 3 codes take the matrix dimension as a single command line argument.
In particular the last program generates indeces in a different order than
the original code. Using `sort` the output can reordered after which the outputs
from all 3 codes should be identical.
