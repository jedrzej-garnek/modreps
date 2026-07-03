Examples
========

Constructing a representation is similar to Magma's ``RModule`` command::

   from sage.all import GF, matrix
   from modreps import ModularRepresentation

   k = GF(3)
   A = matrix(k, [[1, 1], [0, 1]])
   B = matrix(k, [[1, 0], [0, 1]])
   M = ModularRepresentation([A, B])

The main methods mirror the common module operations::

   M.dimension()
   M.generators()
   M.is_indecomposable()
   summands = M.indecomposable_summands()

Summand matrices are written in bases found by the decomposition algorithm, so
the exact displayed matrices can differ from Magma.  The dimensions and
isomorphism classes of the summands are invariant.

The package also includes a loader for the transcript format used by the
provided example files::

   from modreps import load_magma_input

   M = load_magma_input("example1.txt")
   [N.dimension() for N in M.indecomposable_summands()]
