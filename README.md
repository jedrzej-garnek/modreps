# modreps

`modreps` is a small SageMath package for finite-dimensional modular
representations of matrix algebras over finite fields.

The convention matches Magma's matrix modules: vectors are row vectors and
generator matrices act on the right.  A representation is defined by a nonempty
list of square matrices over a common finite field.

## Installation

From this directory, inside a SageMath environment:

```bash
sage -python -m pip install -e .
```

If SageMath is installed in a conda environment named `sage`, use:

```bash
conda run -n sage sage -python -m pip install -e .
```

## Usage

```python
from sage.all import GF, matrix
from modreps import ModularRepresentation

k = GF(3)
A = matrix(k, [[1, 1], [0, 1]])
B = matrix(k, [[1, 0], [0, 1]])

M = ModularRepresentation([A, B])
M.dimension()
M.generators()
M.is_indecomposable()
M.indecomposable_summands()
```

This is analogous to the Magma workflow:

```magma
M := RModule(RSpace(GF(q), dimension), A);
L := IndecomposableSummands(M);
```

The matrices for summands are basis-dependent, so they need not match Magma's
displayed matrices.  The summand dimensions and isomorphism types are the stable
quantities.

## Running tests

```bash
conda run -n sage env DOT_SAGE=/tmp/codex-sage PYTHONPATH=src sage -python -m unittest discover -s tests
```

The tests include the supplied `example1.txt` to `example5.txt` transcripts and
assert the expected indecomposable summand dimensions.
