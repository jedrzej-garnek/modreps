"""Representations of finitely generated matrix algebras over finite fields.

The convention in this module follows Magma's matrix-module convention: vectors
are row vectors and generator matrices act on the right.
"""

from itertools import islice, product

from sage.all import (
    MatrixSpace,
    VectorSpace,
    block_diagonal_matrix,
    matrix,
    set_random_seed,
)
from sage.structure.sage_object import SageObject


def _field_label(field):
    """Return a compact ``GF(q)`` label for a Sage finite field."""

    p = field.characteristic()
    d = field.degree()
    if d == 1:
        return f"GF({p})"
    return f"GF({p}^{d})"


def _is_scalar_matrix(mat):
    """Return whether ``mat`` is a scalar matrix."""

    if mat.nrows() != mat.ncols():
        return False
    n = mat.nrows()
    if n == 0:
        return True
    scalar = mat[0, 0]
    for i in range(n):
        for j in range(n):
            if mat[i, j] != (scalar if i == j else 0):
                return False
    return True


def _matrix_from_vector(field, nrows, ncols, vector):
    """Convert a flattened row-major vector to a matrix."""

    return matrix(field, nrows, ncols, list(vector))


def _product_power_size(base, exponent, limit):
    """Return whether ``base ** exponent`` is at most ``limit``.

    The computation stops once the partial product exceeds ``limit``.
    """

    size = 1
    for _ in range(exponent):
        size *= base
        if size > limit:
            return False
    return True


class ModularRepresentation(SageObject):
    r"""A finite-dimensional module for a matrix algebra over a finite field.

    INPUT:

    - ``generators`` -- nonempty iterable of square matrices of the same size.
      The matrices are the right actions of the algebra generators.
    - ``base_field`` -- optional finite field used to coerce ``generators``.
    - ``name`` -- optional display name.

    EXAMPLES::

        sage: from modreps import ModularRepresentation
        sage: k = GF(3)
        sage: A = matrix(k, [[1, 1], [0, 1]])
        sage: M = ModularRepresentation([A])
        sage: M.dimension()
        2
        sage: M.generators()[0] == A
        True
    """

    def __init__(self, generators, base_field=None, name=None, check=True):
        generators = list(generators)
        if not generators:
            raise ValueError("at least one generator matrix is required")

        if base_field is None:
            try:
                base_field = generators[0].base_ring()
            except AttributeError as exc:
                raise TypeError(
                    "base_field is required when generators are not Sage matrices"
                ) from exc

        if not base_field.is_field() or not base_field.is_finite():
            raise ValueError("base_field must be a finite field")

        converted = []
        dim = None
        for gen in generators:
            mat = matrix(base_field, gen)
            if not mat.is_square():
                raise ValueError("generator matrices must be square")
            if dim is None:
                dim = mat.nrows()
            elif mat.nrows() != dim:
                raise ValueError("all generator matrices must have the same size")
            if check and mat.base_ring() != base_field:
                raise ValueError("could not coerce a generator to the base field")
            mat.set_immutable()
            converted.append(mat)

        self._field = base_field
        self._dim = dim
        self._generators = tuple(converted)
        self._name = name
        self._endomorphism_basis = None

    def _repr_(self):
        label = f"RModule of dimension {self._dim} over {_field_label(self._field)}"
        if self._name:
            field = _field_label(self._field)
            return f"RModule {self._name} of dimension {self._dim} over {field}"
        return label

    def __eq__(self, other):
        if not isinstance(other, ModularRepresentation):
            return False
        return (
            self._field == other._field
            and self._generators == other._generators
            and self._name == other._name
        )

    def base_field(self):
        """Return the finite base field."""

        return self._field

    def dimension(self):
        """Return the vector-space dimension."""

        return self._dim

    def ngens(self):
        """Return the number of algebra generators."""

        return len(self._generators)

    def generators(self):
        """Return the defining generator matrices."""

        return self._generators

    gens = generators

    def generator(self, i):
        """Return the ``i``-th defining generator matrix."""

        return self._generators[i]

    def homomorphism_basis(self, other):
        r"""Return a basis for ``Hom(self, other)``.

        Homomorphisms are matrices ``X`` satisfying ``A_i X = X B_i`` for all
        paired generators ``A_i`` of ``self`` and ``B_i`` of ``other``.
        """

        self._check_compatible_target(other)
        return self._homomorphism_basis_to(other)

    def endomorphism_basis(self):
        """Return a basis for the endomorphism algebra of this representation."""

        if self._endomorphism_basis is None:
            self._endomorphism_basis = self._homomorphism_basis_to(self)
        return self._endomorphism_basis

    def is_indecomposable(self, max_tries=200, exhaustive_limit=20000, seed=None):
        """Return whether no nontrivial direct decomposition was found.

        For small endomorphism rings this is exhaustive.  For larger rings the
        method uses deterministic candidates followed by random endomorphisms;
        increase ``max_tries`` for difficult examples.
        """

        return self._find_direct_split(max_tries, exhaustive_limit, seed) is None

    def indecomposable_summands(
        self, max_tries=200, exhaustive_limit=20000, seed=None, sort=True
    ):
        r"""Return a list of indecomposable direct summands.

        The summands are again :class:`ModularRepresentation` objects.  Their
        generator matrices are written in bases found during the decomposition,
        so they need not match Magma's displayed matrices.
        """

        split = self._find_direct_split(max_tries, exhaustive_limit, seed)
        if split is None:
            return [self]

        left, right = split
        summands = []
        for subspace in (left, right):
            submodule = self._subrepresentation(subspace)
            summands.extend(
                submodule.indecomposable_summands(
                    max_tries=max_tries,
                    exhaustive_limit=exhaustive_limit,
                    seed=seed,
                    sort=False,
                )
            )
        if sort:
            summands.sort(key=lambda rep: rep.dimension())
        return summands

    def is_isomorphic(
        self,
        other,
        max_tries=200,
        exhaustive_limit=20000,
        seed=None,
        certificate=False,
    ):
        """Test whether two representations are isomorphic.

        If ``certificate`` is true, return ``(answer, X)`` where ``X`` is an
        intertwining matrix when an isomorphism is found.
        """

        if not isinstance(other, ModularRepresentation):
            raise TypeError("other must be a ModularRepresentation")
        if not self._has_same_signature(other):
            return (False, None) if certificate else False

        basis = self.homomorphism_basis(other)
        iso = self._find_invertible_linear_combination(
            basis, max_tries=max_tries, exhaustive_limit=exhaustive_limit, seed=seed
        )
        answer = iso is not None
        if certificate:
            return answer, iso
        return answer

    def change_basis(self, basis_matrix, name=None):
        r"""Return the same representation in a new row-vector basis.

        If the rows of ``basis_matrix`` are the new basis vectors in the old
        coordinates, each generator ``A`` is replaced by ``P A P^{-1}``.
        """

        basis_matrix = matrix(self._field, basis_matrix)
        if basis_matrix.nrows() != self._dim or basis_matrix.ncols() != self._dim:
            raise ValueError("basis_matrix must be square of size dimension()")
        if basis_matrix.det() == 0:
            raise ValueError("basis_matrix must be invertible")
        inverse = basis_matrix.inverse()
        return ModularRepresentation(
            [basis_matrix * gen * inverse for gen in self._generators],
            base_field=self._field,
            name=name,
        )

    @classmethod
    def direct_sum(cls, representations, name=None):
        """Return the direct sum of compatible representations."""

        representations = list(representations)
        if not representations:
            raise ValueError("at least one representation is required")

        field = representations[0].base_field()
        ngens = representations[0].ngens()
        for rep in representations:
            if not isinstance(rep, cls):
                raise TypeError("all summands must be ModularRepresentation objects")
            if rep.base_field() != field:
                raise ValueError("all summands must have the same base field")
            if rep.ngens() != ngens:
                raise ValueError("all summands must have the same number of generators")

        generators = []
        for i in range(ngens):
            generators.append(
                block_diagonal_matrix([rep.generator(i) for rep in representations])
            )
        return cls(generators, base_field=field, name=name)

    def direct_sum_with(self, *others, name=None):
        """Return the direct sum of ``self`` and ``others``."""

        return self.direct_sum((self,) + others, name=name)

    def _has_same_signature(self, other):
        return (
            self._field == other._field
            and self._dim == other._dim
            and self.ngens() == other.ngens()
        )

    def _check_compatible_target(self, other):
        if not isinstance(other, ModularRepresentation):
            raise TypeError("target must be a ModularRepresentation")
        if self._field != other._field:
            raise ValueError("representations must have the same base field")
        if self.ngens() != other.ngens():
            raise ValueError("representations must have the same number of generators")

    def _homomorphism_basis_to(self, other):
        field = self._field
        source_dim = self._dim
        target_dim = other._dim
        nvars = source_dim * target_dim
        rows = []

        for source_gen, target_gen in zip(self._generators, other._generators):
            for r in range(source_dim):
                for c in range(target_dim):
                    equation = [field.zero()] * nvars
                    for t in range(source_dim):
                        equation[t * target_dim + c] += source_gen[r, t]
                    for s in range(target_dim):
                        equation[r * target_dim + s] -= target_gen[s, c]
                    rows.append(equation)

        if rows:
            equations = matrix(field, rows)
            kernel = equations.right_kernel()
        else:
            kernel = VectorSpace(field, nvars)

        basis = tuple(
            _matrix_from_vector(field, source_dim, target_dim, vector)
            for vector in kernel.basis()
        )
        for mat in basis:
            mat.set_immutable()
        return basis

    def _find_direct_split(self, max_tries, exhaustive_limit, seed):
        if self._dim <= 1:
            return None

        end_basis = self.endomorphism_basis()
        if len(end_basis) == 1 and _is_scalar_matrix(end_basis[0]):
            return None

        if seed is not None:
            set_random_seed(seed)

        for candidate in self._deterministic_endomorphism_candidates(
            end_basis, exhaustive_limit
        ):
            split = self._split_from_endomorphism(candidate)
            if split is not None:
                return split

        for _ in range(max_tries):
            candidate = self._random_linear_combination(end_basis)
            split = self._split_from_endomorphism(candidate)
            if split is not None:
                return split
        return None

    def _deterministic_endomorphism_candidates(self, basis, exhaustive_limit):
        field = self._field
        q = int(field.order())
        d = len(basis)

        for mat in basis:
            if not _is_scalar_matrix(mat):
                yield mat

        for i in range(d):
            for j in range(i + 1, d):
                candidate = basis[i] + basis[j]
                candidate.set_immutable()
                if not _is_scalar_matrix(candidate):
                    yield candidate

        if _product_power_size(q, d, exhaustive_limit):
            field_elements = list(field)
            for coefficients in product(field_elements, repeat=d):
                if all(coef == 0 for coef in coefficients):
                    continue
                candidate = self._linear_combination(basis, coefficients)
                if not _is_scalar_matrix(candidate):
                    yield candidate
        else:
            field_elements = list(islice((a for a in field if a != 0), 3))
            for i in range(d):
                for j in range(i + 1, d):
                    for coef in field_elements:
                        candidate = basis[i] + coef * basis[j]
                        candidate.set_immutable()
                        if not _is_scalar_matrix(candidate):
                            yield candidate

    def _linear_combination(self, basis, coefficients):
        space = MatrixSpace(self._field, basis[0].nrows(), basis[0].ncols())
        result = space.zero()
        for coefficient, mat in zip(coefficients, basis):
            if coefficient:
                result += coefficient * mat
        result.set_immutable()
        return result

    def _random_linear_combination(self, basis):
        coefficients = [self._field.random_element() for _ in basis]
        if all(coef == 0 for coef in coefficients):
            coefficients[0] = self._field.one()
        return self._linear_combination(basis, coefficients)

    def _split_from_endomorphism(self, endomorphism):
        if _is_scalar_matrix(endomorphism):
            return None

        split = self._split_from_minimal_polynomial(endomorphism)
        if split is not None:
            return split

        power = endomorphism**self._dim
        rank = power.rank()
        if 0 < rank < self._dim:
            kernel = power.left_kernel()
            image = power.row_space()
            if self._are_complementary(kernel, image):
                return kernel, image
        return None

    def _split_from_minimal_polynomial(self, endomorphism):
        minpoly = endomorphism.minimal_polynomial()
        factors = list(minpoly.factor())
        if len(factors) <= 1:
            return None

        for factor, exponent in factors:
            primary = factor(endomorphism) ** exponent
            left = primary.left_kernel()
            dim_left = left.dimension()
            if dim_left == 0 or dim_left == self._dim:
                continue
            quotient = minpoly // (factor**exponent)
            right = quotient(endomorphism).left_kernel()
            if self._are_complementary(left, right):
                return left, right
        return None

    def _are_complementary(self, left, right):
        return (
            left.dimension() > 0
            and right.dimension() > 0
            and left.dimension() + right.dimension() == self._dim
            and left.intersection(right).dimension() == 0
        )

    def _subrepresentation(self, subspace):
        field = self._field
        basis = list(subspace.basis())
        dim = len(basis)
        if dim == 0 or dim == self._dim:
            raise ValueError("subspace must be nonzero and proper")

        restricted_generators = []
        for gen in self._generators:
            rows = []
            for row in basis:
                image = row * gen
                rows.append(list(subspace.coordinate_vector(image)))
            mat = matrix(field, dim, dim, rows)
            mat.set_immutable()
            restricted_generators.append(mat)

        return ModularRepresentation(restricted_generators, base_field=field)

    def _find_invertible_linear_combination(
        self, basis, max_tries, exhaustive_limit, seed
    ):
        if not basis:
            return None

        for candidate in basis:
            if candidate.is_square() and candidate.det() != 0:
                return candidate

        field = self._field
        q = int(field.order())
        d = len(basis)

        if _product_power_size(q, d, exhaustive_limit):
            for coefficients in product(list(field), repeat=d):
                if all(coef == 0 for coef in coefficients):
                    continue
                candidate = self._linear_combination(basis, coefficients)
                if candidate.det() != 0:
                    return candidate

        if seed is not None:
            set_random_seed(seed)

        for _ in range(max_tries):
            candidate = self._random_linear_combination(basis)
            if candidate.det() != 0:
                return candidate
        return None
