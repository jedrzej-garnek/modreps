import unittest
from sage.all import GF, identity_matrix, matrix

from modreps import ModularRepresentation


class ModularRepresentationTests(unittest.TestCase):
    def test_basic_accessors_and_repr(self):
        field = GF(3)
        generator = matrix(field, [[1, 1], [0, 1]])
        module = ModularRepresentation([generator], name="M")

        self.assertEqual(module.base_field(), field)
        self.assertEqual(module.dimension(), 2)
        self.assertEqual(module.ngens(), 1)
        self.assertEqual(module.generators(), (generator,))
        self.assertEqual(repr(module), "RModule M of dimension 2 over GF(3)")

    def test_direct_sum_decomposes_into_expected_summands(self):
        field = GF(3)
        simple = ModularRepresentation([matrix(field, [[1]])])
        jordan = ModularRepresentation([matrix(field, [[1, 1], [0, 1]])])
        module = ModularRepresentation.direct_sum([simple, jordan])

        self.assertEqual(module.dimension(), 3)
        self.assertEqual(
            sorted(s.dimension() for s in module.indecomposable_summands(seed=2)),
            [1, 2],
        )
        self.assertFalse(module.is_indecomposable(seed=2))
        self.assertTrue(simple.is_indecomposable(seed=2))
        self.assertTrue(jordan.is_indecomposable(seed=2))

    def test_change_basis_produces_isomorphic_module(self):
        field = GF(5)
        generators = [
            matrix(field, [[1, 1, 0], [0, 1, 0], [0, 0, 2]]),
            matrix(field, [[2, 0, 0], [0, 2, 1], [0, 0, 2]]),
        ]
        module = ModularRepresentation(generators)
        basis = matrix(field, [[1, 1, 0], [0, 1, 1], [1, 0, 1]])
        changed = module.change_basis(basis)

        isomorphic, certificate = module.is_isomorphic(
            changed, certificate=True, seed=3
        )

        self.assertTrue(isomorphic)
        self.assertNotEqual(certificate.det(), 0)
        for source, target in zip(module.generators(), changed.generators()):
            self.assertEqual(source * certificate, certificate * target)

    def test_isomorphism_rejects_incompatible_modules(self):
        field = GF(3)
        one = ModularRepresentation([identity_matrix(field, 1)])
        two = ModularRepresentation([identity_matrix(field, 2)])

        self.assertFalse(one.is_isomorphic(two))

    def test_constructor_validation(self):
        field = GF(3)

        with self.assertRaisesRegex(ValueError, "at least one"):
            ModularRepresentation([], base_field=field)

        with self.assertRaisesRegex(ValueError, "square"):
            ModularRepresentation([matrix(field, 1, 2, [1, 0])])

        with self.assertRaisesRegex(ValueError, "same size"):
            ModularRepresentation(
                [identity_matrix(field, 1), identity_matrix(field, 2)]
            )
