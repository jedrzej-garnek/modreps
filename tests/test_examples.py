from pathlib import Path
import unittest

from modreps import ModularRepresentation, load_magma_input


ROOT = Path(__file__).resolve().parents[1]


class ExampleTests(unittest.TestCase):
    def test_magma_examples_have_expected_summand_dimensions(self):
        expected = {
            "example1.txt": [1, 3, 6],
            "example2.txt": [6, 14],
            "example3.txt": [12],
            "example4.txt": [1, 3, 6],
            "example5.txt": [2, 9, 9],
        }

        for filename, dimensions in expected.items():
            with self.subTest(filename=filename):
                module = load_magma_input(ROOT / filename)
                summands = module.indecomposable_summands(seed=1, max_tries=500)

                self.assertIsInstance(module, ModularRepresentation)
                self.assertEqual(
                    sorted(summand.dimension() for summand in summands), dimensions
                )
                self.assertEqual(
                    sum(summand.dimension() for summand in summands),
                    module.dimension(),
                )
