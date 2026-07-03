"""Utilities for loading the Magma-style examples shipped with this project."""

import re
from pathlib import Path

from sage.all import GF, matrix
from sage.misc.sage_eval import sage_eval

from modreps.representation import ModularRepresentation


_PRIME_FIELD_RE = re.compile(r"Finite Field of size\s+(\d+)")
_EXTENSION_FIELD_RE = re.compile(r"Finite Field in\s+(\w+)\s+of size\s+(\d+)\^(\d+)")


def load_magma_input(path, name=None):
    """Load the input representation from a Magma-style transcript.

    The parser understands the transcript format used in ``example1.txt`` to
    ``example5.txt``: a base-field line followed by generator matrices.
    """

    path = Path(path)
    text = path.read_text(encoding="utf-8")
    input_text = text.split("Output:", 1)[0]
    field, variable = _parse_base_field(input_text)
    matrices = _parse_matrices(input_text, field, variable)
    return ModularRepresentation(matrices, base_field=field, name=name)


def _parse_base_field(text):
    for line in text.splitlines():
        if "Base field:" not in line:
            continue
        descriptor = line.split("Base field:", 1)[1].strip()
        match = _EXTENSION_FIELD_RE.search(descriptor)
        if match:
            variable, prime, degree = match.groups()
            return GF(int(prime) ** int(degree), name=variable), variable
        match = _PRIME_FIELD_RE.search(descriptor)
        if match:
            return GF(int(match.group(1))), None
        raise ValueError(f"unsupported base field descriptor: {descriptor!r}")
    raise ValueError("could not find a base-field line")


def _parse_matrices(text, field, variable):
    matrices = []
    current = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            if current:
                matrices.append(_rows_to_matrix(current, field, variable))
                current = []
            continue
        if stripped.startswith("[") and stripped.endswith("]"):
            current.append(stripped[1:-1].strip())
    if current:
        matrices.append(_rows_to_matrix(current, field, variable))
    if not matrices:
        raise ValueError("no matrices found")
    return matrices


def _rows_to_matrix(rows, field, variable):
    nrows = len(rows)
    parsed_rows = [_parse_row(row, field, variable) for row in rows]
    widths = {len(row) for row in parsed_rows}
    if widths != {nrows}:
        raise ValueError(f"expected a square matrix, found row widths {sorted(widths)}")
    return matrix(field, parsed_rows)


def _parse_row(row, field, variable):
    tokens = row.split()
    entries = []
    i = 0
    while i < len(tokens):
        parts = [tokens[i]]
        i += 1
        while i < len(tokens) and tokens[i] == "+":
            if i + 1 >= len(tokens):
                raise ValueError(f"dangling plus sign in matrix row {row!r}")
            parts.extend([tokens[i], tokens[i + 1]])
            i += 2
        entries.append(_parse_entry(" ".join(parts), field, variable))
    return entries


def _parse_entry(entry, field, variable):
    entry = entry.replace("^", "**")
    locals_dict = {}
    if variable is not None:
        locals_dict[variable] = field.gen()
    return field(sage_eval(entry, locals=locals_dict))
