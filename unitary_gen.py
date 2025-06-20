"""
Generate the Clifford group for 1 qubit via composition of H and S unitary matrices until closure
"""

import numpy as np


class GroupElement:
    def __init__(self, name, unitary):
        self.name = name
        self.unitary = unitary


I_mat = GroupElement("I", np.array([[1, 0], [0, 1]], dtype=complex))
S = GroupElement("S", np.array([[1, 0], [0, 1j]], dtype=complex))
H = GroupElement("H", np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2))

generators = [H, S]


def equivalent_up_to_global_phase(U1, U2, atol=1e-5):
    product = (
        U1 @ U2.conj().T
    )  # if equiv up to global phase, the product should be a scalar multiple of the identity
    global_phase = product[0, 0]  # diagonal elements will be the global phase factor

    # Zero case
    if np.abs(global_phase) < 1e-8:
        return False

    normalized = product / global_phase  # this will be identity if equiv up to phase
    return np.allclose(
        normalized, np.eye(U1.shape[0]), atol=atol
    )  # compare to identity


def generate_clifford_group(generators, global_phase=False):
    clifford_list = [I_mat]
    queue = [I_mat]

    while queue:
        current = queue.pop(0)
        for gen in generators:
            new_element_name = (
                (gen.name + current.name) if (current.name != "I") else gen.name
            )
            new_element_unitary = gen.unitary @ current.unitary

            if global_phase:
                # Check nothing close including global phase
                is_new = not any(
                    np.allclose(existing.unitary, new_element_unitary, atol=1e-5)
                    for existing in clifford_list
                )
            else:
                # Check nothing close excluding global phase
                is_new = not any(
                    equivalent_up_to_global_phase(existing.unitary, new_element_unitary)
                    for existing in clifford_list
                )

            if is_new:
                new_element = GroupElement(new_element_name, new_element_unitary)
                clifford_list.append(new_element)
                queue.append(new_element)

    return clifford_list


def is_associative(group):
    """
    Check all triples in the group for associativity
    """
    results = []

    for a in group:
        for b in group:
            for c in group:
                left = a.unitary @ (b.unitary @ c.unitary)
                right = (a.unitary @ b.unitary) @ c.unitary
                results.append(np.allclose(left, right, atol=1e-5))

    return all(results)


def make_cayley_table(group, n):
    """
    Make a cayley table; only shows group element composition up to global phase
    """
    # Use first n elements for rows and columns
    n = min(n, len(group))
    table = [["" for _ in range(n)] for _ in range(n)]

    unitaries = [el.unitary for el in group]

    for i in range(n):
        for j in range(n):
            product = unitaries[i] @ unitaries[j]

            # Reduce to unique element name (in S, H) by global phase equivalence
            found_name = None
            for k in range(len(group)):
                if equivalent_up_to_global_phase(product, unitaries[k]):
                    found_name = group[k].name
                    break

            if found_name is None:
                raise ValueError(
                    f"Product of {group[i].name} and {group[j].name} not found in group"
                )

            table[i][j] = found_name

    return table


def print_cayley_table(group, table):
    n = len(table)
    names = [el.name for el in group[:n]]

    max_width = max(
        max(len(name) for name in names),
        max(len(cell) for row in table for cell in row),
    )

    header = " " * (max_width + 3) + " ".join(f"{name:>{max_width}}" for name in names)
    print(header)
    print("-" * len(header))

    # Print each row with row label
    for i, row in enumerate(table):
        row_str = f"{names[i]:>{max_width}} | " + " ".join(
            f"{cell:>{max_width}}" for cell in row
        )
        print(row_str)


if __name__ == "__main__":
    print("Generating Clifford group on 1 qubit up to global phases")
    group_up_to_phases = generate_clifford_group(generators, global_phase=True)
    print("Group generated:")
    names = [element.name for element in group_up_to_phases]
    print(names)
    print("Order:")
    print(len(names))
    print("Checking associativity for this group")
    if is_associative(group_up_to_phases):
        print("It is associative ✅ ")
    else:
        print("Not associative ❌ ")
    print("Generating Cayley table for this group")
    table = make_cayley_table(group_up_to_phases, 5)
    print_cayley_table(group_up_to_phases, table)
