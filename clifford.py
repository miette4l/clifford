"""
Generation of the order-24 Clifford group C/U(1) on 1 qubit by Clifford tableau multiplication.
Demonstration of the nontrivial split extension:
        1 → Z₂xZ₂ → C/U(1) → Sp(2, Z₂) ≅ S₃ → 1     i.e. C/U(1) ≅ Z₂xZ₂⋊Sp(2, Z₂) (semidirect product)
by identification of 4 'sections', group homomorphisms s: Sp(2, Z₂) → C/U(1).
Implementations of group laws with and without 2-cocycle ω: Sp(2, Z₂) × Sp(2, Z₂) → Z₂xZ₂.
Assigning of words to elements.

Group labels:
Z₂xZ₂ -> 'P' where p ∈ P
C/U(1) -> 'C' where c ∈ C
Sp(2, Z₂) -> 'Q' where q ∈ Q
"""

import itertools

from sage.all import GF, Matrix  # type: ignore

from carry_table import build_carry_table

F2 = GF(2)

CARRY_TABLE = build_carry_table()


class CliffordElement:
    def __init__(self, matrix, name=None):
        self.tableau = matrix
        self.name = name

    def __repr__(self):
        return str(self.tableau)


Q_ELEMS = {
    "h": Matrix(F2, [[0, 1], [1, 0]]),
    "s": Matrix(F2, [[1, 0], [1, 1]]),
    "i": Matrix(F2, [[1, 0], [0, 1]]),
    "hs": Matrix(F2, [[1, 1], [1, 0]]),
    "sh": Matrix(F2, [[0, 1], [1, 1]]),
    "shs": Matrix(F2, [[1, 1], [0, 1]]),
}


def mat_key(mat):
    return tuple(map(tuple, mat.rows()))


REV_Q_ELEMS = {mat_key(mat): name for name, mat in Q_ELEMS.items()}

TRIVIAL_SECTION = {
    "i": CliffordElement(Matrix(F2, [[1, 0, 0], [0, 1, 0]]), "I"),
    "h": CliffordElement(Matrix(F2, [[0, 1, 0], [1, 0, 0]]), "H"),
    "s": CliffordElement(Matrix(F2, [[1, 0, 0], [1, 1, 0]]), "S"),
    "hs": CliffordElement(Matrix(F2, [[1, 1, 0], [1, 0, 0]]), "HS"),
    "sh": CliffordElement(Matrix(F2, [[0, 1, 0], [1, 1, 0]]), "SH"),
    "shs": CliffordElement(Matrix(F2, [[1, 1, 0], [0, 1, 0]]), "SHS"),
}


def group_action(q: Matrix, p: Matrix) -> Matrix:
    """
    Nontrivial action.
    """
    return q * p


def cocycle(q1: Matrix, q2: Matrix) -> Matrix:
    """
    2-cocycle ω: Q × Q → P, satisfying:
    s(q1)·s(q2) = ω(q1, q2) · s(q1q2)
    """
    q1_name = REV_Q_ELEMS[mat_key(q1)].upper()
    q2_name = REV_Q_ELEMS[mat_key(q2)].upper()
    if (q1_name, q2_name) in CARRY_TABLE:
        p = CARRY_TABLE[(q1_name, q2_name)]
        return Matrix(F2, [[int(p[0])], [int(p[1])]])
    else:
        raise KeyError(f"Carry table missing entry for ({q1_name}, {q2_name})")


def multiply_with_cocycle(c1: CliffordElement, c2: CliffordElement) -> CliffordElement:
    """
    Reference: Conrad Math 210B. Group cohomology and group extensions.
    (p1, q1)(p2, q2) = (p1 + q1·p2 + ω(q1,q2), q1q2)
    """
    g, v = c1.tableau[:, 0:2], c1.tableau[:, 2]
    h, w = c2.tableau[:, 0:2], c2.tableau[:, 2]

    gh = g * h
    new_v = v + g * w + cocycle(g, h)

    return CliffordElement(gh.augment(new_v))


def multiply_without_cocycle(
    c1: CliffordElement, c2: CliffordElement
) -> CliffordElement:
    """
    (p1, q1)(p2, q2) = (p1 + q1·p2, q1q2)
    """
    g, v = c1.tableau[:, 0:2], c1.tableau[:, 2]
    h, w = c2.tableau[:, 0:2], c2.tableau[:, 2]

    gh = g * h
    new_v = v + g * w

    return CliffordElement(gh.augment(new_v), c1.name + c2.name)


def two_cocycle_condition(r_name, s_name, t_name):
    """
    Reference: Conrad Math 210B. Group cohomology and group extensions (notation simplified).
    Derived from inserting associativity condition into group law.
    r·ω(s, t) + ω(rs, t) + ω(r, st) + ω(r, s) = 0
    where r, s, t ∈ P
    """
    r = Q_ELEMS[r_name.lower()]
    s = Q_ELEMS[s_name.lower()]
    t = Q_ELEMS[t_name.lower()]

    rs = r * s
    st = s * t

    term1 = group_action(r, cocycle(s, t))
    term2 = cocycle(rs, t)
    term3 = cocycle(r, st)
    term4 = cocycle(r, s)

    return term1 + term2 + term3 + term4 == Matrix(F2, [[0], [0]])


def two_cocycle_test():
    for r, s, t in itertools.product(["h", "s", "i"], repeat=3):
        if not two_cocycle_condition(r, s, t):
            print(f"Testing ({r}, {s}, {t})")
            print("FALSE")
            return False
    return True


def generate_group():
    """
    Generate the group from the 'trivial section' generators. Includes cocycle term in group law.
    """
    generators = [TRIVIAL_SECTION["i"], TRIVIAL_SECTION["h"], TRIVIAL_SECTION["s"]]
    all_tableaux = generators

    added_new = True
    while added_new:
        added_new = False
        new_elements = []

        for gen in generators:
            for tab in all_tableaux:
                try:
                    new_tab = multiply_with_cocycle(gen, tab)
                except KeyError:
                    continue

                if all(
                    new_tab.tableau != existing.tableau
                    for existing in all_tableaux + new_elements
                ):
                    new_elements.append(new_tab)
                    added_new = True

        all_tableaux.extend(new_elements)

    return all_tableaux


def generate_group_no_cocycle(section):
    """
    Generate the group from an arbitrary section's generators. No cocycle term in group law.
    """
    generators = [section["i"], section["h"], section["s"]]
    all_tableaux = generators.copy()

    added_new = True
    while added_new:
        added_new = False
        new_elements = []

        for gen in generators:
            for tab in all_tableaux:
                new_tab = multiply_without_cocycle(gen, tab)

                if all(
                    new_tab.tableau != existing.tableau
                    for existing in all_tableaux + new_elements
                ):
                    new_elements.append(new_tab)
                    added_new = True

        all_tableaux.extend(new_elements)

    return all_tableaux


def is_section_homomorphism(section):
    """
    Test whether a section s is a group homomorphism Q → C.
    s(q1)·s(q2) = s(q1·q2)
    """
    for q1, q2 in itertools.product(Q_ELEMS.keys(), repeat=2):
        try:
            lhs = multiply_with_cocycle(section[q1], section[q2])
        except KeyError:
            print(section)
        q1q2 = Q_ELEMS[q1] * Q_ELEMS[q2]
        q1q2_name = REV_Q_ELEMS[mat_key(q1q2)]
        rhs = section[q1q2_name]
        if lhs.tableau != rhs.tableau:
            return False
    return True


def all_p():
    return [Matrix(F2, [[a], [b]]) for a in [0, 1] for b in [0, 1]]


def build_tableau(q, p, name):
    return CliffordElement(q.augment(p), name)


def brute_force_sections():
    """
    Generate all 4^6 = 4096 possible sections and test them to find group homormphisms by brute force.
    """
    p_options = all_p()
    names = ["i", "h", "s", "hs", "sh", "shs"]
    successful_sections = []

    count = 0
    for p_combo in itertools.product(p_options, repeat=6):
        section = {}
        for i, name in enumerate(names):
            section[name] = build_tableau(Q_ELEMS[name], p_combo[i], name.upper())

        if is_section_homomorphism(section):
            print(f"Found homomorphic section! #{len(successful_sections) + 1}")
            successful_sections.append(section)

        count += 1
        if count % 500 == 0:
            print(f"Tried {count} sections...")

    print(f"Total homomorphic sections found: {len(successful_sections)}")
    return successful_sections


if __name__ == "__main__":
    print("Generating group")

    group = generate_group()
    assert len(group) == 24, "Error: Expected order-24 Clifford group"

    print("Group order: ", len(group))

    print("Checking 2-cocycle condition: ", two_cocycle_test())

    print("Searching for sections by brute force...")
    homomorphic_sections = brute_force_sections()

    if homomorphic_sections:
        for i, s in enumerate(homomorphic_sections):
            print(f"Homomorphic section {i + 1}")
            for key, mat in s.items():
                print(f"{key}:")
                print(mat)
                print()
    else:
        print("No homomorphic sections found.")

    print("Can I generate the group with these?")
    print()

    sections_orders = {}

    for i, section_h in enumerate(homomorphic_sections):
        group = generate_group_no_cocycle(section_h)

        sections_orders[i + 1] = len(group)

    print("Sections and their respective group orders")
    print(sections_orders)

    for element in group:
        print(element.name)
        print(element)
        print()

"""
Carry function is a 2-cocycle.

Found 4 alternate sections where multiplication to order-24 group works with no cocycle.
(Also a generator selection (I has [1, 1] phase vector) such that the cocycle is 'encoded'.)

Existence of sections means:
    - The cocycle is “trivial in cohomology".
    - The group extension is split (by 4 distinct sections): G≅Q⋊V, giving a semidirect product structure.
    - The cocycle is actually a *coboundary*.

When I define my section, I am picking my representatives.
When I then generate the full group from this section, I have a complete set of representatives of the group.
I can easily check if any two elements are equal (up to global phase).
We can multiply tableaux according to our group law.
We solve the *word problem* as we have a word for each of our elements that we can reduce any string down to.
"""
