"""
If I change the mapping in SECTION such that I has [1, 1] phase vector,
and generate from this I with H, S, I can generate all 24 without the cocycle.
HOWEVER, this section is not a homomorphism (so it's not a section).
Found 4 sections that are group homomorphisms.
"""

from sage.all import Matrix, GF # type: ignore
from carry_table import build_carry_table
import itertools
from pprint import pprint

F2 = GF(2)

# CARRY_TABLE = CARRY_TABLE_SIMON
CARRY_TABLE = build_carry_table()

class Tableau:
    def __init__(self, matrix, name=None):
        self.data = matrix
        self.name = name

    def __repr__(self):
        return str(self.data)
    
# Define our 6 symplectic matrices each with a label in {h, s, i} (Simon's notation!!)
# Using lowercase to distinguish from operators
BS_MATRICES = {
    'h': Matrix(F2, [[0, 1], [1, 0]]),
    's': Matrix(F2, [[1, 0], [1, 1]]),
    'i': Matrix(F2, [[1, 0], [0, 1]]),
    'hs': Matrix(F2, [[1, 1], [1, 0]]),
    'sh': Matrix(F2, [[0, 1], [1, 1]]),
    'shs': Matrix(F2, [[1, 1], [0, 1]]),
    }

def mat_key(mat):
    return tuple(map(tuple, mat.rows()))

REV_BS_MATRICES = {mat_key(mat): name for name, mat in BS_MATRICES.items()}

# Choice of mapping for BS_matrices (6) by their BS_matrix name back to full tableaux
# This is the 'trivial section'
SECTION = {
    'i': Tableau(Matrix(F2, [[1, 0, 0], [0, 1, 0]]), 'I'),
    'h': Tableau(Matrix(F2, [[0, 1, 0], [1, 0, 0]]), 'H'),
    's': Tableau(Matrix(F2, [[1, 0, 0], [1, 1, 0]]), 'S'),
    'hs': Tableau(Matrix(F2, [[1, 1, 0], [1, 0, 0]]), 'HS'),
    'sh': Tableau(Matrix(F2, [[0, 1, 0], [1, 1, 0]]), 'SH'),
    'shs': Tableau(Matrix(F2, [[1, 1, 0], [0, 1, 0]]), 'SHS')
}

def group_action(g_matrix, phase_vector):
    """
    Nontrivial
    """
    return g_matrix * phase_vector

def carry(bs_1, bs_2):
    name1 = REV_BS_MATRICES[mat_key(bs_1)].upper()
    name2 = REV_BS_MATRICES[mat_key(bs_2)].upper()
    if (name1, name2) in CARRY_TABLE:
        phase_vector = CARRY_TABLE[(name1, name2)]
        return Matrix(F2, [[int(phase_vector[0])], [int(phase_vector[1])]])
    else:
        raise KeyError(f"Carry table missing entry for ({name1}, {name2})")
    
def multiply_tableaux(tab_1, tab_2):
    out = Matrix(F2, [[0, 0, 0], [0, 0, 0]])

    bs_1, m_1 = tab_1.data[:, 0:2], tab_1.data[:, 2]
    bs_2, m_2 = tab_2.data[:, 0:2], tab_2.data[:, 2]

    out[:, 0:2] = bs_1 * bs_2
    out[:, 2] = m_1 + group_action(bs_1, m_2) + carry(bs_1, bs_2) 

    result = Tableau(out)
    return result

def multiply_tableaux_no_carry(tab_1, tab_2):
    out = Matrix(F2, [[0, 0, 0], [0, 0, 0]])

    bs_1, m_1 = tab_1.data[:, 0:2], tab_1.data[:, 2]
    bs_2, m_2 = tab_2.data[:, 0:2], tab_2.data[:, 2]

    out[:, 0:2] = bs_1 * bs_2
    out[:, 2] = m_1 + group_action(bs_1, m_2)  # <- NO carry term

    return Tableau(out, tab_1.name+tab_2.name)

def two_cocycle_condition(r_name, s_name, t_name):
    """ r·c(s, t) + c(rs, t) + c(r, st) + c(r, s) = 0 """
    r = BS_MATRICES[r_name.lower()]
    s = BS_MATRICES[s_name.lower()]
    t = BS_MATRICES[t_name.lower()]

    term1 = group_action(r, carry(s, t))
    
    rs_mat  = r * s
    rs_name = REV_BS_MATRICES[mat_key(rs_mat)]

    term2 = carry(BS_MATRICES[rs_name], t)
    
    st_mat  = s * t
    st_name = REV_BS_MATRICES[mat_key(st_mat)]
    term3 = carry(r, BS_MATRICES[st_name])
    
    term4 = carry(r, s)
    
    result = term1 + term2 + term3 + term4
    
    zero_vector = Matrix(F2, [[0], [0]])
    satisfied = (result == zero_vector)

    return satisfied

def two_cocycle_test():
    for r, s, t in itertools.product(['h', 's', 'i'], repeat=3):
        if not two_cocycle_condition(r, s, t):
            print(f"Testing ({r}, {s}, {t})")
            print("FALSE")
            return False
    return True

def generate_group():
    generators = [SECTION['i'], SECTION['h'], SECTION['s']]
    all_tableaux = generators

    added_new = True
    while added_new:
        added_new = False
        new_elements = []

        for gen in generators:
            for tab in all_tableaux:
                try:
                    new_tab = multiply_tableaux(gen, tab)
                except KeyError:
                    continue

                if all(new_tab.data != existing.data for existing in all_tableaux + new_elements):
                    new_elements.append(new_tab)
                    added_new = True
        
        all_tableaux.extend(new_elements)

    return all_tableaux

def generate_group_no_carry(section):
    generators = [section['i'], section['h'], section['s']]
    all_tableaux = generators.copy()

    added_new = True
    while added_new:
        added_new = False
        new_elements = []

        for gen in generators:
            for tab in all_tableaux:
                new_tab = multiply_tableaux_no_carry(gen, tab)

                if all(new_tab.data != existing.data for existing in all_tableaux + new_elements):
                    new_elements.append(new_tab)
                    added_new = True
        
        all_tableaux.extend(new_elements)

    return all_tableaux

def label_section(group, section):
    group_dict = {}

    for element in group:
        group_dict[element] = None
        for key, value in section.items():
            if element.data == value.data:
                group_dict[element] = key

    return group_dict

def is_section_homomorphism(section):
    for g1, g2 in itertools.product(BS_MATRICES.keys(), repeat=2):
        try:
            lhs = multiply_tableaux(section[g1], section[g2])
        except KeyError:
            print(section)
        g1g2 = BS_MATRICES[g1] * BS_MATRICES[g2]
        g1g2_name = REV_BS_MATRICES[mat_key(g1g2)]
        rhs = section[g1g2_name]
        if lhs.data != rhs.data:
            return False
    return True

def all_phase_vectors():
    return [Matrix(F2, [[a], [b]]) for a in [0, 1] for b in [0, 1]]

def build_tableau(symplectic_matrix, phase_vector, name):
    return Tableau(symplectic_matrix.augment(phase_vector), name)

def brute_force_sections():
    phase_options = all_phase_vectors()
    names = ['i', 'h', 's', 'hs', 'sh', 'shs']
    successful_sections = []

    count = 0
    for phase_vector_combo in itertools.product(phase_options, repeat=6):
        section = {}
        for i, name in enumerate(names):
            section[name] = build_tableau(BS_MATRICES[name], phase_vector_combo[i], name.upper())

        if is_section_homomorphism(section):
            print(f"Found homomorphic section! #{len(successful_sections)+1}")
            successful_sections.append(section)

        count += 1
        if count % 500 == 0:
            print(f"Tried {count} sections...")

    print(f"Total homomorphic sections found: {len(successful_sections)}")
    return successful_sections


if __name__ == "__main__":

    print("Generating group")

    group = generate_group()

    print("Group order: ", len(group))

    pprint(label_section(group, SECTION))

    print("Checking two cocycle condition: ", two_cocycle_test())

    print("Brute-forcing sections...")
    homomorphic_sections = brute_force_sections()

    if homomorphic_sections:
        for i, s in enumerate(homomorphic_sections):
            print(f"Homomorphic section {i+1}")
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

        group = generate_group_no_carry(section_h)

        sections_orders[i+1] = len(group)

    print(sections_orders)

    for element in group:
        print(element.name)
        print(element)
        print()

"""
Carry function is a 2-cocycle.

Found 4 alternate sections where multiplication works with no carry.

This means the cocycle is “trivial in cohomology” i.e. it's not a real obstruction.

The group extension is split: G≅Q⋉V, a semidirect product.

The cocycle is actually a *coboundary*.

When I define my section, I am picking my representatives.
When I then generate the full group from this section, I have a complete set of representatives of the group.
I can easily check if any two elements are equal (up to global phase).
We can multiply tableaux according to our group law.
We solve the *word problem* when we have a word for each of these; I'll do that below vvvv
"""

