"""
Run using `sage -python -m pytest <path_to_this_file>`
"""
import pytest
import random
from itertools import combinations, product
from sage.all import Matrix, GF  # type: ignore
from clifford import Operator, CliffordTableau, StabilizerTableau, CLIFFORD_LOOKUP, STABILIZER_LOOKUP, generate_Clifford_group


F2 = GF(2)

# Test for Z under X
def test_evolution_Z_X():
    stab_tab = StabilizerTableau("0")
    stab_tab.print_bin()
    x_tab = CliffordTableau("X")
    x_tab.print_bin()
    stab_tab.conjugate(x_tab)
    assert stab_tab.tableau == Matrix(F2, [[0, 1, 1, 0]])

def test_evolution_Z_Z():
    stab_tab = StabilizerTableau("0")
    stab_tab.print_bin()
    x_tab = CliffordTableau("Z")
    x_tab.print_bin()
    stab_tab.conjugate(x_tab)
    assert stab_tab.tableau == Matrix(F2, [[0, 1, 0, 0]])

def test_evolution_X_X():
    stab_tab = StabilizerTableau("+")
    stab_tab.print_bin()
    x_tab = CliffordTableau("X")
    x_tab.print_bin()
    stab_tab.conjugate(x_tab)
    assert stab_tab.tableau == Matrix(F2, [[1, 0, 0, 0]])

def test_evolution_Xneg_X():
    stab_tab = StabilizerTableau("-")
    stab_tab.print_bin()
    x_tab = CliffordTableau("X")
    x_tab.print_bin()
    stab_tab.conjugate(x_tab)
    assert stab_tab.tableau == Matrix(F2, [[1, 0, 1, 0]])

def test_evolution_X_H():
    stab_tab = StabilizerTableau("+")
    stab_tab.print_bin()
    x_tab = CliffordTableau("H")
    x_tab.print_bin()
    stab_tab.conjugate(x_tab)
    assert stab_tab.tableau == Matrix(F2, [[0, 1, 0, 0]])

def test_evolution_Y_H():
    stab_tab = StabilizerTableau("i")
    stab_tab.print_bin()
    x_tab = CliffordTableau("H")
    x_tab.print_bin()
    stab_tab.conjugate(x_tab)
    assert stab_tab.tableau == Matrix(F2, [[1, 1, 1, 1]])

def test_evolution_Yneg_H():
    stab_tab = StabilizerTableau("-i")
    stab_tab.print_bin()
    x_tab = CliffordTableau("H")
    x_tab.print_bin()
    stab_tab.conjugate(x_tab)
    assert stab_tab.tableau == Matrix(F2, [[1, 1, 0, 1]])

def test_evolution_Z_H():
    stab_tab = StabilizerTableau("0")
    stab_tab.print_bin()
    x_tab = CliffordTableau("H")
    x_tab.print_bin()
    stab_tab.conjugate(x_tab)
    assert stab_tab.tableau == Matrix(F2, [[1, 0, 0, 0]])

def test_evolution_X_S():
    stab_tab = StabilizerTableau("+")
    stab_tab.print_bin()
    x_tab = CliffordTableau("S")
    x_tab.print_bin()
    stab_tab.conjugate(x_tab)
    assert stab_tab.tableau == Matrix(F2, [[1, 1, 0, 1]])

def test_evolution_Z_S():
    stab_tab = StabilizerTableau("0")
    stab_tab.print_bin()
    x_tab = CliffordTableau("S")
    x_tab.print_bin()
    stab_tab.conjugate(x_tab)
    assert stab_tab.tableau == Matrix(F2, [[0, 1, 0, 0]])

def test_evolution_i_S():
    stab_tab = StabilizerTableau("i")
    stab_tab.print_bin()
    x_tab = CliffordTableau("S")
    x_tab.print_bin()
    stab_tab.conjugate(x_tab)
    assert stab_tab.tableau == Matrix(F2, [[1, 0, 1, 0]])

def test_Clifford_compose_HS():
    h = CliffordTableau("H")
    s = CliffordTableau("S")
    hs = h.compose(s)
    # Try one by one
    stab_tab_1 = StabilizerTableau("+")
    stab_tab_1.conjugate(s)
    stab_tab_1.conjugate(h)
    # Try together
    stab_tab = StabilizerTableau("+")
    stab_tab.conjugate(hs)
    # Compare: is composition function same as sequential operations?
    assert stab_tab == stab_tab_1

def test_Clifford_compose_XZ():
    # Test by action one by one then together on a stabilizer state
    # Depends on conjugation working correctly
    x = CliffordTableau("X")
    z = CliffordTableau("Z")
    xz = x.compose(z)
    # Try one by one
    stab_tab_one_by_one = StabilizerTableau("+")
    stab_tab_one_by_one.conjugate(z)
    stab_tab_one_by_one.conjugate(x)
    # Try together
    stab_tab_together = StabilizerTableau("+")
    stab_tab_together.conjugate(xz)
    # Compare: is composition function same as sequential operations?
    assert stab_tab_together == stab_tab_one_by_one

def test_Clifford_compose_SS():
    # Test by equality with another Clifford operator
    sone = CliffordTableau("S")
    print("Tableau for first S")
    print(sone.tableau)
    print("Tableau for 2nd S")
    stwo = CliffordTableau("S")
    print(stwo.tableau)
    print("Tableau for Z")
    z = CliffordTableau("Z")
    print(z.tableau)
    ss = sone.compose(stwo)
    print("Tableau for SS")
    print(ss.tableau)
    assert z == ss

def test_Clifford_compose_XX():
    id = CliffordTableau("I")
    xone = CliffordTableau("X")
    xtwo = CliffordTableau("X")
    xx = xone.compose(xtwo)
    print("Tableau for XX")
    print(xx.tableau)
    assert id == xx


# Group theory tests

@pytest.fixture
def c1():
    return generate_Clifford_group()

class TestCliffordGroupStructure:
    
    def test_identity_element(self, c1):
        """
        Test that identity element exists and behaves correctly
        """
        identity = CliffordTableau("I")

        assert identity in c1
        
        for gen_name in CLIFFORD_LOOKUP.keys():
            gen = CliffordTableau(gen_name)
            
            left_mult = identity.compose(gen)
            assert left_mult == gen, f"Left identity failed for {gen_name}"
            
            right_mult = gen.compose(identity)
            assert right_mult == gen, f"Right identity failed for {gen_name}"
    
    def test_associativity(self, c1):
        """
        Test associativity: (a*b)*c = a*(b*c)
        This is only testing one set of 3, need to increase spread
        Plus I think this is sometimes failing so there is still a composition bug...
        """
        test_elements = random.sample(list(c1), 3)
        print(test_elements)
        
        a, b, c = test_elements[0], test_elements[1], test_elements[2]
        # (a*b)*c
        ab = a.compose(b)
        abc_left = ab.compose(c)
        print("(a*b)*c")
        abc_left.print_bin()
        
        # a*(b*c)  
        bc = b.compose(c)
        abc_right = a.compose(bc)
        print("a*(b*c)")
        abc_right.print_bin()
        
        assert abc_left == abc_right, f"Associativity failed for ({a}*{b})*{c} vs {a}*({b}*{c})"
    
    def test_closure_under_composition(self, c1):
        """
        Test that the group is closed under composition
        """
        test_elements = random.sample(list(c1), 3)
        print(test_elements)
        
        for a in test_elements:
            for b in test_elements:
                product = a.compose(b)
                assert product is not None, f"Composition {a} * {b} failed"
                assert product in c1

    def test_paulis_normalized(self, c1):
        """
        Test that Pauli matrices are normalized by some random elements of the group
        To normalize I need a check for containment in P_{1}
        """
        paulis = [StabilizerTableau(st_name) for st_name in STABILIZER_LOOKUP.keys()]
        test_elements = random.sample(list(c1), 3)
        print(test_elements)
        
        for p1 in test_elements:
            for p2 in test_elements:
                product = p1.compose(p2)

                for pauli in paulis:
                    result = pauli.conjugate(product)

                    # Assert that each result is in the pauli group - need to generate group 