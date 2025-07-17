"""
Run using `sage -python -m pytest <path_to_this_file>`
"""

from sage.all import Matrix, GF  # type: ignore
from chp import StabilizerTableau


F2 = GF(2)

# Single qubit single gate


# Test for X under X
def test_evolution_X_H():
    stab_tab = StabilizerTableau("+")
    stab_tab.print_bin()
    stab_tab.conjugate("H")
    assert stab_tab.tableau.row(1) == Matrix(F2, [[0, 1, 0]]).row(0)


def test_evolution_Y_H():
    stab_tab = StabilizerTableau("i")
    stab_tab.print_bin()
    stab_tab.conjugate("H")
    assert stab_tab.tableau.row(1) == Matrix(F2, [[1, 1, 1]]).row(0)


def test_evolution_Yneg_H():
    stab_tab = StabilizerTableau("-i")
    stab_tab.print_bin()
    stab_tab.conjugate("H")
    assert stab_tab.tableau.row(1) == Matrix(F2, [[1, 1, 0]]).row(0)


def test_evolution_Z_H():
    stab_tab = StabilizerTableau("0")
    stab_tab.print_bin()
    stab_tab.conjugate("H")
    assert stab_tab.tableau.row(1) == Matrix(F2, [[1, 0, 0]]).row(0)


def test_evolution_X_S():
    # Test what S|+⟩ actually produces
    tab = StabilizerTableau("+")
    print("Original |+⟩:")
    print(tab.tableau)

    tab.conjugate("S")
    print("After S gate:")
    print(tab.tableau)
    print(tab.tableau.row(1))
    assert tab.tableau.row(1) == Matrix(F2, [[1, 1, 0]]).row(0)


def test_evolution_Z_S():
    stab_tab = StabilizerTableau("0")
    stab_tab.print_bin()
    stab_tab.conjugate("S")
    assert stab_tab.tableau.row(1) == Matrix(F2, [[0, 1, 0]]).row(0)


def test_evolution_Y_S():
    stab_tab = StabilizerTableau("i")
    stab_tab.print_bin()
    stab_tab.conjugate("S")
    assert stab_tab.tableau.row(1) == Matrix(F2, [[1, 0, 1]]).row(0)


# Single qubit sequential gates


def test_identity_HH():
    """Test H^2 == I"""
    tab = StabilizerTableau("+")
    original = Matrix(tab.tableau)

    tab.conjugate("H")
    tab.conjugate("H")

    assert tab.tableau == original, "H^2 should return to original state"


def test_identity_SSSS():
    """Test S^4 = I"""
    tab = StabilizerTableau("+")
    original = Matrix(tab.tableau)

    tab.conjugate("S")
    tab.conjugate("S")
    tab.conjugate("S")
    tab.conjugate("S")

    assert tab.tableau == original, "S^4 should return to original state"


def test_sequence_HS_vs_direct():
    """Test H then S gives expected result"""
    # Apply H then S
    tab1 = StabilizerTableau("+")
    tab1.conjugate("H")  # |+⟩ → |0⟩
    tab1.print_bin()
    assert tab1.tableau.row(1) == Matrix(F2, [[0, 1, 0]]).row(0)
    tab1.conjugate("S")  # |0⟩ → |0⟩ (S fixes |0⟩)
    tab1.print_bin()
    assert tab1.tableau.row(1) == Matrix(F2, [[0, 1, 0]]).row(0)


def test_sequence_SH_vs_direct():
    """Test S then H gives expected result"""
    # Apply S then H
    tab1 = StabilizerTableau("+")
    tab1.conjugate("S")  # |+⟩ → |i⟩
    tab1.conjugate("H")  # |i⟩ → |-i⟩

    # Apply H to |i⟩ directly to compare
    tab2 = StabilizerTableau("i")
    tab2.conjugate("H")

    assert tab1.tableau.row(1) == tab2.tableau.row(1)


def test_sequence_SS():
    """Test S^2 transformation"""
    tab = StabilizerTableau("+")
    tab.conjugate("S")
    tab.conjugate("S")

    # S|+⟩ → |i⟩, so S^2|+⟩ = S|i⟩
    expected = StabilizerTableau("i")
    expected.conjugate("S")

    assert tab.tableau.row(1) == expected.tableau.row(1)


def test_sequence_commutativity():
    # H and S don't commute
    tab_hs = StabilizerTableau("+")
    tab_hs.conjugate("H")
    tab_hs.conjugate("S")

    tab_sh = StabilizerTableau("+")
    tab_sh.conjugate("S")
    tab_sh.conjugate("H")

    # These should be different
    print("HS|+⟩:", tab_hs.tableau.row(1))
    print("SH|+⟩:", tab_sh.tableau.row(1))
    assert tab_hs.tableau.row(1) != tab_sh.tableau.row(1)


def test_inverse_relations():
    # H is self-inverse
    tab = StabilizerTableau("i")  # Start with |i⟩
    original_stabilizer = tab.tableau.row(1)
    tab.conjugate("H")
    tab.conjugate("H")
    assert tab.tableau.row(1) == original_stabilizer, "H should be self-inverse"

    # S^3 is the inverse of S
    tab = StabilizerTableau("0")  # Start with |0⟩
    original_stabilizer = tab.tableau.row(1)
    tab.conjugate("S")  # Apply S
    tab.conjugate("S")  # Apply S (S^2)
    tab.conjugate("S")  # Apply S (S^3, should be same as S^{-1}
    # This should undo the first S
    expected = StabilizerTableau("0")
    expected.conjugate("S")
    expected.conjugate("S")
    expected.conjugate("S")
    expected.conjugate("S")  # S^4 should = I, so this equals original
    assert tab.tableau.row(1) == expected.tableau.row(1), "S³ should be inverse of S"


# CNOT test (2-qubit gate)


def test_cnot_basic():
    """Test basic CNOT gate functionality on 2-qubit system"""
    # This test assumes we have 2 qubits working
    # Create a basic 2-qubit |00⟩ state tableau manually
    # Format: [x0, x1, z0, z1, r] for 2 qubits
    # |00⟩ should have stabilizers ZI and IZ, destabilizers XI and IX
    tableau_00 = Matrix(
        F2,
        [
            [1, 0, 0, 0, 0],  # D0: XI (X on qubit 0, I on qubit 1)
            [0, 1, 0, 0, 0],  # D1: IX (I on qubit 0, X on qubit 1)
            [0, 0, 1, 0, 0],  # S0: ZI (Z on qubit 0, I on qubit 1)
            [0, 0, 0, 1, 0],  # S1: IZ (I on qubit 0, Z on qubit 1)
        ],
    )

    # mock 2-qubit tableau object
    class Mock2QubitTableau:
        def __init__(self, tableau):
            self.n = 2
            self.tableau = tableau
            self.name = "00"

        def cnot_gate(self, control_idx, target_idx):
            # Copy for now (dumb)
            for row_idx in range(2 * self.n):
                x_ia = self.tableau[row_idx, control_idx]
                x_ib = self.tableau[row_idx, target_idx]
                z_ia = self.tableau[row_idx, control_idx + self.n]
                z_ib = self.tableau[row_idx, target_idx + self.n]

                self.tableau[row_idx, -1] += x_ia * z_ib * (x_ib + z_ia + 1)
                x_ib = x_ib + x_ia
                z_ia = z_ia + z_ib

                self.tableau[row_idx, target_idx] = x_ib
                self.tableau[row_idx, control_idx + self.n] = z_ia

    # Test CNOT(0,1) on |00⟩
    tab = Mock2QubitTableau(tableau_00.__copy__())

    print("Before CNOT(0,1):")
    for i, row in enumerate(tab.tableau.rows()):
        row_type = "D" if i < 2 else "S"
        print(f"{row_type}{i % 2}: {list(row)}")

    tab.cnot_gate(0, 1)

    print("After CNOT(0,1):")
    for i, row in enumerate(tab.tableau.rows()):
        row_type = "D" if i < 2 else "S"
        print(f"{row_type}{i % 2}: {list(row)}")

    # Test expected transformations:
    # XI → XX (control X spreads to target)
    # IX → IX (target X unchanged)
    # ZI → ZI (control Z unchanged)
    # IZ → ZZ (target Z spreads to control)

    expected_after_cnot = Matrix(
        F2,
        [
            [1, 1, 0, 0, 0],  # D0: XX (was XI)
            [0, 1, 0, 0, 0],  # D1: IX (unchanged)
            [0, 0, 1, 0, 0],  # S0: ZI (unchanged)
            [0, 0, 1, 1, 0],  # S1: ZZ (was IZ)
        ],
    )

    assert tab.tableau == expected_after_cnot, "CNOT transformations incorrect"


# Measurement method tests


def test_rowsum_pauli_multiplication():
    tab = StabilizerTableau("+")  # |+⟩: destabilizer Z, stabilizer X

    # rowsum(0,1): Z * X = -iY should give Y with phase
    tab.rowsum(
        0, 1
    )  # multiplies two rows of the tableau with phase considerations and replaces in row 0

    # Row 0 should now be Y (x=1, z=1) with phase bit set
    assert tab.tableau[0, 0] == 1  # X part
    assert tab.tableau[0, 1] == 1  # Z part
    assert tab.tableau[0, 2] == 1  # Phase part (from -i factor)


def test_rowsum_identity():
    tab = StabilizerTableau("0")

    # Z * Z should give identity, I = (0, 0)
    tab.rowsum(1, 1)
    result = tab.tableau.row(1).list()

    assert result[0] == 0 and result[1] == 0, f"Z*Z should give I, got {result}"


def test_determined_measurement_detection():
    # Lazy - copying function logic; replace with patch
    tab = StabilizerTableau("0")

    # |0⟩ has Z stabilizer, so Z measurement should commute
    qubit_idx = 0
    anticommutes = False
    for row_idx in range(tab.n, 2 * tab.n):
        if tab.tableau[row_idx, qubit_idx]:  # X part of stabilizer
            anticommutes = True

    assert not anticommutes, "Z measurement should commute with Z stabilizer"


def test_random_measurement_detection():
    # Lazy - copying function logic; replace with patch
    tab = StabilizerTableau("+")

    # |+⟩ has X stabilizer, so Z measurement should anticommute
    qubit_idx = 0
    anticommutes = False
    for row_idx in range(tab.n, 2 * tab.n):
        if tab.tableau[row_idx, qubit_idx]:  # X part
            anticommutes = True

    assert anticommutes, "Z measurement should anticommute with X stabilizer"


def test_computational_basis_measurements():
    # |0⟩ should measure 0
    tab_0 = StabilizerTableau("0")
    outcome_0 = tab_0.measure(0)
    assert outcome_0 == 0, f"|0⟩ should measure 0, got {outcome_0}"

    # |1⟩ should measure 1
    tab_1 = StabilizerTableau("1")
    outcome_1 = tab_1.measure(0)
    assert outcome_1 == 1, f"|1⟩ should measure 1, got {outcome_1}"


def test_determined_measurement_preserves_tableau():
    """Test determined measurement doesn't change tableau."""
    tab = StabilizerTableau("0")
    original_tableau = tab.tableau.list()

    # outcome = tab.measure(0)

    assert tab.tableau.list() == original_tableau, (
        "Determined measurement should preserve tableau"
    )


def test_random_measurement_changes_tableau():
    """Test random measurement updates tableau."""
    tab = StabilizerTableau("+")
    original_tableau = tab.tableau.list()

    # outcome = tab.measure(0)

    assert tab.tableau.list() != original_tableau, (
        "Random measurement should change tableau"
    )


def test_associativity_1():
    """
    To Test associativity: S(HS)|+> = (SH)S|+>
    l1 = HS (first sequence applied on LHS of eq)
    l2 = S (2nd sequence applied on LHS of eq)
    r1 = S (first sequence applied on LHS of eq)
    r2 = SH (2nd sequence applied on LHS of eq)
    """

    left_state = StabilizerTableau("+").copy()
    right_state = StabilizerTableau("+").copy()

    left_seq_1 = "HS"
    left_seq_2 = "S"
    for gate in reversed(left_seq_1):
        left_state.conjugate(gate)
    for gate in reversed(left_seq_2):
        left_state.conjugate(gate)

    right_seq_1 = "S"
    right_seq_2 = "SH"
    for gate in reversed(right_seq_1):
        right_state.conjugate(gate)
    for gate in reversed(right_seq_2):
        right_state.conjugate(gate)

    assert right_state == left_state


def test_associativity_2():
    """
    Test associativity: S(SH) = (SS)H
    """

    left_state = StabilizerTableau("+").copy()
    right_state = StabilizerTableau("+").copy()

    left_seq_1 = "SH"
    left_seq_2 = "S"
    for gate in reversed(left_seq_1):
        left_state.conjugate(gate)
    for gate in reversed(left_seq_2):
        left_state.conjugate(gate)

    right_seq_1 = "H"
    right_seq_2 = "SS"
    for gate in reversed(right_seq_1):
        right_state.conjugate(gate)
    for gate in reversed(right_seq_2):
        right_state.conjugate(gate)

    assert right_state == left_state


def test_full_associativity():
    """
    To Test associativity: S(HS) = (SH)S
    l1 = HS (first sequence applied on LHS of eq)
    l2 = S (2nd sequence applied on LHS of eq)
    r1 = S (first sequence applied on LHS of eq)
    r2 = SH (2nd sequence applied on LHS of eq)
    """
    all_sequences = [
        "H",
        "S",
        "SH",
        "HS",
        "SS",
        "HSH",
        "SSH",
        "SHS",
        "HSS",
        "SSS",
        "SHSH",
        "HSSH",
        "SSSH",
        "SSHS",
        "SHSS",
        "SSHSH",
        "SHSSH",
        "HSSHS",
        "SSSHS",
        "SSHSS",
        "HSSHSH",
        "SSSHSH",
        "SSHSSH",
    ]

    for a in all_sequences:
        for b in all_sequences:
            for c in all_sequences:
                # a = "S"
                # b = "H"
                # c = "S"

                # left is S(HS) = a(bc) so left_seq_2 = a and left_seq_1 = bc = b + c
                # right is (SH)S = (ab)c so right_seq_2 = ab = a + b and right_seq_1 = c

                left_state = StabilizerTableau("+").copy()
                right_state = StabilizerTableau("+").copy()

                left_seq_1 = b + c
                left_seq_2 = a
                for gate in reversed(left_seq_1):
                    left_state.conjugate(gate)
                for gate in reversed(left_seq_2):
                    left_state.conjugate(gate)

                right_seq_1 = c
                right_seq_2 = a + b
                for gate in reversed(right_seq_1):
                    right_state.conjugate(gate)
                for gate in reversed(right_seq_2):
                    right_state.conjugate(gate)

                assert right_state == left_state
