"""
Builds the carry table for the 2-cocycle ω: Q × Q → P
used in the split extension C ≅ P ⋊ Q (Clifford group on 1 qubit).
Each carry is a pair of bits representing a Pauli correction.
"""

import numpy as np

I = np.array([[1, 0], [0, 1]], dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

H = (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)
S = np.array([[1, 0], [0, 1j]], dtype=complex)
SH = S @ H
HS = H @ S
SHS = S @ H @ S

# This is actually a (non-homomorphic) 'section'. The carry table depends on this choice
UNITARIES = {
    "I": I,
    "H": H,
    "S": S,
    "SH": SH,
    "HS": HS,
    "SHS": SHS,
}

PAULIS = {"I": I, "X": X, "Y": Y, "Z": Z}
PAULI_TO_CARRY = {"I": [0, 0], "X": [1, 0], "Z": [0, 1], "Y": [1, 1]}


def remove_global_phase(U):
    """return U with global phase stripped (det(U)=1)."""
    phase = np.angle(np.linalg.det(U)) / U.shape[0]
    return U * np.exp(-1j * phase)


def same_up_to_phase(A, B, atol=1e-10):
    """True if A ≈ e^{iφ} B."""
    return np.allclose(remove_global_phase(A), remove_global_phase(B), atol=atol)


def identify_pauli(U):
    """Which Pauli (incl I) is U up to phase?"""
    for name, P in PAULIS.items():
        if same_up_to_phase(U, P):
            return name
    return None


PHASES = [1, -1, 1j, -1j]

# This is the normal form kinda thing


def canonical_decomposition(prod):
    """
    Decompose a unitary as:
        prod ≈ phase · P · Q
    with:
        - phase ∈ {±1, ±i}
        - P ∈ Pauli group
        - Q ∈ canonical Clifford reps (Q_GENERATORS)
    Returns:
        (pauli_label, Q_label)
    """
    for phase in PHASES:
        M = remove_global_phase(prod / phase)
        for p_name, P in PAULIS.items():
            for q_name, Q in UNITARIES.items():
                if np.allclose(M, remove_global_phase(P @ Q), atol=1e-10):
                    return p_name, q_name
    raise ValueError("Input matrix not decomposable into phase * Pauli * Clifford")


def compute_carry(name1, name2):
    U1 = UNITARIES[name1]
    U2 = UNITARIES[name2]
    raw_prod = U1 @ U2

    pauli_name, canon_name = canonical_decomposition(raw_prod)
    carry_bits = PAULI_TO_CARRY[pauli_name]
    return carry_bits, canon_name


def build_carry_table():
    table = {}
    for u1 in UNITARIES:
        for u2 in UNITARIES:
            carry, _ = compute_carry(u1, u2)
            table[(u1, u2)] = carry
    return table


if __name__ == "__main__":
    carry_table = build_carry_table()
    for k, v in carry_table.items():
        print(f"{k}: {v}")