"""
Generate a carry table to use as the 2-cocycle
"""

import numpy as np

I = np.array([[1, 0], [0, 1]], dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

H = (1 / np.sqrt(2)) * np.array([[1, 1],
                                 [1, -1]], dtype=complex)
S = np.array([[1, 0],
              [0, 1j]], dtype=complex)
SH  = S @ H
HS  = H @ S
SHS = S @ H @ S

unitaries = {
    'I': I,
    'H': H,
    'S': S,
    'SH': SH,
    'HS': HS,
    'SHS': SHS,
}

PAULIS = {'I': I, 'X': X, 'Y': Y, 'Z': Z}
PAULI_TO_CARRY = {'I': [0, 0], 'X': [1, 0], 'Z': [0, 1], 'Y': [1, 1]}

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
    Return (pauli_name, canonical_name) such that
        prod ≈ phase * Pauli * canonical_unitary
    with phase ∈ {±1, ±i}.
    """
    for phase in PHASES:
        M = remove_global_phase(prod / phase)   # divide by phase, then strip any leftover
        for pname, P in PAULIS.items():
            for gname, U in unitaries.items():
                if np.allclose(M, remove_global_phase(P @ U), atol=1e-10):
                    return pname, gname
    raise ValueError("Matrix does not match any phase*Pauli*canonical element.")

def compute_carry(name1, name2):
    U1 = unitaries[name1]
    U2 = unitaries[name2]
    raw_prod = U1 @ U2

    pauli_name, canon_name = canonical_decomposition(raw_prod)
    carry_bits = PAULI_TO_CARRY[pauli_name]
    return carry_bits, canon_name

def build_carry_table():
    table = {}
    for g1 in unitaries:
        for g2 in unitaries:
            carry, _ = compute_carry(g1, g2)
            table[(g1, g2)] = carry
    return table


if __name__ == "__main__":
    carry_table = build_carry_table()
    for k, v in carry_table.items():
        print(f"{k}: {v}")