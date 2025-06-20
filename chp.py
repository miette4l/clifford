"""
An implementation of Aaronson and Gottesman's CHP program:
https://arxiv.org/pdf/quant-ph/0406196
"""

from sage.all import Matrix, GF, UniversalCyclotomicField, copy  # type: ignore
import numpy as np
import random


F2 = GF(2)
roots = UniversalCyclotomicField().gen(8)


# Stabilizer tableau for the 1-qubit stabilizer states in [ X | Z | r ]
# r is the phase bit: ri is 1 if Ri has negative phase and 0 if ri has positive phase
# First row is the destabilizer; 2nd row is the stabilizer
# N.b. multiple tableaux can represent a single state (different destabilizers)

STABILIZER_LOOKUP = {
    "0": Matrix(
        F2,
        [
            [1, 0, 0],  # Destabilizer: X
            [0, 1, 0],
        ],
    ),  # Stabilizer: Z
    "1": Matrix(F2, [[1, 0, 0], [0, 1, 1]]),  # -Z
    "+": Matrix(
        F2,
        [
            [0, 1, 0],  # Z
            [1, 0, 0],
        ],
    ),  # X
    "-": Matrix(F2, [[0, 1, 0], [1, 0, 1]]),  # -X
    "i": Matrix(F2, [[0, 1, 0], [1, 1, 0]]),  # Y
    "-i": Matrix(F2, [[0, 1, 0], [1, 1, 1]]),  # -Y
    "I": Matrix(F2, [[0, 0, 0], [0, 0, 0]]),
}


class StabilizerTableau:
    def __init__(self, stabilizer_state: str):
        self.name = stabilizer_state
        self.n = 1  # only 1 for now
        try:
            self.tableau = STABILIZER_LOOKUP[stabilizer_state].__copy__()
        except KeyError:
            raise ValueError(f"Unknown stabilizer state '{stabilizer_state}'")

    def __hash__(self):
        """
        Hash based only on stabilizer generators (rows n to 2n-1)
        This ensures that tableaux representing the same quantum state
        have the same hash, even if they have different destabilizers.
        """
        # Extract only the stabilizer rows (second half of tableau)
        # stabilizer_rows = []
        # for row_idx in range(self.n, 2 * self.n):
        #     stabilizer_rows.extend(self.tableau.row(row_idx).list())

        # return hash(tuple(stabilizer_rows))

        # Full tableaux hash:

        rows = []
        for row_idx in range(2 * self.n):
            rows.extend(self.tableau.row(row_idx).list())

        return hash(tuple(rows))

        # Interestingly this returns 24 states in the orbit - showing degeneracy 4?

    def __eq__(self, other):
        # This isn't very precise...
        return self.tableau == other.tableau

    def copy(self):
        new_tableau = StabilizerTableau("+")
        new_tableau.tableau = self.tableau.__copy__()
        new_tableau.name = self.name if self.name else "None"
        return new_tableau

    def print_bin(self):
        print(self.tableau)

    def print_string(self):
        if self.tableau is None:
            print("No tableau initialized")
            return

        for i in range(0, self.n):
            print(self.row_to_pauli_string(i), "Destabilizer")
        print(f"{'-' * 20}")
        for i in range(self.n, (2 * self.n)):
            print(self.row_to_pauli_string(i), "Stabilizer")

    def row_to_pauli_string(self, row_idx: int) -> str:
        """
        Tableau row [ X | Z | r ] -> Pauli string like '+XZY
        """
        row = self.tableau[row_idx, :]
        phase = "-" if row[0, -1] else "+"
        pauli_str = ""

        for qubit_idx in range(self.n):
            x = int(row[0, qubit_idx])  # X part: columns 0 to n-1
            z = int(row[0, qubit_idx + self.n])  # Z part: columns n to 2n-1

            if x == 0 and z == 0:
                pauli_str += "I"
            elif x == 1 and z == 0:
                pauli_str += "X"
            elif x == 0 and z == 1:
                pauli_str += "Z"
            else:  # x == 1 and z == 1
                pauli_str += "Y"

        return f"{phase}{pauli_str}"

    def pauli_string_to_row(self, pauli_str):
        """
        Pauli string like '+XZY' -> tableau row [ X | Z | r ]
        Returns a SageMath Matrix row
        """
        if len(pauli_str) < 2:
            raise ValueError(
                "Pauli string must have at least phase and one Pauli operator"
            )

        # Extract phase
        phase_char = pauli_str[0]
        if phase_char == "+":
            phase_bit = 0
        elif phase_char == "-":
            phase_bit = 1
        else:
            raise ValueError("Pauli string must start with '+' or '-'")

        # Extract n-qubit Pauli operator
        pauli_op = pauli_str[1:]
        if not len(pauli_op) == self.n:
            raise ValueError("Wrong number of Pauli operators")

        # Initialize row: [x0, x1, ..., xn-1, z0, z1, ..., zn-1, r]
        row_data = [0] * (2 * self.n + 1)

        for i, pauli_char in enumerate(pauli_op):
            if pauli_char == "I":
                x_bit, z_bit = 0, 0
            elif pauli_char == "X":
                x_bit, z_bit = 1, 0
            elif pauli_char == "Z":
                x_bit, z_bit = 0, 1
            elif pauli_char == "Y":
                x_bit, z_bit = 1, 1
            else:
                raise ValueError(f"Unknown Pauli operator: {pauli_char}")

            row_data[i] = x_bit
            row_data[i + self.n] = z_bit

        row_data[-1] = phase_bit

        row = Matrix(F2, [row_data])

        return row

    def conjugate(self, clifford: str, **kwargs):
        """
        Apply h_gate, p_gate, etc. via individual functions
        """
        if len(kwargs) == 0:
            qubit_idx = 0
        elif len(kwargs) == 1:
            qubit_idx = kwargs["qubit_idx"]
        elif len(kwargs) == 2:
            control_idx, target_idx = kwargs["control_idx"], kwargs["target_idx"]

        if clifford == "H":
            self.h_gate(qubit_idx)
        elif clifford == "S":
            self.s_gate(qubit_idx)
        elif clifford == "CNOT":
            self.cnot_gate(control_idx, target_idx)
        else:
            raise ValueError(f"Unknown Clifford gate: {clifford}")

        self.name = self.tableau_to_state_name()

        return self.tableau

    def h_gate(self, qubit_idx: int = 0):
        """
        Apply H gate via conjugation HCH of Pauli-operator stabilizers

        For all i ∈ {1, . . ., 2n}, set r_i = r_i ⊕ x_{ia}*z_{ia} and swap x_{ia} with z_{ia}
        (only n=1 for now so qubit_idx always 0)

        Introduces a global phase of +1 or -1 depending on the input state.
        """
        for row_idx in range(2 * self.n):
            x = self.tableau[row_idx, qubit_idx]
            z = self.tableau[row_idx, qubit_idx + self.n]

            self.tableau[row_idx, -1] += x * z

        self.tableau.swap_columns(qubit_idx, self.n + qubit_idx)

    def s_gate(self, qubit_idx: int):
        """
        For all i ∈ {1, . . ., 2n}, set r_i = r_i ⊕ x_{ia}*z_{ia}
        and then set z_{ia} = z_{ia} ⊕ x_{ia}

        Introduces a global phase of i or -i depending on the input state.
        But this isn't tracked... huh?
        """
        for row_idx in range(2 * self.n):
            x = self.tableau[row_idx, qubit_idx]
            z = self.tableau[row_idx, qubit_idx + self.n]

            # Update phase
            self.tableau[row_idx, -1] += x * z
            # I don't understand where i goes here...

            # Update Z column: z := z ⊕ x (do this row by row)
            self.tableau[row_idx, qubit_idx + self.n] += x

    def cnot_gate(self, control_idx: int, target_idx: int):
        """
        Apply CNOT with control and target

        For all i ∈ {1, . . ., 2n},
        set r_i = r_i ⊕ x_{ia}*z_{ib}*(x_{ib} ⊕ z_{ia} ⊕ 1),
        x_{ib} = x_{ib} ⊕ x_{ia},
        z_{ia} = z_{ia} ⊕ z_{ib}
        """
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

    def generate_tableau(self, n_qubits: int, state: str):
        """
        Generate a multi-qubit tableau (not implemented yet)
        I want this to dynamically generate not just check the lookup
        Plus move this outside the class then have it inject

        Should work like:
        Given state, do orbit trick to find stabilizers AND destabilizers as Pauli strings
        convert Pauli strings to rows
        construct tableau from rows of both types in correct order
        [ X | Z | r ]
        create StabilizerTableau object from this tableau
        """
        if n_qubits != 1:
            raise NotImplementedError("Only 1-qubit states supported currently.")
        self.tableau = copy(STABILIZER_LOOKUP[state])
        self.n = n_qubits
        return self.tableau

    def orbit_trick(self, n_qubits, state: str):
        """
        NOT DONE
        Should return a list of stabilizers and destabilizers for a known state as Pauli strings
        Should be generative not a lookup
        """
        if n_qubits > 1:
            raise NotImplementedError("n > 1 not supported")

        if state == "0":
            return ["Z"]
        elif state == "+":
            return ["X"]
        elif state == "i":
            return ["Y"]
        else:
            return ["Unknown"]

    def tableau_to_state_name(self):
        """
        NOT DONE
        I was rather this gave some kind of canonical form/unique repr
        Given a tableau, find which stabilizer state it corresponds to
        """
        for state_name, ref_tableau in STABILIZER_LOOKUP.items():
            if self.tableau.list() == ref_tableau.list():
                return state_name
        return "Unknown"

    def measure(self, qubit_idx):
        n = self.n

        # Check if Z on qubit_idx anticommutes with any stabilizer
        for row_idx in range(n, 2 * n):
            x_ia = self.tableau[row_idx, qubit_idx]
            if x_ia:
                return self.random_measurement(row_idx)  # Pass the row index as p
        return self.determined_measurement(qubit_idx)

    def random_measurement(self, p):
        """
        Coin toss first to actually simulate a random state, then perform:
            First call rowsum for all i ∈ {1, . . 2n},
            Set entire (p-n)th row equal to pth
            Set pth row to 0, after saving r
            Return r as measurement outcome
        """
        # Generate random outcome (coin flip)
        heads = random.randint(0, 1)
        n = self.n

        if heads:
            self.tableau[p, 2 * n] += 1

        # For all rows i in [0 .. 2n - 1], update row i by rowsum(i, p)
        for i in range(2 * n):
            if i != p:
                self.rowsum(i, p)

        # Set entire (p-n)th row equal to pth
        self.tableau[p - n, :] = self.tableau[p, :].__copy__()

        # Store measurement outcome (phase bit of p before zeroing)
        outcome = self.tableau[p, 2 * n]

        # Set pth row to 0
        self.tableau[p, :] = 0

        return outcome

    def determined_measurement(self, qubit_idx):
        """
        Set (2n+1)st row to 0
        Call rowsum
        Return r_{2n+1} as measurement outcome
        """
        n = self.n

        # Temporarily extend tableau by one row for the computation
        zero_row = Matrix(F2, [[0] * (2 * n + 1)])
        extended_tableau = self.tableau.stack(zero_row)

        # Temporarily replace tableau
        original_tableau = self.tableau
        self.tableau = extended_tableau

        try:
            extra_row = 2 * n  # The "2n+1"-th row

            # The extra row is already zero from stack operation
            # Apply rowsum with stabilizer row for measured qubit
            stabilizer_row = n + qubit_idx  # Row (n+p) in A&G notation
            self.rowsum(extra_row, stabilizer_row)

            # Return the phase bit of the extra row as the measurement outcome
            outcome = int(self.tableau[extra_row, 2 * n])

        finally:
            # Restore original tableau
            self.tableau = original_tableau

        return outcome

    def rowsum(self, i, p):
        """
        Replace row i with the product of rows i and p
        """

        n = self.n
        # extract rs
        r_i = self.tableau[i, 2 * n]  # phase bit for i
        r_p = self.tableau[p, 2 * n]  # phase bit for p

        phase = 0

        # extract xs and zs
        for k in range(n):
            x_i = self.tableau[i, k]
            z_i = self.tableau[i, k + n]
            x_p = self.tableau[p, k]
            z_p = self.tableau[p, k + n]

            # Turn (x, z) of i and p each into an int from 0 to 3 to index table
            # 0: I, 1: X, 2: Z, 3: Y
            mapping = {
                (0, 0): 0,  # I
                (1, 0): 1,  # X
                (0, 1): 2,  # Z
                (1, 1): 3,  # Y
            }
            a = mapping[(x_i, z_i)]
            b = mapping[(x_p, z_p)]

            # Phase table from Aaronson-Gottesman, defines number of i factors
            # Matrix entry: phase_table[a][b]
            phase_contrib = [
                [0, 0, 0, 0],
                [0, 0, 2, 2],
                [0, 2, 0, 2],
                [0, 2, 2, 0],
            ][a][b]

            phase += phase_contrib

        # Update phase bit of ith row
        self.tableau[i, 2 * n] = (r_i + r_p + (phase % 4) // 2) % 2

        # Now do XOR of X and Z parts of ith row
        for k in range(2 * n):
            self.tableau[i, k] += self.tableau[p, k]


def rowsum_by_row(a, b):
    """
    Given two rows a and b (numpy arrays of length 2n + 1),
    return the product row (a * b) as a new numpy array.
    This is used for testing associativity under rowsum.
    """

    n = (len(a) - 1) // 2  # number of qubits
    import numpy as np

    result = np.zeros_like(a)

    r_a = a[2 * n]  # phase bit of a
    r_b = b[2 * n]  # phase bit of b

    phase = 0

    mapping = {
        (0, 0): 0,  # I
        (1, 0): 1,  # X
        (0, 1): 2,  # Z
        (1, 1): 3,  # Y
    }

    phase_table = [
        [0, 0, 0, 0],
        [0, 0, 2, 2],
        [0, 2, 0, 2],
        [0, 2, 2, 0],
    ]

    for k in range(n):
        x_a = a[k]
        z_a = a[k + n]
        x_b = b[k]
        z_b = b[k + n]

        a_type = mapping[(x_a, z_a)]
        b_type = mapping[(x_b, z_b)]

        phase += phase_table[a_type][b_type]

    # Phase bit calculation same as original function
    new_phase = (r_a + r_b + (phase % 4) // 2) % 2
    result[2 * n] = new_phase

    # XOR of X and Z parts
    for k in range(2 * n):
        result[k] = a[k] ^ b[k]  # XOR operation

    return result


def is_associative(a, b, c):
    """
    Take 3 "rows" (np arrays not actual tableaux matrix rows) to test associativity of rowsum_by_row
    (a * b) * c == a * (b * c)
    """
    left = rowsum_by_row(rowsum_by_row(a, b), c)
    right = rowsum_by_row(a, rowsum_by_row(b, c))

    return np.array_equal(left, right)


def random_row(n):
    """
    Generate a random row of length 2n + 1:
    - First n bits: X bits (0 or 1)
    - Next n bits: Z bits (0 or 1)
    - Last bit: phase bit (0 or 1)
    """
    row = np.random.randint(0, 2, size=2 * n + 1, dtype=np.uint8)
    return row


def test_associativity_random(n, trials=1000):
    for _ in range(trials):
        a = random_row(n)
        b = random_row(n)
        c = random_row(n)
        if not is_associative(a, b, c):
            print("Associativity failed for:")
            print("a =", a)
            print("b =", b)
            print("c =", c)
            return False
    print(f"Associativity holds for all {trials} random triples under rowsum.")
    return True


if __name__ == "__main__":
    outcomes = []

    for i in range(0, 1000):
        tab = StabilizerTableau("-")
        outcome = tab.measure(0)
        outcomes.append(outcome)

    distrib = {"0": outcomes.count(0), "1": outcomes.count(1)}

    print(distrib)

    # Verify associativity

    n = 5  # number of qubits
    test_associativity_random(n, trials=10000)
