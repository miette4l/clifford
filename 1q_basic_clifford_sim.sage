"""
We want to:
- Represent arbitrary 1-qubit Pauli stabilizers as binary symplectic vectors in a tableau
- Represent any 1-qubit Clifford gate as a Clifford tableau with phase
- Apply that Clifford gate to a stabilizer (multiply tableau)
"""

from abc import ABC, abstractmethod

F2 = GF(2)

CLIFFORD_LOOKUP = {
            "I": Matrix(F2, [[1, 0, 0, 0], [0, 1, 0, 0]]),  # X -> X, Z -> Z
            "X": Matrix(F2, [[1, 0, 0, 0], [0, 1, 1, 0]]),  # Z -> -Z
            "Z": Matrix(F2, [[1, 0, 1, 0], [0, 1, 0, 0]]),  # X -> -X
            "H": Matrix(F2, [[0, 1, 0, 0], [1, 0, 0, 0]]),  # X <-> Z
            "S": Matrix(F2, [[1, 1, 0, 0], [0, 1, 0, 0]]),  # X -> Y
        }


class Operator():
    """
    Objects to use mainly for testing and conversion for printing. 
    Here I will store the matrices for testing the simulator is equivalent to usual matrix multiplication methods.
    And the strings so I can pretty print.
    """

    def __init__(self, name, matrix):
        self.name = name # Pauli string or name of gate
        self.matrix = matrix
        self.pauli = False
        self.clifford = False


class Tableau(ABC):
    """
    [ x₁ z₁ r₁ i₁ ] where 1==presence and 0==absence
    [ x₂ z₂ r₂ i₂ ]
    """

    def __init__(self):
        self.tableau = None
        self.xz = None
        self.r = None
        self.i = None

    def split_tableau(self):
        self.xz = self.tableau[:, :-2]
        self.r = self.tableau[:, -2]
        self.i = self.tableau[:, -1]

    def recombine_tableau(self):
        self.tableau = self.xz.augment(self.r).augment(self.i)
    
    def print_bin(self):
        """
        Print the binary symplectic form tableau
        """
        print(self.tableau)

    def print_string(self):
        """
        Print as tensor product pauli operators for more readability
        """
        pass


class CliffordTableau(Tableau):
    def __init__(self, clifford: str):
        super().__init__()
        self.n = 1

        try:
            self.tableau = CLIFFORD_LOOKUP.get(clifford)
        except KeyError:
            print("Not (yet?) supported")
            return False
        
        self.split_tableau()

    def print_bin(self):
        print("Clifford tableau (binary)")
        super().print_bin()


class StabilizerTableau(Tableau):
    def __init__(self, stabilizer_state: str):
        super().__init__()
        self.n = 1
        if stabilizer_state == "0":
            self.tableau = Matrix(F2, [[0, 1, 0, 0], [0, 0, 0, 0]])
        else:
            print("Not ready")

        self.split_tableau()

    def print_bin(self):
        print("Stabilizer tableau (binary)")
        super().print_bin()

    def conjugate(self, clifford_tableau: CliffordTableau):

        # Do the signs and is
        # New r: r' = r + (x_old * r_C[0]) + (z_old * r_C[1])
        x_old = self.xz[0, 0]
        z_old = self.xz[0, 1]

        for row in range(self.xz.nrows()):
            x = self.xz[row, 0]
            z = self.xz[row, 1]

            # Binary mod-2 addition of original and gate-induced phases (bleh this is so long)
            self.r[row] = (self.r[row] + x * clifford_tableau.r[0] + z * clifford_tableau.r[1]) % 2
            self.i[row] = (self.i[row] + x * clifford_tableau.i[0] + z * clifford_tableau.i[1]) % 2

        # Do the xs and zs
        self.xz = self.xz * clifford_tableau.xz

        self.recombine_tableau()

        return self.tableau

class StabilizerState():
    def __init__(self, state: str):
        self.r = None
        self.i = None


if __name__ == "__main__":
    stab_tab = StabilizerTableau("0")
    stab_tab.print_bin()
    x_tab = CliffordTableau("X")
    x_tab.print_bin()
    stab_tab.conjugate(x_tab)
    stab_tab.print_bin()