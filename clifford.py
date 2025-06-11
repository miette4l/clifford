"""
Run using `sage clifford.py`

1) We want to:
- Represent arbitrary 1-qubit stabilizers as binary symplectic vectors in a tableau
- Represent any 1-qubit Clifford gate as a Clifford tableau with stabilizer phases tracked
- Conjugate the arbitrary stabilizer under the arbitrary Clifford

2) We then want to:
- Generate the Clifford group for 1 qubit as Clifford tableaus via composition
- This means being able to compose Clifford tableaus and properly handle stabilizer phases

3) Next I might want to:
- Generate the traditional Clifford group for 2 qubits by adding a C-gate?
"""

from abc import ABC, abstractmethod
from sage.all import Matrix, GF, UniversalCyclotomicField  # type: ignore


F2 = GF(2)
roots = UniversalCyclotomicField().gen(8)

# Each Clifford gate is written in tableau format expressing its conjugations of X and Z
CLIFFORD_LOOKUP = {
            "I": Matrix(F2, [[1, 0, 0, 0], [0, 1, 0, 0]]),  # X -> X, Z -> Z
            "X": Matrix(F2, [[1, 0, 0, 0], [0, 1, 1, 0]]),  # X -> X, Z -> -Z
            "Y": Matrix(F2, [[1, 0, 0, 0], [0, 1, 0, 0]]),  # X -> X, Z -> Z
            "Z": Matrix(F2, [[1, 0, 1, 0], [0, 1, 0, 0]]),  # X -> -X, Z -> Z
            "H": Matrix(F2, [[0, 1, 0, 0], [1, 0, 0, 0]]),  # X <-> Z
            "S": Matrix(F2, [[1, 1, 0, 1], [0, 1, 0, 0]]),  # X -> Y, Z -> Z
        }

# Each stabilizer state is expressed by its stabilizer
STABILIZER_LOOKUP = {
    "0": Matrix(F2, [[0, 1, 0, 0]]),   # Stabilizer: Z
    "1": Matrix(F2, [[0, 1, 1, 0]]),   # Stabilizer: -Z
    "+": Matrix(F2, [[1, 0, 0, 0]]),   # Stabilizer: X
    "-": Matrix(F2, [[1, 0, 1, 0]]),   # Stabilizer: -X
    "i": Matrix(F2, [[1, 1, 0, 1]]),   # Stabilizer: Y
    "-i": Matrix(F2, [[1, 1, 1, 1]]),  # Stabilizer: -Y
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
    [ x₁ z₁ r₁ i₁ ]
    [ x₂ z₂ r₂ i₂ ]
    """

    def __init__(self):
        self.tableau = None
        self.xz = None
        self.r = None
        self.i = None

    def __eq__(self, other):
        """
        For direct equality checking between Tableau objects
        """
        return self.tableau == other.tableau

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
        Print as Pauli strings for more readability
        """
        pass


class CliffordTableau(Tableau):
    def __init__(self, clifford: str):
        super().__init__()
        self.n = 1

        try:
            self.tableau = CLIFFORD_LOOKUP[clifford]
        except KeyError:
            print("Not (yet?) supported")
            return False
        
        self.split_tableau()

    def __hash__(self):
        """
        Returns the hash value of a Clifford tableau.
        """
        return hash(tuple(tuple(row) for row in self.tableau.rows()))
    
    def compose(self, other):
        """
        Compose two Clifford tableaux: self ∘ other (other applied first, then self so R->L)
        """
        new_tableau = CliffordTableau("I")
        
        new_tableau.xz = other.xz * self.xz
        
        # For each generator in other (X and Z)
        for i in range(2):
            # Start with base phases from other
            phase_r = other.r[i, 0]
            phase_i = other.i[i, 0]
            
            # Get what 'other' transforms this generator to
            other_x = other.xz[i, 0] 
            other_z = other.xz[i, 1]
            
            # Add phases from self's transformation
            if other_x:  # If other has X component
                phase_r = (phase_r + self.r[0, 0]) % 2
                phase_i = (phase_i + self.i[0, 0]) % 2
                
            if other_z:  # If other has Z component  
                phase_r = (phase_r + self.r[1, 0]) % 2
                phase_i = (phase_i + self.i[1, 0]) % 2
            
            # Special case for anticommutation: if other maps to Y and self is Hadamard
            if other_x and other_z and self.xz == Matrix(F2, [[0,1],[1,0]]):
                phase_r = (phase_r + 1) % 2  # H(Y) = -Y
                
            # Another special case: i^2=-1
            if other.i[i, 0] and self.i[i, 0]:
                phase_r = (phase_r + self.i[i, 0]) % 2

            # Set final phases
            new_tableau.r[i, 0] = phase_r
            new_tableau.i[i, 0] = phase_i
        
        new_tableau.recombine_tableau()
        return new_tableau


    def print_bin(self):
        super().print_bin()


class StabilizerTableau(Tableau):
    def __init__(self, stabilizer_state: str):
        super().__init__()
        self.n = 1
        try:
            self.tableau = STABILIZER_LOOKUP[stabilizer_state]
        except KeyError:
            print(f"Unknown stabilizer state '{stabilizer_state}'")
            return False

        self.split_tableau()

    def print_bin(self):
        super().print_bin()

    def conjugate(self, clifford_tableau: CliffordTableau):
        
        for row in range(self.xz.nrows()):
            x = self.xz[row, 0]
            z = self.xz[row, 1]

            # Do the rs and is
            # First combining any is (i^2=-1 case)
            if self.i[row, 0] and clifford_tableau.i[row, 0]:
                self.r[row, 0] += clifford_tableau.i[row, 0]
            
            # Then general stuff
            self.r[row, 0] = (self.r[row, 0] + x * clifford_tableau.r[0, 0] + z * clifford_tableau.r[1, 0]) % 2
            self.i[row, 0] = (self.i[row, 0] + x * clifford_tableau.i[0, 0] + z * clifford_tableau.i[1, 0]) % 2  
            
            # Finally X and Z anticommutation
            if clifford_tableau.xz == Matrix(F2, [[0,1],[1,0]]) and x and z:
                self.r[row, 0] = (self.r[row, 0] + 1) % 2

        # Do the xs and zs
        self.xz = self.xz * clifford_tableau.xz

        self.recombine_tableau()

        return self.tableau


class StabilizerState():
    def __init__(self, stab_tab: StabilizerTableau):
        self.tableau = stab_tab
        self.global_r = 0
        self.global_i = 0


def generate_Clifford_group():
    """
    Generate the Clifford group for 1 qubit from H, S and I^{1/8}.
    """
    generators = [CliffordTableau("H"), CliffordTableau("S")]

    group = set(generators)
    added = True

    while added:
        added = False
        new_elements = set()

        for g1 in group:
            for g2 in generators:
                composed = g1.compose(g2)
                if composed not in group:
                    new_elements.add(composed)

        if new_elements:
            group.update(new_elements)
            added = True

    return group


if __name__ == "__main__":
    grwp = generate_Clifford_group()
    print(len(grwp))
    eighth_roots = [roots**k for k in range(8)]
    print(eighth_roots)

    big_group = set()

    # Each object will have to be an eighth ROU or scalar 1->8
    
    for element in grwp:
        for i, root in enumerate(eighth_roots):
            scaled_element = (element, i+1)
            big_group.add(scaled_element)

    print(len(big_group))