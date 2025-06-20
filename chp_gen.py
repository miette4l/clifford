"""
Generate all distinct states reachable from |+⟩ by applying H and S gates with a CHP simulator
This finds the orbit of |+⟩ under the Clifford group action
So is generating the group action as opposed to the abstract group

NB the stabilizer only hash -> 6 states
the full tableau hash -> 24 tableaux (with destabilizer->state degeneracy)
"""

from chp import StabilizerTableau
from itertools import product
from sage.all import Matrix, GF # type: ignore

F2 = GF(2)


class StateOrbitGenerator:
    def __init__(self, initial_state="+"):
        self.initial_state = initial_state
        self.generators = ["H", "S"]
        self.found_states = set()
        self.state_sequences = {}

    def apply_gate_sequence(self, gate_sequence):
        state = StabilizerTableau(self.initial_state)

        for gate in reversed(gate_sequence):
            state.conjugate(gate)

        # Hello rowsum??? Not needed until measurement??

        return state

    def find_stabilizer_state_name(self, tableau):
        """Try to identify if this tableau corresponds to a stabilizer state"""
        stabilizer_states = {
            "0": StabilizerTableau("0"),
            "1": StabilizerTableau("1"),
            "+": StabilizerTableau("+"),
            "-": StabilizerTableau("-"),
            "i": StabilizerTableau("i"),
            "-i": StabilizerTableau("-i"),
        }

        target_hash = hash(tableau)

        for name, stabilizer_tableau in stabilizer_states.items():
            if hash(stabilizer_tableau) == target_hash:
                return name

        return None

    def generate_orbit_brute_force(self, max_length=8):
        """
        Generate all reachable states by brute force up to max_length
        """
        print(f"Generating orbit of |{self.initial_state}⟩ under H and S gates...")

        found_states = set()
        state_info = []  # List of (sequence, tableau, state_name)

        # Start with initial state (empty sequence)
        initial_tableau = StabilizerTableau(self.initial_state)
        initial_hash = hash(initial_tableau)
        found_states.add(initial_hash)
        initial_name = self.find_stabilizer_state_name(initial_tableau)
        state_info.append(([], initial_tableau, initial_name))
        print(f"Initial state: |{self.initial_state}⟩ -> {initial_name}")

        # Generate sequences of increasing length
        for length in range(1, max_length + 1):
            print(f"\nChecking sequences of length {length}...")
            new_states_found = 0

            for sequence in product(self.generators, repeat=length):
                tableau = self.apply_gate_sequence(sequence)
                tableau_hash = hash(tableau)

                if tableau_hash not in found_states:
                    found_states.add(tableau_hash)
                    state_name = self.find_stabilizer_state_name(tableau)
                    state_info.append((list(sequence), tableau, state_name))

                    seq_str = "".join(sequence)
                    print(
                        f"  New state #{len(state_info)}: {seq_str} -> |{state_name}⟩"
                        if state_name
                        else f"  New state #{len(state_info)}: {seq_str} -> Unknown"
                    )
                    new_states_found += 1

            print(f"Found {new_states_found} new states at length {length}")
            print(f"Total unique states found: {len(state_info)}")

            # If no new states found, we've probs found the complete orbit
            if new_states_found == 0:
                print(
                    f"No new states found at length {length}. Orbit appears complete."
                )
                break

        return state_info

    def analyze_orbit_structure(self, state_info):
        """Print the structure of the found orbit"""
        print(f"Total states in orbit: {len(state_info)}")

        # Group by state name
        by_state_name = {}
        unknown_count = 0

        for sequence, tableau, state_name in state_info:
            if state_name:
                if state_name not in by_state_name:
                    by_state_name[state_name] = []
                by_state_name[state_name].append(sequence)
            else:
                unknown_count += 1

        print(f"Known stabilizer states found: {len(by_state_name)}")
        print(f"Unknown states: {unknown_count}")

        print("\nWays to reach each stabilizer state:")
        for state_name in sorted(by_state_name.keys()):
            sequences = by_state_name[state_name]
            print(f"  |{state_name}⟩: {len(sequences)} ways")

            # Show shortest sequences
            shortest = min(sequences, key=len)
            longest = max(sequences, key=len)
            shortest_str = "".join(shortest) if shortest else "I"
            longest_str = "".join(longest) if longest else "I"

            if len(sequences) <= 3:
                # Show all if few
                all_seqs = [("".join(seq) if seq else "I") for seq in sequences]
                print(f"    Sequences: {all_seqs}")
            else:
                # Show range if many
                print(f"    Shortest: {shortest_str}, Longest: {longest_str}")

        # Check if we found all 6 stabilizer single-qubit states
        expected_states = {"0", "1", "+", "-", "i", "-i"}
        # Corresponding to conjugation to Z, -Z, X, -X, Y, -Y
        found_states = set(by_state_name.keys())

        print(f"\nExpected stabilizer states: {expected_states}")
        print(f"Found stabilizer states: {found_states}")

        if found_states == expected_states:
            print("All 6 stabilizer single-qubit states are reachable!")
        else:
            missing = expected_states - found_states
            extra = found_states - expected_states
            if missing:
                print(f"Missing states: {missing}")
            if extra:
                print(f"Extra states: {extra}")

    def explore_complete_orbit(self):
        """Main method to explore the complete orbit"""
        print(f"Orbit exploration from |{self.initial_state}⟩")

        # Generate the complete orbit
        state_info = self.generate_orbit_brute_force()

        # Analyze the results
        self.analyze_orbit_structure(state_info)

        return state_info


class GroupActionTableGenerator:
    def __init__(self):
        self.states = ["0", "1", "+", "-", "i", "-i"]
        self.generators = ["H", "S"]

    def apply_sequence(self, sequence, state):
        """Apply gate sequence to a state, return result state"""
        tableau = StabilizerTableau(state)
        for gate in reversed(sequence):
            tableau.conjugate(gate)

        return tableau

        # Find which state this is
        for s in self.states:
            if hash(tableau) == hash(StabilizerTableau(s)):
                return s
        return None

    def sequence_to_permutation(self, sequence):
        """Convert sequence to how it permutes all 6 states"""
        return tuple(self.apply_sequence(sequence, s) for s in self.states)

    def generate_group(self, max_length=6):
        """Generate all distinct Clifford elements"""
        seen_perms = set()
        elements = []

        # Try all sequences up to max_length
        for length in range(max_length + 1):
            for seq in product(self.generators, repeat=length):
                perm = self.sequence_to_permutation(list(seq))

                if perm not in seen_perms:
                    seen_perms.add(perm)
                    seq_str = "".join(seq) if seq else "I"
                    elements.append((seq_str, perm))
                    print(f"{len(elements):2d}. {seq_str:8s} -> {perm}")

        print(f"\nFound {len(elements)} Clifford group elements")
        return elements

    def sequences_equal(self, seq1, seq2):
        """Check if two sequences produce the same result"""
        state = "+"
        s1 = self.apply_sequence(seq1, state)
        s2 = self.apply_sequence(seq2, state)
        return s1 == s2

    def check_relations(self):
        """Check H^2 = I and S^4 = I"""
        print("\nChecking group relations:")

        h2_perm = self.sequence_to_permutation(["H", "H"])
        identity_perm = self.sequence_to_permutation([])
        print(f"H^2 = I? {h2_perm == identity_perm}")

        s4_perm = self.sequence_to_permutation(["S", "S", "S", "S"])
        print(f"S^4 = I? {s4_perm == identity_perm}")

        # H and S don't commute
        hs_eq_sh = self.sequences_equal(["H", "S"], ["S", "H"])
        print(f"HS = SH: {hs_eq_sh} (should be False)")

    def find_order(self, sequence):
        for n in range(1, 25):
            seq = sequence * n
            if self.sequences_equal(seq, []):  # Empty sequence is identity
                print(f"Order of {sequence}: {n}")
                break

        return n

    def print_orders(self, elements):
        orders = {order: 0 for order in range(1, 5)}

        for element in elements:
            if element[0] == "I":
                orders[1] += 1
            else:
                order = gen.find_order(element[0])
                orders[order] += 1

    def print_carrying_table(self, state):
        sequences = ["H", "S", "SH", "HS", "HSH"]
        seq_to_bits = {}

        for seq_1 in sequences:  # seq_1 is row
            for seq_2 in sequences:  # seq_2 is column
                tab = self.apply_sequence(seq_1 + seq_2, state)
                seq_to_bits[seq_1, seq_2] = tab.tableau[0, -1], tab.tableau[1, -1]

        print("")
        print("")
        # Now make it into a table
        print(f"Carry Table: [r_d, r_s] from state |{state}> as row ∘ column")

        # Header
        header = f"{'σ':<6}"
        for seq in sequences:
            header += f"{seq:<8}"
        print(header)
        print("-" * len(header))

        phase_vectors = {}
        # Rows
        for seq_1 in sequences:
            row = f"{seq_1:<6}"
            for seq_2 in sequences:
                x_bit, z_bit = seq_to_bits[seq_1, seq_2]
                row += f"{z_bit}{x_bit}      "  # Inverted notation right now to match carry_table.py, swap to {x_bit}{z_bit} preferred
                phase_vectors[seq_1, seq_2] = Matrix(F2, [[z_bit, x_bit]])
            print(row)

        print("")
        print("")

        # Check the associativity of the subgroup of phase vectors

        results = []
        for a in phase_vectors.values():
            for b in phase_vectors.values():
                for c in phase_vectors.values():
                    results.append((((a + b)) +c) == (a +(b + c)))
        print(all(results))


if __name__ == "__main__":
    # Find the orbit of [H, S]

    explorer = StateOrbitGenerator("+")
    orbit = explorer.explore_complete_orbit()

    print(f"\nComplete orbit has {len(orbit)} distinct states.")

    # Orbit-stabilizer theorem - stabilizer order 4??
    # |Orbit| * |stabilizer| = |G|
    # 6 * 4 = 24

    # Generate a group action table
    # Representing Clifford operations as *permutations* of stabilizer states
    # Cayley's theorem??

    gen = GroupActionTableGenerator()

    elements = gen.generate_group()

    gen.check_relations()

    gen.print_carrying_table("+")
