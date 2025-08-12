import stim

GEN_TO_STIM = {
    'h': "H 0",
    's': "S 0",
    'i': "I 0"
}

def word_to_stim_circuit(word):
    """
    Convert a generator word like ['h','s','h'] into a stim.Circuit.
    """
    circuit = stim.Circuit()
    for g in word:
        circuit.append_from_stim_program_text(GEN_TO_STIM[g])
    return circuit

def stim_tableau_for_word(word):
    sim = stim.TableauSimulator()
    circuit = word_to_stim_circuit(word)
    tab = circuit.to_tableau()
    sim.do_circuit(circuit)
    print("The tableau for the circuit:")
    print(tab)
    print("The tableau for the state:")
    print(sim.current_inverse_tableau())
    assert tab == sim.current_inverse_tableau()
    return sim.current_inverse_tableau()

def words_equivalent_stim(word1, word2):
    return stim_tableau_for_word(word1) == stim_tableau_for_word(word2)

if __name__ == "__main__":
    print(words_equivalent_stim("ssss", "i")) # equiv ZZ
    print(words_equivalent_stim("hsshhssh", "i")) # equiv XX
    print(words_equivalent_stim("sshsshsshssh", "i")) # equiv YY