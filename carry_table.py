"""
Automate generation of carry table from subsequent conjugations
Use Pauli string rules
"""

conjugation_rules = {
    "HXH": "Z",
    "HYH": "-Y",
    "HZH": "X",
    "SXS": "Y",
    "SYS": "-X",
    "SZS": "Z",
}


def conjugate(p: str, c: str):
    """
    p can be "X", "Y" or "Z"
    c can be "H" or "S"
    """
    conjugation = str(c + p + c[::])
    result = conjugation_rules[conjugation]
    if "-" in result:
        return result[1:], 1
    return result, 0


def chain(p: str, c_chain: str):
    """
    p can be "X", "Y" or "Z"
    c can be a chain of "H"s and/or "S"s
    """
    label = p
    overall_sign = 0
    for c in reversed(c_chain):
        label, sign = conjugate(label, c)
        overall_sign ^= sign
    return label, overall_sign


def generate_table_element(c_chain: str):
    label, phase_x = chain("X", c_chain)
    label_z, phase_z = chain("Z", c_chain)
    return phase_x, phase_z


def generate_table():
    headers = [
        "H",
        "S",
        "SH",
        "HS",
        "HH",
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
    headers = headers[:7]
    table = {}

    for row, row_header in enumerate(headers):
        for col, column_header in enumerate(headers):
            c_chain = row_header + column_header  # row dot column
            table[row, col] = generate_table_element(c_chain)

    # Print header
    print("     " + "  ".join(f"{h:>4}" for h in headers))
    print("   +" + "------" * len(headers))

    # Print rows
    for row, row_header in enumerate(headers):
        row_str = f"{row_header:>2} |"
        for col in range(len(headers)):
            phase_x, phase_z = table[(row, col)]
            cell = f"{phase_x}{phase_z}"
            row_str += f"  {cell:>4}"
        print(row_str)

    return table


if __name__ == "__main__":
    table = generate_table()
    group = [element for element in table.values()]

    # Cross-testing with unitary_gen.py to see if chaining is right

    def track_evolution(c_chain: str):
        label_x, phase_x = chain("X", c_chain)
        label_z, phase_z = chain("Z", c_chain)
        return (label_x, phase_x), (label_z, phase_z)

    pauli_frame = track_evolution("HS")
    # HS here corresponds to matrix HS unitarily acting on state vectors X and Z

    print(pauli_frame)
