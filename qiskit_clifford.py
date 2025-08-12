from qiskit import QuantumCircuit
from qiskit.quantum_info import Clifford

import numpy as np

n=2

qc = QuantumCircuit(n)
qc.h(0)
qc.cx(0, 1)
cliff = Clifford(qc)

# print(cliff)

# Print the Clifford destabilizer rows
# print(cliff.to_labels(mode="D"))
 
# # Print the Clifford stabilizer rows
# print(cliff.to_labels(mode="S"))

# print(cliff.tableau.astype(int))

qc = QuantumCircuit(1)
qc.x(0)
cliff = Clifford(qc)

print(cliff)
print(cliff.tableau.astype(int))