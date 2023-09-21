import matplotlib.pyplot as plt
import numpy as np
from math import pi, atan2, sqrt
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, transpile
from qiskit.tools.visualization import circuit_drawer
from qiskit.quantum_info import state_fidelity, OneQubitEulerDecomposer
from qiskit import BasicAer
from qiskit.circuit.library import QFT
from qiskit.visualization import plot_histogram
from qiskit_aer import AerSimulator
from qiskit import Aer
from qiskit.extensions import UnitaryGate
from callGridSynth import *
from qubitdecomp import *
from CTcompiler import *
import scipy

import warnings
warnings.filterwarnings("ignore")

q = QuantumRegister(3)
qc = QuantumCircuit(q)
qc.x(0)
qc.crx(pi/3,0,1)
#qc.cx(1,2)
qc2 = compile(qc,1e-8)




measurement = QuantumCircuit(3,3)
measurement.barrier(range(3))
measurement.measure(range(3),range(3))

circuit1 = qc.compose(measurement)
circuit2 = qc2.compose(measurement)


aer_sim = Aer.get_backend('aer_simulator')
circuit1_compiled = transpile(circuit1, aer_sim)
runs = 2000;
job_sim = aer_sim.run(circuit1_compiled, shots=runs)
result_sim = job_sim.result()
counts = result_sim.get_counts(circuit1_compiled)

print("Results for exact rotation: "+str(counts))

aer_sim = Aer.get_backend('aer_simulator')
circuit2_compiled = transpile(circuit2, aer_sim)
runs = 2000;
job_sim = aer_sim.run(circuit2_compiled, shots=runs)
result_sim = job_sim.result()
counts = result_sim.get_counts(circuit2_compiled)

print("Results for approx rotation: "+str(counts))
