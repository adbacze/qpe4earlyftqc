import matplotlib.pyplot as plt
import numpy as np
from math import pi, atan2, sqrt
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, transpile
from qiskit.quantum_info.operators import Operator, Pauli
from qiskit.extensions import RXGate, XGate, CXGate
from qiskit.tools.visualization import circuit_drawer
from qiskit.quantum_info import state_fidelity, TwoQubitBasisDecomposer, Operator, OneQubitEulerDecomposer
from qiskit import BasicAer
from qiskit.circuit.library import QFT
from qiskit.visualization import plot_histogram
from qiskit_aer import AerSimulator
from qiskit import Aer
from qiskit_aer.noise import (NoiseModel, QuantumError, ReadoutError,
    pauli_error, depolarizing_error, thermal_relaxation_error)
from qiskit.extensions import UnitaryGate
from qiskit.quantum_info.random import random_unitary
import scipy
import math
from CliffordTCompiler.callGridSynth import *
from CliffordTCompiler.CTcompiler import *
from qpe_protocols.aRPE import *
from qpe_protocols.qcels import *
from qpe_protocols.MMqcels import *
