from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit.quantum_info.operators import Operator, Pauli
from qiskit.extensions import RXGate, XGate, CXGate
from qiskit.quantum_info import TwoQubitBasisDecomposer, Operator, OneQubitEulerDecomposer
from callGridSynth import *

def reduce(circuit): #Take circuit object and reduce it to cliffords and single qubit rotations
    cx = Operator(CXGate())
    UZXZ = TwoQubitBasisDecomposer(cx, euler_basis="ZXZ")
    ZXZ = OneQubitEulerDecomposer(basis = "ZXZ")

    G = circuit.data
    newCircuit = QuantumCircuit(circuit.width())

    for i in range(len(G)):

        if (G[i][0].num_qubits == 1):
            qubit1 = G[i][1][0].index
            if (G[i][0].name == "h"): #If circuit instruction is already a clifford leave as is
                newCircuit.h(qubit1)
            elif (G[i][0].name == "s"):
                newCircuit.s(qubit1)
            elif (G[i][0].name == "t"):
                newCircuit.t(qubit1)
            else:
                temp = ZXZ(G[i][0]) #Non clifford single qubit gates get Euler decomposed
                newCircuit = newCircuit.compose(temp,[qubit1])

        elif (G[i][0].num_qubits == 2):
            temp = UZXZ(G[i][0])
            qubit1 = G[i][1][0].index
            qubit2 = G[i][1][1].index
            newCircuit = newCircuit.compose(temp,[qubit1,qubit2])

    return newCircuit

def compileCT(circuit, precision):
    ctCircuit = QuantumCircuit(circuit.width())
    temp = reduce(circuit)
    temp = temp.data

    for i in range(len(temp)):
        if (temp[i][0].name == "h"):
            qubit = temp[i][1][0].index
            ctCircuit.h(qubit)
        elif (temp[i][0].name == "s"):
            qubit = temp[i][1][0].index
            ctCircuit.s(qubit)
        elif (temp[i][0].name == "t"):
            qubit = temp[i][1][0].index
            ctCircuit.t(qubit)
        elif (temp[i][0].name == 'rz'):
            angle = temp[i][0].params
            qubit = temp[i][1][0].index
            if (qubit == ' '):
                print(temp[i])
            #ctCircuit.rz(angle[0],int(qubit))
            epsilonRZ(ctCircuit,int(qubit),angle[0],precision/pow(len(temp),1/2))
        elif (temp[i][0].name == 'rx'):
            angle = temp[i][0].params
            qubit = temp[i][1][0].index
            #print(qubit)
            ctCircuit.h(int(qubit))
            epsilonRZ(ctCircuit,int(qubit),angle[0],precision/pow(len(temp),1/2))
            ctCircuit.h(int(qubit))
        elif (temp[i][0].name == 'unitary'):
            control = temp[i][1][0].index
            target = temp[i][1][1].index
            ctCircuit.cx(int(control),int(target))

    return ctCircuit
