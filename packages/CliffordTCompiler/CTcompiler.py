from callGridSynth import *

#This script provides functions to take qiskit circuit objects and compile them into circuits containing only Clifford+T gates
#Note that since this utilizes qiskit's built in Decomposer functions input circuits cannot contain classical registers.
#If classical registers are needed (i.e. for measurement), you can first compile a quantum circuit containing no classical registers,
#then merge the compiled circuit with a quantum circuit containing the classical register using qiskit's "compose" function
#Additionally this compler currently only works for gates acting on no more than two qubits, gates acting on more than 2 qubits
#are simply ignored and functions will output compiled circuits without those gates

#There are two compile functions: compileCT and compileCT2. Each take a quantum circuit object and return the equivalent
#Clifford+T within desired precision. The difference is in how synthesis precision is set, compileCT will compile the entire
#circuit within set precision assuming linear growth in synthesis error amongst all Z-rotations.

#compileCT2 synthesizes all rotations to the specified precision. In order to ensure the correct global error, the total number
#of Z-rotations needs to known a-priori so the correct precision can be specified.
#compileCT2 is useful when compiling circuits containing mid-circuit measurements, classical controls or circuits that are part
#of a algorithm containing many circuits.

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
            elif (G[i][0].name == "z"):
                newCircuit.z(qubit1)
            elif (G[i][0].name == "x"):
                newCircuit.x(qubit1)
            else:
                temp = ZXZ(G[i][0]) #Non clifford single qubit gates get Euler decomposed
                newCircuit = newCircuit.compose(temp,[qubit1])

        elif (G[i][0].num_qubits == 2):
            temp = UZXZ(G[i][0])
            qubit1 = G[i][1][0].index
            qubit2 = G[i][1][1].index
            newCircuit = newCircuit.compose(temp,[qubit1,qubit2])

    return newCircuit

def compileCT(circuit, precision): #precision sets precision for entire circuit assuming linear growth in synthesis error
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
        elif (temp[i][0].name == "x"):
            qubit = temp[i][1][0].index
            ctCircuit.x(qubit)
        elif (temp[i][0].name == "z"):
            qubit = temp[i][1][0].index
            ctCircuit.z(qubit)
        elif (temp[i][0].name == 'rz'):
            angle = temp[i][0].params
            qubit = temp[i][1][0].index
            epsilonRZ(ctCircuit,int(qubit),angle[0],precision/len(temp))
        elif (temp[i][0].name == 'rx'):
            angle = temp[i][0].params
            qubit = temp[i][1][0].index
            ctCircuit.h(int(qubit))
            epsilonRZ(ctCircuit,int(qubit),angle[0],precision/len(temp))
            ctCircuit.h(int(qubit))
        elif (temp[i][0].name == 'unitary'):
            control = temp[i][1][0].index
            target = temp[i][1][1].index
            ctCircuit.cx(int(control),int(target))

    return ctCircuit

def compileCT2(circuit, precision): #precision sets synthesis for each individual Z-rotation
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
        elif (temp[i][0].name == "x"):
            qubit = temp[i][1][0].index
            ctCircuit.x(qubit)
        elif (temp[i][0].name == "z"):
            qubit = temp[i][1][0].index
            ctCircuit.z(qubit)
        elif (temp[i][0].name == 'rz'):
            angle = temp[i][0].params
            qubit = temp[i][1][0].index
            epsilonRZ(ctCircuit,int(qubit),angle[0],precision)
        elif (temp[i][0].name == 'rx'):
            angle = temp[i][0].params
            qubit = temp[i][1][0].index
            ctCircuit.h(int(qubit))
            epsilonRZ(ctCircuit,int(qubit),angle[0],precision)
            ctCircuit.h(int(qubit))
        elif (temp[i][0].name == 'unitary'):
            control = temp[i][1][0].index
            target = temp[i][1][1].index
            ctCircuit.cx(int(control),int(target))

    return ctCircuit
