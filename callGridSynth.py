import subprocess
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, transpile
import matplotlib.pyplot as plt
import numpy as np

#https://www.mathstat.dal.ca/~selinger/newsynth/#downloading

def decimalGridSynth(angle, precision):

    cmd = f"gridsynth "+str(angle)+" -d "+str(precision);
    command_output = subprocess.check_output(cmd, shell=True).decode('utf-8')
    command_output = command_output[:-1] # Drop the trailing newline
    sequence = command_output.replace('W','')# drop the global phase
    sequence = sequence[::-1] # reverse to go from operaor to circuit form

    return sequence

def decimalRZ(circuit,qubit,angle,precision): #

    gateSequence = decimalGridSynth(angle,precision)

    for i in range(1,len(gateSequence)):
        if (gateSequence[i] == "H"):
            circuit.h(qubit)
        elif (gateSequence[i] == "T"):
            circuit.t(qubit)
        elif (gateSequence[i] == "S"):
            circuit.s(qubit)
        elif (gateSequence[i] == "X"):
            circuit.x(qubit)
        elif (gateSequence[i] == "Z"):
            circuit.z(qubit)

def binaryGridSynth(angle, precision):

    cmd = f"gridsynth "+str(angle)+" -b "+str(precision);
    command_output = subprocess.check_output(cmd, shell=True).decode('utf-8')
    command_output = command_output[:-1] # Drop the trailing newline
    sequence = command_output.replace('W','')# drop the global phase
    sequence = sequence[::-1] # reverse to go from operaor to circuit form

    return sequence

def binaryRZ(circuit,qubit,angle,precision): #

    gateSequence = binaryGridSynth(angle,precision)

    for i in range(1,len(gateSequence)):
        if (gateSequence[i] == "H"):
            circuit.h(qubit)
        elif (gateSequence[i] == "T"):
            circuit.t(qubit)
        elif (gateSequence[i] == "S"):
            circuit.s(qubit)
        elif (gateSequence[i] == "X"):
            circuit.x(qubit)
        elif (gateSequence[i] == "Z"):
            circuit.z(qubit)

def epsilonGridSynth(angle, precision):
    #print(precision)
    if (abs(angle) < precision):
        angle = 0

    if (angle < 0):
        angle = '"('+str(angle)+')"'

    cmd = f"gridsynth "+str(angle)+" -e "+str(precision);
    command_output = subprocess.check_output(cmd, shell=True).decode('utf-8')
    command_output = command_output[:-1] # Drop the trailing newline
    sequence = command_output.replace('W','')# drop the global phase
    sequence = sequence[::-1] # reverse to go from operaor to circuit form

    return sequence

def epsilonRZ(circuit,qubit,angle,precision): #

    gateSequence = epsilonGridSynth(angle,precision)

    for i in range(1,len(gateSequence)):
        if (gateSequence[i] == "H"):
            circuit.h(qubit)
        elif (gateSequence[i] == "T"):
            circuit.t(qubit)
        elif (gateSequence[i] == "S"):
            circuit.s(qubit)
        elif (gateSequence[i] == "X"):
            circuit.x(qubit)
        elif (gateSequence[i] == "Z"):
            circuit.z(qubit)

def epsilonClCRZ(circuit,qubit,angle,precision,cbit,cvalue): #Add classical control by adding control
                                                             #to all operations
    gateSequence = epsilonGridSynth(angle,precision)

    for i in range(1,len(gateSequence)):
        if (gateSequence[i] == "H"):
            circuit.h(qubit).c_if(cbit,cvalue)
        elif (gateSequence[i] == "T"):
            circuit.t(qubit).c_if(cbit,cvalue)
        elif (gateSequence[i] == "S"):
            circuit.s(qubit).c_if(cbit,cvalue)
        elif (gateSequence[i] == "X"):
            circuit.x(qubit).c_if(cbit,cvalue)
        elif (gateSequence[i] == "Z"):
            circuit.z(qubit).c_if(cbit,cvalue)
