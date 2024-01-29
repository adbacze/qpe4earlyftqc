import subprocess
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, transpile
import matplotlib.pyplot as plt
import numpy as np

#gridsynth can be downloaded here: https://www.mathstat.dal.ca/~selinger/newsynth/#downloading

#Purpose of this file is to call gridsynth, which is written in haskell, using python functions.
#gridsynth provides multiple means for specifying precision of synthesis, including specifying digits of accuracy in binary or decimal
#or specifying an epsilon such that the synthesised angle is am epsilon approximation of the true angle

#Functions with GridSynth suffix take an angle and a desired precision as inputs and return a string specifying the
#sequence of gates which approximates an Z-rotation (up to a global phase) by specified angle to specified precision.
#The prefix on these functions indicates how the precision is specified.

#Interface with qiskit is done in functions with "RZ" suffix
#These functions take a qiskit circuit object, an int specifying a qubit in said circuit object, as well as an angle and a desired precision.
#The function calls gridsynth to synthesize rotation then adds all gates to circuit
#The prefix on these functions indicates how the precision is specified.

#There is an additional function named "epsilonClCRZ". This function allows for the synthesis of a classically controlled rotation
#as used in IPE. This is done simply by applying a classical control to every gate

def decimalGridSynth(angle, precision):
    if (angle < 0): #correct formatting for gridsynth with negative angle
        angle = '"('+str(angle)+')"'
    cmd = f"~/mystuff/python_scripts/gridsynth "+str(angle)+" -d "+str(precision);
    command_output = subprocess.check_output(cmd, shell=True).decode('utf-8')
    command_output = command_output[:-1] # Drop the trailing newline
    sequence = command_output.replace('W','')# drop the global phase
    sequence = sequence[::-1] # reverse to go from operator to circuit form

    return sequence

def decimalRZ(circuit,qubit,angle,precision): #

    gateSequence = decimalGridSynth(angle,precision) #obtain string specifying gate sequence

    for i in range(0,len(gateSequence)): #parse through string, applying each gate to specified qubit within qiskit circuit
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
    if (angle < 0): #correct formatting for gridsynth with negative angle
        angle = '"('+str(angle)+')"'
    cmd = f"~/mystuff/python_scripts/gridsynth "+str(angle)+" -b "+str(precision);
    command_output = subprocess.check_output(cmd, shell=True).decode('utf-8')
    command_output = command_output[:-1] # Drop the trailing newline
    sequence = command_output.replace('W','')# drop the global phase
    sequence = sequence[::-1] # reverse to go from operaor to circuit form

    return sequence

def binaryRZ(circuit,qubit,angle,precision): #

    gateSequence = binaryGridSynth(angle,precision) #obtain string specifying gate sequence

    for i in range(0,len(gateSequence)): #parse through string, applying each gate to specified qubit within qiskit circuit
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

    if (angle < 0): #correct formatting for gridsynth with negative angle
        angle = '"('+str(angle)+')"'

    cmd = f"~/mystuff/python_scripts/gridsynth "+str(angle)+" -e "+str(precision);
    command_output = subprocess.check_output(cmd, shell=True).decode('utf-8')
    command_output = command_output[:-1] # Drop the trailing newline
    sequence = command_output.replace('W','')# drop the global phase
    sequence = sequence[::-1] # reverse to go from operaor to circuit form

    return sequence

def epsilonRZ(circuit,qubit,angle,precision): #

    gateSequence = epsilonGridSynth(angle,precision) #obtain string specifying gate sequence

    for i in range(0,len(gateSequence)): #parse through string, applying each gate to specified qubit within qiskit circuit
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

def epsilonClCRZ(circuit,qubit,angle,precision,cbit,cvalue): #Add classical control by adding control to all operations
                                                             #cbit/cvalue designate classical bit(s) and value to controll on
    gateSequence = epsilonGridSynth(angle,precision) #obtain string specifying gate sequence

    for i in range(0,len(gateSequence)): #parse through string, applying each gate controlled on specified ClassicalRegister to specified qubit within qiskit circuit
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
