This folder contains files to take qiskit QuantumCircuit objects and "compile" them into their equivalent Clifford+T circuits.
Workflow is parse through a QuantumCircuit and decompose all 1 and 2 qubit gates into Clifford gates and Z rotations (currently does not work with gates acting on more than 2 qubits).
All Z rotations are then passed to gridsynth (https://www.mathstat.dal.ca/~selinger/newsynth/) and gate sequence approximating the rotation is added to the circuit.
These scripts are built to run on Linux, minor changes may be needed to run on Windows as well as the Windows version of gridsynth.

The script "callGridSynth" uses the OS to call gridsynth.
It provides the path to gridsynth starting from the parent directory (qpe4earlyftqc) to this subdirectory.
Need to be careful that path is correct.
