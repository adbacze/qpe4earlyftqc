# qpe4earlyftqc
This is part of a project to evaluate resource requirements of quantum phase estimation (qpe) protocols in an early fault tolerant setting. 
Files inculde scripts to run various qpe protocols using IBMs qiskit package, along with a "Clifford+T Compilier" with takes a quantum circuit and synthesizes it into the universal gate set of Cliffords + Ts.
Additionally, subdirectories H2_simulations and estimator_bias_tests contain the scripts used to run all simulations described in **ADD LINK TO ARXIV**

In order for the Clifford+T Compilier to work, gridsynth must be downloaded (https://www.mathstat.dal.ca/~selinger/newsynth/) and placed in file ./packages/C

Files are structured such that scripts should be run from parent directory (e.g. to run "H2ipe.py" either run as "$py ./H2_simulations/H2ipe" ot first move the script from /H2_simulations into parent directory)  

All code was developed using Python 3.9.12. As well as following qiskit versions: 

{'qiskit-terra': '0.24.1', 'qiskit-aer': '0.12.0', 'qiskit-ignis': None, 'qiskit-ibmq-provider': '0.20.2', 'qiskit': '0.43.1', 'qiskit-nature': None, 'qiskit-finance': None, 'qiskit-optimization': None, 'qiskit-machine-learning': None}
