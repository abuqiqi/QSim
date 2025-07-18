#pragma once

#include "qcircuit.h"
#include "omsim.h"

#define random(x) (rand()%x)

//
// Quantum circuit generators
//
vector<QGate> testkctrls(int k);
vector<QGate> testkqubits(int k);
vector<QGate> testkfusion(int k);

QCircuit test(int numQubits);
QCircuit QFT(int numQubits);
QCircuit QFT_Quirk(int numQubits);
QCircuit Adder(int numQubits);
QCircuit QAOA(int numQubits);
QCircuit Grover(int numQubits, int k);
QCircuit IQP(int numQubits);
QCircuit VQC(int numQubits);
QCircuit VQC1(int numQubits);
QCircuit VQC2(int numQubits); // circuit 5
QCircuit VQC_NN(int numQubits);
QCircuit VQC_CB(int numQubits);
QCircuit VQC_AA(int numQubits);
QCircuit Hyperbolic(int numQubits);
QCircuit QPE(int numQubits);

// Generate a circuit with regular state vector using H, Z, X, CX gates
QCircuit RandomRegular(int numQubits, int numDepths);

// Generate a circuit with relatively random state vector using H, Z, X, RY, CX gates
QCircuit RandomMedium(int numQubits, int numDepths);

// Generate a circuit with random state vector using H, Z, X, RZ, CX gates
QCircuit RandomRX(int numQubits, int numDepths);
QCircuit RandomRZ(int numQubits, int numDepths);

// Generate a circuit with random state vector using RY, CX gates
QCircuit RandomRandom(int numQubits, int numDepths);

// SupermarQ
QCircuit GHZ(int numQubits);
// Mermin-Bell only 3q and 4q
QCircuit BitCode(int numQubits, int numRounds=2);
QCircuit PhaseCode(int numQubits, int numRounds=2);
QCircuit VQE(int numQubits, int numLayers=2);
QCircuit Hamiltonian(int numQubits, int numSteps=2);
QCircuit QAOAFermionicSwapProxy(int numQubits);
