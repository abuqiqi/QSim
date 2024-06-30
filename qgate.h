#pragma once

#include "matrix.h"

class QGate {
public:
    string gname; // gate name
    vector<int> controlQubits; // the control qubits of the gate
    vector<int> targetQubits; // the target qubits of the gate
    shared_ptr<Matrix<DTYPE>> gmat; // the gate matrix
    
    QGate(string gname_, vector<int> controls_, vector<int> targets_);
    QGate(string gname_, vector<int> controls_, vector<int> targets_, double theta);
    QGate(const QGate& other);

    QGate& operator=(const QGate& other);

    int numQubits(); // the number of input/output qubits of the gate
    int numControlQubits(); // the number of control qubits of the gate
    int numTargetQubits(); // the number of target qubits of the gate

    bool isIDE(); // check if the gate is an identity gate
    bool isSingle(); // check if the gate is a single-qubit gate
    bool is2QubitControl(int qid); // check if qubit[qid] is the control qubit of the 2-qubit gate
    bool is2QubitTarget(int qid); // check if qubit[qid] is the target qubit of the 2-qubit gate
    bool isControlQubit(int qid); // check if qubit[qid] is a control qubit of the gate
    bool isTargetQubit(int qid); // check if qubit[qid] is a target qubit of the gate

    bool is2x2GateMatrix(); // check if the gate matrix is a 2x2 matrix
    bool is4x4GateMatrix(); // check if the gate matrix is a 4x4 matrix
    
    void print(); // print the gate information

    ~QGate();
};