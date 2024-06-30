#include "svsim.h"

/**
 * @brief State vector simulation of a gate on the state vector
 * 
 * @param sv    the initial state vector
 * @param svlen the length of the state vector
 * @param qc    a quantum circuit
 */
void StateVectorSim(Matrix<DTYPE>& sv, long long svlen, QCircuit& qc) {
    for (int j = 0; j < qc.numDepths; ++ j) {
        for (int qid = 0; qid < qc.numQubits; ++ qid) {
            QGate& gate = qc.gates[j][qid];
            if (gate.isIDE() || gate.is2QubitControl(qid)) {
                continue;
            }
            if (gate.is2x2GateMatrix()) {
                svsim2x2GateMatrix(sv, svlen, gate);
            } else if (gate.is4x4GateMatrix()) {
                svsim4x4GateMatrix(sv, svlen, gate);
            } else {
                cout << "[ERROR] Unimplemented gate: " << gate.gname << endl;
                exit(1);
            }
        }
    }
}

/**
 * @brief SVSim for 2x2 gates, including CU gates
 * 
 * @param sv    the state vector
 * @param svlen the length of the state vector
 * @param gate  the processing gate
 */
void svsim2x2GateMatrix(Matrix<DTYPE>& sv, long long svlen, QGate& gate) {
    bool isAccessed[svlen];
    memset(isAccessed, 0, svlen*sizeof(bool));

    int qid = gate.targetQubits[0];
    long long stride = (1 << qid); // stride = 2^qid for gate applied on q[qid]

    for (long long j = 0; j < svlen; ++ j) {
        if (isAccessed[j]) continue;

        if (j + stride >= svlen) {
            cout << "[ERROR] LocalComputing exceeds the length of state vector." << endl;
            exit(1);
        }

        isAccessed[j] = isAccessed[j + stride] = true;

        if (! gate.isSingle()) { // deal with controlled gates
            int ctrl = gate.controlQubits[0]; // the ctrl bit whose range is [0, n-1]
            long long ctrlmask = (1 << ctrl); // ctrl qubit = 1, else 0
            if ((j & ctrlmask) == 0) { // j is the binary index formed by low-order qubits
                continue; // the control qubit is 0, do nothing
            }
        }

        apply2x2GateMatrix(gate.gmat, sv.data[j][0], sv.data[j + stride][0]);
    }
}

/**
 * @brief Multiplication of a 2x2 gate matrix
 * 
 * @param gmat 2x2 gate matrix
 * @param a    amplitude 0
 * @param b    amplitude 1
 */
void apply2x2GateMatrix(shared_ptr<Matrix<DTYPE>> gmat, DTYPE &a, DTYPE &b) {
    DTYPE temp0 = a*gmat->data[0][0] + b*gmat->data[0][1];
    DTYPE temp1 = a*gmat->data[1][0] + b*gmat->data[1][1];
    a = temp0;
    b = temp1;
}

/**
 * @brief SVSim for 4x4 gates (SWAP)
 * 
 * @param sv    local state-vector
 * @param svlen length of local state-vector, i.e., 2^l
 * @param gate  the processing gate
 */
void svsim4x4GateMatrix(Matrix<DTYPE>& sv, long long svlen, QGate& gate) {
    int targ1 = gate.targetQubits[0];
    int targ2 = gate.targetQubits[1];

    bool isAccessed[svlen];
    memset(isAccessed, 0, svlen*sizeof(bool));

    /*  Since 2 input qubits are local, the amplitude starts from |0...0>
           c  t
        |..0..0..> -+------------------+
                    | strideLow = 2^t  |
        |..0..1..> -+                  |
                                       | strideHigh = 2^c
        |..1..0..> -+------------------+
                    | strideLow = 2^t
        |..1..1..> -+
    */
    long long strideLow = (1 << targ1);
    long long strideHigh = (1 << targ2);

    for (long long j = 0; j < svlen; ++ j) {
        if (isAccessed[j]) continue;

        if (j + strideHigh + strideLow >= svlen) {
            cout << "[ERROR] LocalComputing exceeds the length of state vector." << endl;
            exit(1);
        }

        isAccessed[j] = isAccessed[j + strideLow] = isAccessed[j + strideHigh] = isAccessed[j + strideHigh + strideLow] = true;

        apply4x4GateMatrix(gate.gmat, sv.data[j][0], sv.data[j + strideLow][0], sv.data[j + strideHigh][0], sv.data[j + strideHigh + strideLow][0]);
    }
}

/**
 * @brief Multiplication of a 4x4 gate matrix
 * 
 * @param gmat 4x4 gate matrix
 * @param a    amplitude 0
 * @param b    amplitude 1
 * @param c    amplitude 2
 * @param d    amplitude 3
 */
void apply4x4GateMatrix(shared_ptr<Matrix<DTYPE>> gmat, DTYPE &a, DTYPE &b, DTYPE &c, DTYPE &d) {
    DTYPE temp0 = a*gmat->data[0][0] + b*gmat->data[0][1] + c*gmat->data[0][2] + d*gmat->data[0][3];
    DTYPE temp1 = a*gmat->data[1][0] + b*gmat->data[1][1] + c*gmat->data[1][2] + d*gmat->data[1][3];
    DTYPE temp2 = a*gmat->data[2][0] + b*gmat->data[2][1] + c*gmat->data[2][2] + d*gmat->data[2][3];
    DTYPE temp3 = a*gmat->data[3][0] + b*gmat->data[3][1] + c*gmat->data[3][2] + d*gmat->data[3][3];
    a = temp0;
    b = temp1;
    c = temp2;
    d = temp3;
}