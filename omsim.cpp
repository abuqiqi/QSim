#include "omsim.h"

/**
 * @brief Operation matrix simulation of a quantum circuit
 * 
 * @param sv the state vector
 * @param qc a quantum circuit
 * @return Matrix<DTYPE> the operation matrix
 */
Matrix<DTYPE> OperationMatrixSim(Matrix<DTYPE>& sv, QCircuit& qc) {
    Matrix<DTYPE> opmat, levelmat;
    opmat.identity(sv.row);
    levelmat.identity(2);

    // calculate the operation matrix of the quantum circuit
    for (int j = 0; j < qc.numDepths; ++ j) {

        int qid = qc.numQubits-1;

        // get the highest gate matrix
        while (qc.gates[j][qid].isMARK()) {
            -- qid;
        }
        if (qid < 0) {
            cout << "[ERROR] Invalid level with no target gate: " << j << endl;
            exit(1);
        }
        levelmat = move(getCompleteMatrix(qc.gates[j][qid]));

        for (int i = qid - 1; i >= 0; -- i) {
            if (qc.gates[j][i].isMARK()) {
                continue;
            }
            Matrix<DTYPE> tmpmat = move(getCompleteMatrix(qc.gates[j][i]));
            levelmat = move(levelmat.tensorProduct(tmpmat));
        }

        opmat = levelmat * opmat;
    }

    sv = opmat * sv;
    return opmat;
}

/**
 * @brief Get a complete gate matrix according to the applied qubits
 * 
 * @param gate the processing gate
 * @return Matrix<DTYPE> a complete gate matrix
 */
Matrix<DTYPE> getCompleteMatrix(QGate& gate) {
    if (gate.isIDE() || gate.isMARK() || gate.isSingle()) {
        return * gate.gmat;
    }
    if (gate.gname == "SWAP") {
        return genSwapGate(gate);
    }
    if (gate.numTargets() == 1) {
        return genControlGate(gate);
    }
    cout << "[ERROR] getCompleteMatrix: " << gate.gname << " not implemented" << endl;
    exit(1);
}

/**
 * @brief Generate the gate matrix of a swap gate
 * 
 * @param gate the processing SWAP gate
 * @return Matrix<DTYPE> a complete gate matrix
 */
Matrix<DTYPE> genSwapGate(QGate& gate) {
    int span = abs(gate.targetQubits[0] - gate.targetQubits[1]);
    Matrix<DTYPE> mat;
    mat.identity(1 << span);

    ll mask1 = (1 << gate.targetQubits[0]);
    ll mask2 = (1 << gate.targetQubits[1]);
    ll row;

    for (ll i = 0; i < (1 << span); ++ i) {
        if ((i & mask1) == 0 && (i & mask2) == mask2) {
            // suppose qid1 > qid2
            // i   := |..0..1..>
            // row := |..1..0..>
            row = i ^ mask1 ^ mask2;
            swapRow(i, row, mat);
        }
    }
    return mat;
}

/**
 * @brief Swap two rows of a gate matrix
 * 
 * @param r1   row index 1
 * @param r2   row index 2
 * @param gate return value
 */
void swapRow(ll r1, ll r2, Matrix<DTYPE>& gate) {
    DTYPE tmp[gate.row];
    memcpy(tmp, gate.data[r1], gate.row * sizeof(DTYPE));
    memcpy(gate.data[r1], gate.data[r2], gate.row * sizeof(DTYPE));
    memcpy(gate.data[r2], tmp, gate.row * sizeof(DTYPE));
}

/**
 * @brief Generate the gate matrix of a control gate
 *
 * @param gate the processing gate
 * @return Matrix<DTYPE> a complete gate matrix
*/
Matrix<DTYPE> genControlGate(QGate& gate) {
    cout << "[ERROR] Not implemented" << endl;
    return Matrix<DTYPE>();
}