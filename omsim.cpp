#include "omsim.h"

/**
 * @brief Operation matrix simulation of a quantum circuit
 * 
 * @param sv the state vector
 * @param qc a quantum circuit
 * @return Matrix<DTYPE> the operation matrix
 */
Matrix<DTYPE> OMSim(Matrix<DTYPE>& sv, QCircuit& qc) {
    Matrix<DTYPE> opmat, levelmat;
    opmat.identity(sv.row);
    levelmat.identity(2);

    // calculate the operation matrix of the quantum circuit
    for (int j = 0; j < qc.numDepths; ++ j) {
        int qid = qc.numQubits-1;

        // get the highest gate matrix
        while (qc.gates[j][qid].isMARK()) {
            // skip the pseudo placeholder MARK gates placed at control positions
            -- qid;
        }
        if (qid < 0) {
            cout << "[ERROR] Invalid level with no target gate: " << j << endl;
            exit(1);
        }
        // [TODO] Calculate the operation matrix of level j //////////////////////////
        // cout << "[TODO] Calculate the operation matrix of level j" << endl;
        // exit(1);
        // [TODO] Step 1. Let levelmat be the complete gate matrix of the highest gate
        levelmat = move(getCompleteMatrix(qc.gates[j][qid]));
        // [TODO] Step 2. Get the complete gate matrices of the remaining gates
        //        Step 2.1. Skip the Mark gates
        //        Step 2.2. Calculate the tensor product of the gate matrices
        for (int i = qid - 1; i >= 0; -- i) {
            if (qc.gates[j][i].isMARK()) {
                continue;
            }
            Matrix<DTYPE> tmpmat = move(getCompleteMatrix(qc.gates[j][i]));
            levelmat = move(levelmat.tensorProduct(tmpmat));
        }
        // [TODO] Step 3. Update the operation matrix opmat
        opmat = levelmat * opmat;
    }
    // [TODO] Update the state vector sv
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
    if (gate.isIDE() || gate.isSingle()) {
        return * gate.gmat;
    }
    if (gate.is2QubitControlled()) {
        // [TODO] Return the complete matrix of a 2-qubit controlled gate
        // cout << "[TODO] return the complete matrix of a 2-qubit controlled gate" << endl;
        // exit(1);
        return genControlledGateMatrix(gate);
    }
    if (gate.gname == "SWAP") {
        // [TODO] Return the complete matrix of a SWAP gate
        // cout << "[TODO] return the complete matrix of a SWAP gate" << endl;
        // exit(1);
        return genSwapGateMatrix(gate);
    }
    cout << "[ERROR] getCompleteMatrix: " << gate.gname << " not implemented" << endl;
    exit(1);
}

/**
 * @brief Generate the gate matrix of a 2-qubit controlled gate
 *
 * @param gate the processing gate
 * @return Matrix<DTYPE> a complete gate matrix
*/
Matrix<DTYPE> genControlledGateMatrix(QGate& gate) {
    int ctrl = gate.controlQubits[0];
    int targ = gate.targetQubits[0];
    Matrix<DTYPE> ctrlmat, basismat, IDE;
    ctrlmat.zero(1 << (abs(ctrl - targ) + 1)); // initialize the complete matrix with all zeros
    IDE.identity(2);
    ll mask = ((1 << abs(ctrl - targ)) - 1); // mask the control qubit from qubit[targ+-1] to qubit[ctrl]

    for (ll i = 0; i < (1 << abs(ctrl-targ)); ++ i) {
        // basis |i> has abs(ctrl-targ) qubits and length (1 << abs(ctrl-targ))
        basismat.zero(1 << abs(ctrl-targ));
        basismat.data[i][i] = 1;  // basismat = | i >< i |

        if ((i & mask) == mask) { // control qubit = 1
            if (ctrl > targ) { // ctrlmat += | i >< i | \otimes gate
                ctrlmat += basismat.tensorProduct(*gate.gmat);
            } else { // ctrlmat += gate \otimes | i >< i |
                ctrlmat += gate.gmat->tensorProduct(basismat);
            }
        } else {
            if (ctrl > targ) { // ctrlmat += | i >< i | \otimes I
                ctrlmat += basismat.tensorProduct(IDE);
            } else { // ctrlmat += I \otimes | i >< i |
                ctrlmat += IDE.tensorProduct(basismat);
            }
        }
    }
    return ctrlmat;
}

/**
 * @brief Generate the gate matrix of a swap gate
 * 
 * @param gate the processing SWAP gate
 * @return Matrix<DTYPE> a complete gate matrix
 */
Matrix<DTYPE> genSwapGateMatrix(QGate& gate) {
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
