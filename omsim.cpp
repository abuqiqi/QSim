#include "omsim.h"

Matrix<DTYPE> OMSim(Matrix<DTYPE>& sv, QCircuit& qc) {
    Matrix<DTYPE> opmat = getOperationMatrix(qc);
    sv = opmat * sv;
    return opmat;
}

/**
 * @brief Construct the operation matrix of a quantum circuit
 * 
 * @param qc a quantum circuit
 * @return Matrix<DTYPE> the operation matrix
 */
Matrix<DTYPE> getOperationMatrix(QCircuit& qc) {
    // qc.print();

    Matrix<DTYPE> opmat, levelmat;
    opmat.identity(1 << qc.numQubits);
    levelmat.identity(2);

    // calculate the operation matrix of the quantum circuit
    for (int j = 0; j < qc.numDepths; ++ j) {
        // cout << "[DEBUG] [OMSim] level " << j << endl;
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
        // Calculate the operation matrix of level j //////////////////////////
        // Step 1. Let levelmat be the complete gate matrix of the highest gate
        levelmat = move(getCompleteMatrix(qc.gates[j][qid]));
        // levelmat.print();
        // ///////////////////////////////////////////////////////////////////////////
        // Step 2. Get the complete gate matrices of the remaining gates
        //      Step 2.1. Skip the MARK gates
        //      Step 2.2. Calculate the tensor product of the gate matrices
        for (int i = qid - 1; i >= 0; -- i) {
            if (qc.gates[j][i].isMARK()) {
                continue;
            }
            Matrix<DTYPE> tmpmat = move(getCompleteMatrix(qc.gates[j][i]));
            // qc.gates[j][i].print();
            // cout << "[DEBUG] tmpmat: " << endl;
            // tmpmat.print();
            levelmat = move(levelmat.tensorProduct(tmpmat));
        }
        // levelmat.print();
        // cout << endl;
        // ///////////////////////////////////////////////////////////////////////////
        // Step 3. Update the operation matrix opmat for the entire circuit
        opmat = levelmat * opmat;
        // ///////////////////////////////////////////////////////////////////////////
    }
    // cout << endl;
    return opmat;
}

//
// Utility functions
//

/**
 * @brief Get a complete gate matrix according to the applied qubits
 * 
 * @param gate the processing gate
 * @return Matrix<DTYPE> a complete gate matrix
 */
Matrix<DTYPE> getCompleteMatrix(QGate& gate) {
    if (gate.isMARK() || gate.isIDE() /*|| gate.isSingle()*/) {
        // return * gate.gmat;
        return * Matrix<DTYPE>::MatrixDict[gate.gname];
    }
    if (gate.is2QubitControlled()) {
        // Return the complete matrix of a 2-qubit controlled gate
        return genControlledGateMatrix(gate);
    }
    if (gate.gname == "SWAP") {
        // Return the complete matrix of a SWAP gate
        return genSwapGateMatrix(gate);
    }
    if (gate.is2QubitNonControlled() && (gate.targetQubits[0] + 1 != gate.targetQubits[1])) { // non-adjacent
        return gen2QubitGateMatrix(gate);
    }
    return * gate.gmat;
    // cout << "[ERROR] getCompleteMatrix: " << gate.gname << " not implemented" << endl;
    // exit(1);
}

/**
 * @brief Generate the gate matrix of a controlled gate
 *
 * @param gate the processing gate
 * @return Matrix<DTYPE> a complete gate matrix
 */
Matrix<DTYPE> genControlledGateMatrix(QGate& gate) {
    int ctrl = gate.controlQubits[0];
    int targ = gate.targetQubits[0];
    Matrix<DTYPE> ctrlmat, basismat, IDE;
    ctrlmat.zero(1 << (abs(ctrl - targ) + 1), 1 << (abs(ctrl - targ) + 1)); // initialize the complete matrix with all zeros
    IDE.identity(2);
    ll mask = ctrl > targ ? (1 << (ctrl - targ - 1)) : 1; // mask the control qubit

    for (ll i = 0; i < (1 << abs(ctrl-targ)); ++ i) {
        // basismat = | i >< i |
        basismat.zero(1 << abs(ctrl-targ), 1 << abs(ctrl-targ));
        basismat.data[i][i] = 1;

        // Calculate the complete gate matrix of a 2-qubit controlled gate
        //  Case 1. If ctrl = 1 and ctrl > targ, ctrlmat += | i >< i | \otimes gate
        //  Case 2. If ctrl = 1 and ctrl < targ, ctrlmat += gate \otimes | i >< i |
        //  Case 3. If ctrl = 0 and ctrl > targ, ctrlmat += | i >< i | \otimes IDE
        //  Case 4. If ctrl = 0 and ctrl < targ, ctrlmat += IDE \otimes | i >< i |
        if ((i & mask) == mask) { // control qubit = 1
            if (ctrl > targ) { // ctrlmat += | i >< i | \otimes gate
                ctrlmat += basismat.tensorProduct(*gate.gmat);
            } else { // ctrlmat += gate \otimes | i >< i |
                ctrlmat += gate.gmat->tensorProduct(basismat);
            }
        } else {
            if (ctrl > targ) { // ctrlmat += | i >< i | \otimes IDE
                ctrlmat += basismat.tensorProduct(IDE);
            } else { // ctrlmat += IDE \otimes | i >< i |
                ctrlmat += IDE.tensorProduct(basismat);
            }
        }
    }
    return ctrlmat;
}

/**
 * @brief Generate the gate matrix of a SWAP gate
 * 
 * @param gate the processing SWAP gate
 * @return Matrix<DTYPE> a complete gate matrix
 */
Matrix<DTYPE> genSwapGateMatrix(QGate& gate) {
    // when adding a SWAP, the target qubits are sorted in ascending order
    int span = gate.targetQubits[1] - gate.targetQubits[0] + 1;

    Matrix<DTYPE> mat;
    mat.identity(1 << span);

    ll mask0 = (1 << (span - 1));
    ll mask1 = 1;
    ll row;

    for (ll i = 0; i < (1 << span); ++ i) {
        if ((i & mask0) == 0 && (i & mask1) == mask1) {
            // i   := |0..1>
            // row := |1..0>
            row = i ^ mask0 ^ mask1;
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

// Swap two columns of a gate matrix
void swapCol(ll c1, ll c2, Matrix<DTYPE>& gate) {
    DTYPE tmp[gate.col];
    for (ll i = 0; i < gate.row; ++ i) {
        tmp[i] = gate.data[i][c1];
        gate.data[i][c1] = gate.data[i][c2];
        gate.data[i][c2] = tmp[i];
    }
}

Matrix<DTYPE> gen2QubitGateMatrix(QGate& gate) {
    int span = gate.targetQubits[1] - gate.targetQubits[0] + 1;
    ll n = (1 << span);
    ll block_size = n / 4;
    Matrix<DTYPE> mat;
    mat.identity(n);

    // MAT \otimes I \otimes I ...
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            ll base_row = i * block_size; // the left-upper element index
            ll base_col = j * block_size;
            for (ll k = 0; k < block_size; ++k) { // diagonal
                mat.data[base_row + k][base_col + k] = gate.gmat->data[i][j];
            }
        }
    }

    // swap rows and cols
    ll mask0 = (1 << (span - 2));
    ll mask1 = 1;
    ll row;
    for (ll i = 0; i < (1 << span); ++ i) {
        if ((i & mask0) == 0 && (i & mask1) == mask1) {
            // i   := |.0..1>
            // row := |.1..0>
            row = i ^ mask0 ^ mask1;
            swapRow(i, row, mat);
            swapCol(i, row, mat);
        }
    }
    return mat;
}
