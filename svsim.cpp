#include "svsim.h"

/**
 * @brief State vector simulation of a quantum circuit on the state vector
 * 
 * @param sv    the initial state vector
 * @param qc    a quantum circuit
 */
void SVSim(Matrix<DTYPE>& sv, QCircuit& qc) {
    for (int j = 0; j < qc.numDepths; ++ j) {
        for (int qid = 0; qid < qc.numQubits; ++ qid) {
            QGate& gate = qc.gates[j][qid];
            if (gate.isIDE() || gate.isMARK()) {
                continue;
            }
            svsimForGate(sv, gate);
        }
    }
}

/**
 * @brief State vector simulation for a quantum gate
 * 
 * @param sv    the state vector
 * @param gate  the processing gate
 */
void svsimForGate(Matrix<DTYPE>& sv, QGate& gate) {
    bool isAccessed[sv.row];
    memset(isAccessed, 0, sv.row*sizeof(bool));
    
    ll numAmps = (1 << gate.numTargets()); // the number of amplitudes involved in matrix-vector multiplication
    Matrix<DTYPE> amps_vec(numAmps, 1);

    // 1. Calculate the strides for the involved amplitudes
    /*  <e.g.>
        (1) If there is only one target qubit, two amplitudes are involved. 
            If the target qubit is q[k], strides = {0, 2^k}. 
        (2) If there are two target qubits, four amplitudes are involved. 
                      c  t
            idx[0] |..0..0..> -+ strides[0] = 0 --+------------------+
                               |                  |                  |
            idx[1] |..0..1..> -+ strides[1] = 2^t |                  |
                                                  |                  |
            idx[2] |..1..0..> --------------------+ strides[2] = 2^c |
                                                                     |
            idx[3] |..1..1..> ---------------------------------------+ strides[3] = 2^c + 2^t
        (3) If there are 'gate.numTargets()' target qubits, there are 'numAmps' amplitudes. 
    */
    vector<ll> strides;
    // [TODO] //////////////////////////////////////////////////
    for (ll idx = 0; idx < numAmps; ++ idx) {
        ll stride = 0;
        for (int j = 0; j < gate.numTargets(); ++ j) {
            if (idx & (1 << j)) { // if the j-th bit of idx is 1
                stride += (1 << gate.targetQubits[j]);
            }
        }
        strides.push_back(stride);
    }
    ////////////////////////////////////////////////////////////

    // 2. Iterate over all amplitudes
    for (ll ampidx = 0; ampidx < sv.row; ++ ampidx) {
        // 2.1. Skip the amplitude if it is already accessed
        if (isAccessed[ampidx]) continue;

        // 2.2. Get the involved amplitudes and mark them as accessed
        for (ll idx = 0; idx < numAmps; ++ idx) {
            if (ampidx + strides[idx] >= sv.row) {
                cout << "[ERROR] Exceed the length of the state vector." << endl;
                exit(1);
            }
            amps_vec.data[idx][0] = sv.data[ampidx + strides[idx]][0];
            isAccessed[ampidx + strides[idx]] = true;
        }

        // 3. Check the control bits of the current amplitude
        //    If the control bits are not satisfied, skip this amplitude group
        int ret = checkControlMask(ampidx, gate);
        if (ret == 0) {
            continue;
        }

        // 4. [TODO] Apply the gate matrix
        amps_vec = (*gate.gmat) * amps_vec;
        for (ll idx = 0; idx < numAmps; ++ idx) {
            sv.data[ampidx + strides[idx]][0] = amps_vec.data[idx][0];
        }
    }
}

/**
 * @brief Check if the index of a amplitude contains a legal control mask of the gate
 * 
 * @param amp  the amplitude index
 * @param gate the processing gate
 * @return int 0: illegal control mask; 1: legal control mask
 */
int checkControlMask(ll amp, QGate& gate) {
    int ctrl;
    ll ctrlmask = 0;
    for (int i = 0; i < gate.numControls(); ++ i) {
        // [TODO] Check the control qubit of the gate /////////////
        ctrl = gate.controlQubits[i];
        ctrlmask += (1 << ctrl);
        
        // 0-controlled and the control qubit of amp is 1
        if (ctrl < 0 && (amp & ctrlmask) == 1) {
            return 0;
        }

        // 1-controlled and the control qubit of amp is 0
        if (ctrl > 0 && (amp & ctrlmask) == 0) {
            return 0;
        }
        ////////////////////////////////////////////////////////////
    }
    return 1;
}