#include "qcircuit.h"

QCircuit::QCircuit() {}

/**
 * @brief Construct an n-qubit 1-level quantum circuit object
 * 
 * @param numQubits_ #Qubits
 * @param numDepths_ #Depths
 */
QCircuit::QCircuit(int numQubits_, string name_){
    numQubits = numQubits_;
    numDepths = 0;
    add_level(); // numDepths += 1
    name = name_;
}

QCircuit::QCircuit(int numQubits_, int numDepths_, string name_){
    numQubits = numQubits_;
    while (numDepths < numDepths_) {
        add_level(); // numDepths += 1
    }
    name = name_;
}

QCircuit::QCircuit(QCircuit& qc, int qidFrom, int qidTo, int depthFrom, int depthTo) {
    numQubits = qidTo - qidFrom;
    numDepths = depthTo - depthFrom;
    name = qc.name + "_q" + to_string(qidFrom) + "-" + to_string(qidTo) + "_d" + to_string(depthFrom) + "-" + to_string(depthTo);
    for (int j = depthFrom; j < depthTo; ++ j) {
        vector<QGate> level;
        for (int i = qidFrom; i < qidTo; ++ i) {
            QGate gate = qc.gates[j][i];
            if (gate.isMARK() && (gate.targetQubits[0] < qidFrom || gate.targetQubits[0] >= qidTo)) {
                gate.gname = "IDE";
                level.push_back(gate);
            } else {
                // remap qubits
                for (auto& qid : gate.controlQubits) {
                    qid -= qidFrom;
                }
                for (auto& qid : gate.targetQubits) {
                    qid -= qidFrom;
                }
                level.push_back(gate);
            }
        }
        gates.push_back(level);
    }
}

// 
// Single-qubit gates
// 

/**
 * @brief Apply an H gate to qubit[qid]
 * 
 * @param qid   qubit id
 */
void QCircuit::h(int qid) {
    if (! gates[numDepths-1][qid].isIDE()) {
        add_level();
    }
    gates[numDepths-1][qid] = QGate("H", {}, {qid});
}

/**
 * @brief Apply an X gate to qubit[qid]
 * 
 * @param qid   qubit id
 */
void QCircuit::x(int qid) {
    if (! gates[numDepths-1][qid].isIDE()) {
        add_level();
    }
    gates[numDepths-1][qid] = QGate("X", {}, {qid});
}

/**
 * @brief Apply a Y gate to qubit[qid]
 * 
 */
void QCircuit::y(int qid) {
    if (! gates[numDepths-1][qid].isIDE()) {
        add_level();
    }
    gates[numDepths-1][qid] = QGate("Y", {}, {qid});
}

/**
 * @brief Apply a Z gate to qubit[qid]
 * 
 * @param qid   qubit id
 */
void QCircuit::z(int qid) {
    if (! gates[numDepths-1][qid].isIDE()) {
        add_level();
    }
    gates[numDepths-1][qid] = QGate("Z", {}, {qid});
}

/**
 * @brief Apply an RX gate to qubit[qid]
 * 
 * @param theta the gate parameter
 * @param qid   qubit id
 */
void QCircuit::rx(double theta, int qid) {
    if (! gates[numDepths-1][qid].isIDE()) {
        add_level();
    }
    gates[numDepths-1][qid] = QGate("RX", {}, {qid}, theta);
}

/**
 * @brief Apply an RY gate to qubit[qid]
 * 
 * @param theta the gate parameter
 * @param qid   qubit id
 */
void QCircuit::ry(double theta, int qid) {
    if (! gates[numDepths-1][qid].isIDE()) {
        add_level();
    }
    gates[numDepths-1][qid] = QGate("RY", {}, {qid}, theta);
}

/**
 * @brief Apply an RZ gate to qubit[qid]
 * 
 * @param theta the gate parameter
 * @param qid   qubit id
 */
void QCircuit::rz(double theta, int qid) {
    if (! gates[numDepths-1][qid].isIDE()) {
        add_level();
    }
    gates[numDepths-1][qid] = QGate("RZ", {}, {qid}, theta);
}

/**
 * @brief Apply a U1 gate to qubit[qid]
 * 
 * @param lambda the gate parameter
 * @param qid    qubit id
 */
void QCircuit::u1(double lambda, int qid) {
    if (! gates[numDepths-1][qid].isIDE()) {
        add_level();
    }
    gates[numDepths-1][qid] = QGate("U1", {}, {qid}, lambda);
}

// 
// 2-qubit gates
// 

/**
 * @brief Apply a CX gate to qubit[ctrl] and qubit[targ]
 * 
 * @param ctrl  control qubit id
 * @param targ  target qubit id
 */
void QCircuit::cx(int ctrl, int targ) {
    int start = min(ctrl, targ);
    int end = max(ctrl, targ);
    for (int i = start; i <= end; ++ i) {
        if (! gates[numDepths-1][i].isIDE()) {
            add_level();
            break;
        }
    }
    for (int i = start; i <= end; ++ i) {
        gates[numDepths-1][i] = QGate("MARK", {ctrl}, {targ});
    }
    gates[numDepths-1][targ] = QGate("X", {ctrl}, {targ});
}

/**
 * @brief Apply a CY gate to qubit[ctrl] and qubit[targ]
 * 
 * @param ctrl  control qubit id
 * @param targ  target qubit id
 */
void QCircuit::cy(int ctrl, int targ) {
    int start = min(ctrl, targ);
    int end = max(ctrl, targ);
    for (int i = start; i <= end; ++ i) {
        if (! gates[numDepths-1][i].isIDE()) {
            add_level();
            break;
        }
    }
    for (int i = start; i <= end; ++i) {
        gates[numDepths-1][i] = QGate("MARK", {ctrl}, {targ});
    }
    gates[numDepths-1][targ] = QGate("Y", {ctrl}, {targ});
}

/**
 * @brief Apply a CZ gate to qubit[ctrl] and qubit[targ]
 * 
 * @param ctrl  control qubit id
 * @param targ  target qubit id
 */
void QCircuit::cz(int ctrl, int targ) {
    int start = min(ctrl, targ);
    int end = max(ctrl, targ);
    for (int i = start; i <= end; ++ i) {
        if (! gates[numDepths-1][i].isIDE()) {
            add_level();
            break;
        }
    }
    for (int i = start; i <= end; ++i) {
        gates[numDepths-1][i] = QGate("MARK", {ctrl}, {targ});
    }
    gates[numDepths-1][targ] = QGate("Z", {ctrl}, {targ});
}

void QCircuit::cs(int ctrl, int targ) {
    int start = min(ctrl, targ);
    int end = max(ctrl, targ);
    for (int i = start; i <= end; ++ i) {
        if (! gates[numDepths-1][i].isIDE()) {
            add_level();
            break;
        }
    }
    for (int i = start; i <= end; ++i) {
        gates[numDepths-1][i] = QGate("MARK", {ctrl}, {targ});
    }
    gates[numDepths-1][targ] = QGate("S", {ctrl}, {targ});
}

void QCircuit::ct(int ctrl, int targ) {
    int start = min(ctrl, targ);
    int end = max(ctrl, targ);
    for (int i = start; i <= end; ++ i) {
        if (! gates[numDepths-1][i].isIDE()) {
            add_level();
            break;
        }
    }
    for (int i = start; i <= end; ++i) {
        gates[numDepths-1][i] = QGate("MARK", {ctrl}, {targ});
    }
    gates[numDepths-1][targ] = QGate("T", {ctrl}, {targ});
}

void QCircuit::crz(double theta, int ctrl, int targ) {
    int start = min(ctrl, targ);
    int end = max(ctrl, targ);
    for (int i = start; i <= end; ++ i) {
        if (! gates[numDepths-1][i].isIDE()) {
            add_level();
            break;
        }
    }
    for (int i = start; i <= end; ++i) {
        gates[numDepths-1][i] = QGate("MARK", {ctrl}, {targ});
    }
    gates[numDepths-1][targ] = QGate("RZ", {ctrl}, {targ}, theta);
}

void QCircuit::cu1(double lambda, int ctrl, int targ) {
    int start = min(ctrl, targ);
    int end = max(ctrl, targ);
    for (int i = start; i <= end; ++ i) {
        if (! gates[numDepths-1][i].isIDE()) {
            add_level();
            break;
        }
    }
    for (int i = start; i <= end; ++i) {
        gates[numDepths-1][i] = QGate("MARK", {ctrl}, {targ});
    }
    gates[numDepths-1][targ] = QGate("U1", {ctrl}, {targ}, lambda);
}

/**
 * @brief Apply a SWAP gate to qubit[qid1] and qubit[qid2]
 * 
 * @param qid1  qubit id 1
 * @param qid2  qubit id 2
 */
void QCircuit::swap(int qid1, int qid2) {
    int start = min(qid1, qid2);
    int end = max(qid1, qid2);
    for (int i = start; i <= end; ++ i) {
        if (! gates[numDepths-1][i].isIDE()) {
            add_level();
            break;
        }
    }
    for (int i = start; i <= end; ++ i) {
        gates[numDepths-1][i] = QGate("MARK", {}, {start, end});
    }
    gates[numDepths-1][end] = QGate("SWAP", {}, {start, end});
}

/**
 * @brief Apply a controlled-Z gate with multiple control qubits
 * 
 * @param ctrls 
 * @param targs 
 */
void QCircuit::ccz(vector<int> ctrls, int targ) {
    sort(ctrls.begin(), ctrls.end());
    int start = ctrls[0] < targ ? ctrls[0] : targ;
    int end = ctrls[ctrls.size()-1] > targ ? ctrls[ctrls.size()-1] : targ;
    for (int i = start; i <= end; ++ i) {
        if (! gates[numDepths-1][i].isIDE()) {
            add_level();
            break;
        }
    }
    for (int i = start; i <= end; ++ i) {
        gates[numDepths-1][i] = QGate("MARK", ctrls, {targ});
    }
    gates[numDepths-1][targ] = QGate("Z", ctrls, {targ});
}

/**
 * @brief Add a barrier to the quantum circuit
 */
void QCircuit::barrier() {
    add_level();
}

/**
 * @brief Set the circuit depth to numDepths_
 * 
 * @param numDepths_  the target circuit depth
 */
void QCircuit::setDepths(int numDepths_) {
    int range = numDepths_ - numDepths;
    for (int i = 0; i < range; ++ i){
        add_level();
    }
}

/**
 * @brief Apply a sequence of quantum gates to the circuit
 * 
 * @param gateSeq  the sequence of quantum gates
 */
void QCircuit::applyGates(vector<QGate>& gateSeq) {
    for (auto& gate : gateSeq) {
        if (gate.isIDE()) continue;
        for (auto& qid : gate.qubits()) {
            if (! gates[numDepths-1][qid].isIDE()) {
                add_level();
                break;
            }
        }
        int start = gate.qubits()[0];
        int end = gate.qubits()[gate.qubits().size()-1];
        for (int i = start; i <= end; ++ i) {
            gates[numDepths-1][i] = QGate("MARK", gate.controlQubits, gate.targetQubits);
        }
        int targ = gate.targetQubits[gate.targetQubits.size()-1]; // the highest target qubit
        if (gate.params.size() > 0) {
            gates[numDepths-1][targ] = QGate(gate.gname, gate.controlQubits, gate.targetQubits, gate.params[0]);
        } else {
            gates[numDepths-1][targ] = QGate(gate.gname, gate.controlQubits, gate.targetQubits, *gate.gmat);
        }
    }
}

/**
 * @brief Print the structure of the quantum circuit
 */
void QCircuit::print() {
    printInfo();
    int start = 0;
    if (numQubits >= 6) {
        start = numQubits - 6;
    }
    for (int i = numQubits - 1; i >= start; -- i) {
        cout << "q[" << i << "]\t";
        for (int j = 0; j < numDepths; ++ j) {
            if (j > 10) {
                cout << "...";
                break;
            }
            if (gates[j][i].isControlQubit(i)) {
                cout << "C";
            } else if (gates[j][i].isTargetQubit(i)) {
                cout << "T";
            }
            cout << gates[j][i].gname << "\t"; 
        }
        cout << endl;
    }
    cout << endl;
}

/**
 * @brief Print the information of the quantum circuit
 */
void QCircuit::printInfo() {
    cout << "[INFO] [" << name << "] numQubits: [" << numQubits << "] numDepths: [" << numDepths << "] numGates: [" << numQGates() << "]" << endl;
}

/**
 * @brief Add a new level full of IDE gates to the circuit
 */
void QCircuit::add_level() {
    vector<QGate> level;
    for (int i = 0; i < numQubits; ++ i) {
        level.push_back(QGate("IDE", {}, {i}));
    }
    gates.push_back(level);
    numDepths ++;
}

// Clear all gates in the circuit
void QCircuit::clear() {
    gates.clear();
    numDepths = 0;
    add_level(); // numDepths += 1
}

/**
 * @brief Calculate the number of quantum gates in the circuit
 * 
 * @return ll 
 */
ll QCircuit::numQGates() {
    ll numGates = 0;
    for (int j = 0; j < numDepths; ++ j) {
        for (int i = 0; i < numQubits; ++ i) {
            if (! gates[j][i].isIDE() && ! gates[j][i].isMARK()) {
                numGates ++;
            }
        }
    }
    return numGates;
}

string QCircuit::cmatkey() {
    string key = "";
    for (int i = 0; i < numQubits; ++ i) {
        if (i > 0) {
            key += "|";
        }
        for (int j = 0; j < numDepths; ++ j) {
            key += gates[j][i].gmatKey();
        }
    }
    return key;
}
