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
void QCircuit::applyGates(const vector<QGate>& gateSeq) {
    for (auto& gate : gateSeq) {
        if (gate.gname == "H") {
            h(gate.targetQubits[0]);
        } else if (gate.gname == "X") {
            x(gate.targetQubits[0]);
        } else if (gate.gname == "Y") {
            y(gate.targetQubits[0]);
        } else if (gate.gname == "Z") {
            z(gate.targetQubits[0]);
        } else if (gate.gname == "RX") {
            rx(gate.params[0], gate.targetQubits[0]);
        } else if (gate.gname == "RY") {
            ry(gate.params[0], gate.targetQubits[0]);
        } else if (gate.gname == "RZ") {
            rz(gate.params[0], gate.targetQubits[0]);
        } else if (gate.gname == "U1") {
            u1(gate.params[0], gate.targetQubits[0]);
        } else if (gate.gname == "CX") {
            cx(gate.controlQubits[0], gate.targetQubits[0]);
        } else if (gate.gname == "CY") {
            cy(gate.controlQubits[0], gate.targetQubits[0]);
        } else if (gate.gname == "CZ") {
            cz(gate.controlQubits[0], gate.targetQubits[0]);
        } else if (gate.gname == "CS") {
            cs(gate.controlQubits[0], gate.targetQubits[0]);
        } else if (gate.gname == "CT") {
            ct(gate.controlQubits[0], gate.targetQubits[0]);
        } else if (gate.gname == "CRZ") {
            crz(gate.params[0], gate.controlQubits[0], gate.targetQubits[0]);
        } else if (gate.gname == "CU1") {
            cu1(gate.params[0], gate.controlQubits[0], gate.targetQubits[0]);
        } else if (gate.gname == "SWAP") {
            swap(gate.targetQubits[0], gate.targetQubits[1]);
        } else if (gate.gname == "CCZ") {
            ccz(gate.controlQubits, gate.targetQubits[0]);
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
    cout << "[INFO] [" << name << "] numQubits: [" << numQubits << "] numDepths: [" << numDepths << "]" << endl;
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
