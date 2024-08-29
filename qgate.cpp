#include "qgate.h"

QGate::QGate() {
    gname = "NULL";
    controlQubits = {};
    targetQubits = {};
    params = {};
    gmat = nullptr;
}

/**
 * @brief Construct a new QGate::QGate object, initialize the gate matrix with the given name
 * 
 * @param gname_ the gate name
 * @param controls_ control qubits
 * @param targets_ target qubits
 */
QGate::QGate(string gname_, vector<int> controls_, vector<int> targets_) {
    gname = gname_;
    controlQubits = controls_;
    targetQubits = targets_;
    params = {};
    gmat = Matrix<DTYPE>::MatrixDict[gname];
    if (gmat == nullptr) {
        cout << "[ERROR] Gate " << gname << " not found in MatrixDict" << endl;
        exit(1);
    }
}

/**
 * @brief Construct a new QGate::QGate object with a parameter
 * 
 * @param gname_ the gate name
 * @param controls_ control qubits
 * @param targets_ target qubits
 * @param theta a parameter
 */
QGate::QGate(string gname_, vector<int> controls_, vector<int> targets_, double theta) {
    gname = gname_;
    controlQubits = controls_;
    targetQubits = targets_;
    params = {theta};

    string matkey = gname + to_string(theta);
    gmat = Matrix<DTYPE>::MatrixDict[matkey];
    if (gmat != nullptr) { // the gate matrix already exists
        // cout << "[DEBUG] Matrix already exists: " << matkey << ", " << gmat << endl;
        return;
    }

    Matrix<DTYPE> mat;
    if (gname == "RX") {
        mat.rotationX(theta);
    } else if (gname == "RY") {
        mat.rotationY(theta);
    } else if (gname == "RZ") {
        mat.rotationZ(theta);
    } else if (gname == "U1") {
        mat.u1(theta);
    } else {
        cout << "[ERROR] Gate " << gname << " not implemented" << endl;
        exit(1);
    }
    Matrix<DTYPE>::MatrixDict[matkey] = make_shared<Matrix<DTYPE>>(move(mat));
    gmat = Matrix<DTYPE>::MatrixDict[matkey];
}

/**
 * @brief Construct a new QGate::QGate object with a given matrix
 * 
 * @param gname_ the gate name
 * @param controls_ control qubits
 * @param targets_ target qubits
 * @param mat the gate matrix
 */
QGate::QGate(string gname_, vector<int> controls_, vector<int> targets_, Matrix<DTYPE>& mat) {
    gname = gname_;
    controlQubits = controls_;
    targetQubits = targets_;
    params = {};
    gmat = Matrix<DTYPE>::MatrixDict[gname];
    if (gmat == nullptr) {
        Matrix<DTYPE>::MatrixDict[gname] = make_shared<Matrix<DTYPE>>(move(mat));
        gmat = Matrix<DTYPE>::MatrixDict[gname];
    }
}

// Copy construct a new QGate::QGate object
QGate::QGate(const QGate& other) {
    gname = other.gname;
    controlQubits = other.controlQubits;
    targetQubits = other.targetQubits;
    params = other.params;
    gmat = other.gmat;
}

// Copy assignment
QGate& QGate::operator=(const QGate& other) {
    gname = other.gname;
    controlQubits = other.controlQubits;
    targetQubits = other.targetQubits;
    params = other.params;
    gmat = other.gmat;
    return *this;
}

// Compare two gates
bool QGate::operator==(const QGate& other) const {
    return gname == other.gname && controlQubits == other.controlQubits && targetQubits == other.targetQubits && params == other.params;
}

// Return the qubits of the gate
vector<int> QGate::qubits() {
    vector<int> qubits = controlQubits;
    qubits.insert(qubits.end(), targetQubits.begin(), targetQubits.end());
    sort(qubits.begin(), qubits.end());
    return qubits;
}

// Return the number of input/output qubits of the gate
int QGate::numQubits() {
    return controlQubits.size() + targetQubits.size();
}

// Return the number of control qubits of the gate
int QGate::numControls() {
    return controlQubits.size();
}

// Return the number of target qubits of the gate
int QGate::numTargets() {
    return targetQubits.size();
}

// Return the key of the gate matrix in the MatrixDict
string QGate::gmatKey() {
    string matkey = gname;
    for (const auto& param : params) {
        matkey += to_string(param);
    }
    return matkey;
}

// Get the full gate matrix for a controlled gate
shared_ptr<Matrix<DTYPE>> QGate::getFullMatrix() {
    if (is2QubitControlled()) { // generate a 2-qubit controlled gate's matrix
        cout << "[DEBUG] is2QubitControlled" << endl;
        string matkey = gmatKey();
        shared_ptr<Matrix<DTYPE>> fullmat;
        if (controlQubits[0] > targetQubits[0]) { // CU gate
            matkey = "C" + matkey;
            fullmat = Matrix<DTYPE>::MatrixDict[matkey];
            if (fullmat == nullptr) { // insert a new entry to the MatrixDict
                Matrix<DTYPE> mat(4, 4);
                Matrix<DTYPE> basismat(2, 2);
                basismat.data[0][0] = 1;
                mat = basismat.tensorProduct(*Matrix<DTYPE>::MatrixDict["IDE"]);
                basismat.data[0][0] = 0;
                basismat.data[1][1] = 1;
                mat += basismat.tensorProduct(*gmat);
                Matrix<DTYPE>::MatrixDict[matkey] = make_shared<Matrix<DTYPE>>(move(mat));
                fullmat = Matrix<DTYPE>::MatrixDict[matkey];
            }
        } else { // UC gate
            matkey = matkey + "C";
            fullmat = Matrix<DTYPE>::MatrixDict[matkey];
            if (fullmat == nullptr) { // insert a new entry to the MatrixDict
                Matrix<DTYPE> mat(4, 4);
                Matrix<DTYPE> basismat(2, 2);
                basismat.data[0][0] = 1;
                mat = Matrix<DTYPE>::MatrixDict["IDE"]->tensorProduct(basismat);
                basismat.data[0][0] = 0;
                basismat.data[1][1] = 1;
                mat += gmat->tensorProduct(basismat);
                Matrix<DTYPE>::MatrixDict[matkey] = make_shared<Matrix<DTYPE>>(move(mat));
                fullmat = Matrix<DTYPE>::MatrixDict[matkey];
            }
        }
        return fullmat;
    } else if (isControlled()) { // generate a controlled gate's matrix
        cout << "[TODO] Generate the full matrix of a controlled gate for state vector simulation." << endl;
        exit(1);
    }
    if (gmat == nullptr) {
        cout << "[ERROR] Gate " << gname << " not found in MatrixDict" << endl;
        exit(1);
    }
    return gmat;
}

// Check if the gate is an identity gate
bool QGate::isIDE() {
    return gname == "IDE";
}

// Check if the gate is a placeholder gate
bool QGate::isMARK() {
    return gname == "MARK";
}

// Check if the gate is a single-qubit gate
bool QGate::isSingle() {
    return gname != "IDE" && gname != "MARK" && controlQubits.size() == 0 && targetQubits.size() == 1;
}

// Check if the gate is a controlled gate
bool QGate::isControlled() {
    return gname != "MARK" && controlQubits.size() > 0;
}

// Check if the gate is a 2-qubit controlled gate
bool QGate::is2QubitControlled() {
    return gname != "MARK" && controlQubits.size() == 1 && targetQubits.size() == 1;
}

// Check if the gate is hermitian
bool QGate::isHermitian() {
    return gname == "X" || gname == "Y" || gname == "Z" || gname == "H" || gname == "SWAP";
}

// Check if the gate is a phase gate
bool QGate::isPhase() {
    return gname == "Z" || gname == "S" || gname == "T" || gname == "U1";
}

// Check if qubit[qid] is a control qubit of the gate
bool QGate::isControlQubit(int qid) {
    return find(controlQubits.begin(), controlQubits.end(), qid) != controlQubits.end();
}

// Check if qubit[qid] is a target qubit of the gate
bool QGate::isTargetQubit(int qid) {
    return gname != "IDE" && gname != "MARK" && find(targetQubits.begin(), targetQubits.end(), qid) != targetQubits.end();
}

// Print the gate information
void QGate::printInfo() {
    cout << "===== Gate: " << gname << " =====" << endl;
    cout << "Control qubits: ";
    for (const auto& ctrl : controlQubits) {
        cout << ctrl << " ";
    }
    cout << endl;
    cout << "Target qubits: ";
    for (const auto& targ : targetQubits) {
        cout << targ << " ";
    }
    cout << endl;
    cout << "Parameters: ";
    for (const auto& param : params) {
        cout << param << " ";
    }
    cout << endl;
    cout << "=====================" << endl;
}

// Print the gate information and gate matrix
void QGate::print() {
    printInfo();
    if (gmat != nullptr)
        gmat->print();
    cout << "=====================" << endl;
}

// Destructor
QGate::~QGate() {
    return;
}

// Compare two integers by their absolute values
// Control qubits can be negative to denote 0-controlled
bool compareByAbsoluteValue(int a, int b) {
    return std::abs(a) < std::abs(b);
}