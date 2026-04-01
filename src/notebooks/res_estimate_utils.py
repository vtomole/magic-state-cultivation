from qiskit import transpile
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import HGate, CXGate, SGate, CYGate
from qiskit.circuit.equivalence_library import SessionEquivalenceLibrary as sel

from mqt.core import load
import numpy as np

import stim

def calculate_movements_for_arch(arch, qc, compiler):
    circ = load(qc)
    code = compiler.compile(circ)
    num_moves = code.count("@+ move")
    return num_moves

from qiskit import QuantumCircuit

def remove_mid_circ_meas(qc: QuantumCircuit, return_meas_reset_moves: bool = False):
    num_qubits = qc.num_qubits
    moves = 0

    new_qc = QuantumCircuit(num_qubits)

    for instruction in qc:
        name = instruction.operation.name
        targets = [qc.find_bit(q).index for q in instruction.qubits]

        if name == "reset":
            moves += 1
        elif name == "measure":
            moves += 1
        elif name == "barrier":
            print("skipping barrier")
            # new_qc.barrier()
        else:
            new_qc.append(instruction.operation, targets)

    if return_meas_reset_moves:
        return new_qc, moves
    return new_qc



def stim_to_qiskit(stim_circuit: stim.Circuit, return_meas_reset_moves:bool = False) -> QuantumCircuit:
    num_qubits = stim_circuit.num_qubits
    num_measures = stim_circuit.num_measurements
    qc = QuantumCircuit(num_qubits, num_measures)
    meas_idx = 0
    moves = 0
    for instruction in stim_circuit:
        name = instruction.name
        targets = instruction.targets_copy()
        
        if name == "H":
            for t in targets:
                qc.h(t.value)
        elif name == "X":
            for t in targets:
                qc.x(t.value)
        elif name == "Y":
            for t in targets:
                qc.y(t.value)
        elif name == "Z":
            for t in targets:
                qc.z(t.value)
        elif name == "S":
            for t in targets:
                qc.s(t.value)
        elif name == "S_DAG":
            for t in targets:
                qc.sdg(t.value)
        elif name == "H_XY":
            for t in targets:
                qc.ry(np.pi/4, t.value)
        elif name == "H_YZ":
            for t in targets:
                qc.ry(-np.pi/4, t.value)
        elif name == "T":
            for t in targets:
                qc.t(t.value)
        elif name == "T_DAG":
            for t in targets:
                qc.tdg(t.value)
        elif name == "SQRT_X":
            for t in targets:
                qc.sx(t.value)
        elif name == "SQRT_X_DAG":
            for t in targets:
                qc.sxdg(t.value)
        elif name == "CX" or name == "CNOT" or name == "ZCX":
            for i in range(0, len(targets), 2):
                qc.cx(targets[i].value, targets[i+1].value)
        elif name == "CZ":
            for i in range(0, len(targets), 2):
                qc.cz(targets[i].value, targets[i+1].value)
        elif name == "CY":
            for i in range(0, len(targets), 2):
                qc.cy(targets[i].value, targets[i+1].value)
        elif name == "SWAP":
            for i in range(0, len(targets), 2):
                qc.swap(targets[i].value, targets[i+1].value)
        elif name == "R":
            # for t in targets:
            #     qc.reset(t.value)
            moves += 1
        elif name == "M":
            # for t in targets:
            #     qc.measure(t.value, meas_idx)
            #     meas_idx += 1
            moves += 1
        elif name == "MX":
            # for t in targets:
            #     qc.h(t.value)
            #     qc.measure(t.value, meas_idx)
            #     qc.h(t.value)
            #     meas_idx += 1
            moves += 1
        elif name == "MY":
            # for t in targets:
            #     qc.sdg(t.value)
            #     qc.h(t.value)
            #     qc.measure(t.value, meas_idx)
            #     qc.h(t.value)
            #     qc.s(t.value)
            #     meas_idx += 1
            moves += 1
        elif name == "MR":
            moves += 1
        elif name == "TICK":
            # qc.barrier()
            pass
        elif name in ("DEPOLARIZE1", "DEPOLARIZE2", "X_ERROR", "Z_ERROR",
                       "DETECTOR", "OBSERVABLE_INCLUDE", "QUBIT_COORDS"):
            pass
        else:
            print(f"Warning: unsupported gate '{name}', skipping")
    
    if return_meas_reset_moves:
        return qc, moves
    
    return qc