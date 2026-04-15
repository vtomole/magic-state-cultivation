#NOTE: This requires integration with full_circuit code from MSC_foldedH
import numpy as np
import cirq
import stimcirq

from qiskit import QuantumCircuit

from full_clifford_sim.main_complied_fxns import full_circuit
from full_clifford_sim.cirq_utilities import *
## take DEM from ref_circuit, and decode shots from sample circutis


if __name__ == '__main__':

    num_shots_per_sverr = 20
    num_sverr_shots = 10

    p = 0.001
    errmodel = 'uniform'
    basis =  'HXY'
    if basis == 'HXY':
        true_bloch_rep = [1/np.sqrt(2),1/np.sqrt(2),0]
        magic = True
    ghz_size = 3
    errD = False
    cX =  True
    d2 = 7
    prep = 'hookinj'

    if prep == "unitstab":
        endofinj = 8
    elif prep == "hookinj":
        endofinj = 12
    elif prep == "optunit":
        endofinj = 2
    else:
        raise ValueError("injection stratey not well defined")

    #create a reference circuit used for matching
    ref_circuit = full_circuit(nm=p,
                         dfinal=d2,
                         ps_on_d3=1,
                         prep=prep,
                         neutralatom=False,
                         cultiv_only=True,
                         )

    criq_rep = stimcirq.stim_circuit_to_cirq_circuit(remove_QCs(ref_circuit, p, reverse_Y=False).without_noise())

    if magic: #break down circuit

        #### process Clifford circuit into nonClifford Circ
        edited_cirq_rep = cirq.Circuit()
        sevendiaggate = cirq.MatrixGate(np.array([[0,1],[1j,0]]))
        sodddiaggate = cirq.MatrixGate(np.array([[0,1],[-1j,0]]))

        cy_moments = 0
        ghzh_counter = 0

        for midx, moment in enumerate(criq_rep):
            new_moment = []
            this_moment_is_CY = False

            for op in moment:

                if op.gate == cirq.ops.SingleQubitCliffordGate(_clifford_tableau=cirq.CliffordTableau( #ONLY FOR UNITARY PREP
                    1, rs=np.array([True, False]), xs=np.array([[True], [True]]), zs=np.array([[False], [True]]))):
                    new_moment.append(cirq.PhasedXPowGate(phase_exponent=3/4, exponent=1/2).on(op.qubits[0]))

                elif op.gate == cirq.ControlledGate(cirq.Y): #add a CH gate layer
                    if this_moment_is_CY:
                        pass # print(f"{op} at moment {midx} is removed")
                    else:
                        ghzqub1 = op.qubits[0].x

                        if cy_moments % 3 == 0: #gate layer 1: CCZs
                            ccz1 = cirq.CCZ(cirq.LineQubit(ghzqub1) , cirq.LineQubit(4), cirq.LineQubit(4*d2))
                            ccz2 = cirq.CCZ(cirq.LineQubit(ghzqub1+1) ,cirq.LineQubit(2*d2+6) ,cirq.LineQubit(6*d2+2) )
                            ccz3 = cirq.CCZ(cirq.LineQubit(ghzqub1+2) , cirq.LineQubit(4*d2+8) , cirq.LineQubit(8*d2+4))
                            edited_cirq_rep.append(cirq.Moment([ccz1, ccz2, ccz3]))

                        elif cy_moments % 3 == 1: #gate layer 2 diagonal gates
                            cs1 =  sevendiaggate( cirq.LineQubit(0)).controlled_by(cirq.LineQubit(ghzqub1))
                            csdag2 = sodddiaggate( cirq.LineQubit(2*d2+2)).controlled_by(cirq.LineQubit(ghzqub1+1))
                            cs3 = sevendiaggate( cirq.LineQubit(4*d2+4)).controlled_by(cirq.LineQubit(ghzqub1+2))
                            edited_cirq_rep.append(cirq.Moment([cs1, csdag2, cs3]))

                        elif cy_moments % 3 == 2: #gate layer 3: 2 diags one off diag
                            csdag1 = sodddiaggate( cirq.LineQubit(6*d2+6)).controlled_by(cirq.LineQubit(ghzqub1))
                            cs2 = sevendiaggate( cirq.LineQubit(8*d2+8)).controlled_by(cirq.LineQubit(ghzqub1+1))
                            ccz3last = cirq.CCZ(cirq.LineQubit(ghzqub1+2) , cirq.LineQubit(8) , cirq.LineQubit(8*d2))
                            edited_cirq_rep.append(cirq.Moment([csdag1, cs2, ccz3last]))

                        this_moment_is_CY = True

                elif op.gate == cirq.H and op.qubits[0].x == 2*(d2+1)**2 +1:
                    if ghzh_counter % 2 == 1:
                        edited_cirq_rep.append(cirq.Moment(cirq.ZPowGate(exponent=-0.25)(op.qubits[0])))
                        new_moment.append(op) #("Adding T before meas")
                    else:
                        new_moment.append(op)
                    ghzh_counter+=1

                else:
                    new_moment.append(op)

            if this_moment_is_CY:
                cy_moments+=1
            edited_cirq_rep.append(cirq.Moment(new_moment))
    else: #no need to reprocess circuit
        edited_cirq_rep = criq_rep

    # ── Convert edited_cirq_rep → Qiskit QuantumCircuit ───────────────────────
    _sevendiag_mat = np.array([[0, 1], [1j,  0]])
    _sodddiag_mat  = np.array([[0, 1], [-1j, 0]])

    n_qubits = max(q.x for q in edited_cirq_rep.all_qubits()) + 1
    n_cbits  = sum(len(op.qubits) for moment in edited_cirq_rep
                   for op in moment if isinstance(op.gate, cirq.MeasurementGate))

    qiskit_circ = QuantumCircuit(n_qubits, n_cbits)
    creg_idx = 0

    for moment in edited_cirq_rep:
        for op in moment:
            qs = [q.x for q in op.qubits]

            if op.gate == cirq.CCZ:
                qiskit_circ.ccz(qs[0], qs[1], qs[2])

            elif op.gate == cirq.H:
                qiskit_circ.h(qs[0])

            elif op.gate == cirq.X:
                qiskit_circ.x(qs[0])
            elif op.gate == cirq.Y:
                qiskit_circ.y(qs[0])
            elif op.gate == cirq.Z:
                qiskit_circ.z(qs[0])

            elif op.gate == cirq.S:
                qiskit_circ.s(qs[0])
            elif op.gate == cirq.S**-1:
                qiskit_circ.sdg(qs[0])

            elif op.gate == cirq.T:
                qiskit_circ.t(qs[0])
            elif op.gate == cirq.T**-1 or (
                    isinstance(op.gate, cirq.ZPowGate) and np.isclose(op.gate.exponent, -0.25)):
                qiskit_circ.tdg(qs[0])

            elif op.gate == cirq.CNOT or op.gate == cirq.CX:
                qiskit_circ.cx(qs[0], qs[1])
            elif op.gate == cirq.CZ:
                qiskit_circ.cz(qs[0], qs[1])

            elif op.gate == cirq.I:
                qiskit_circ.id(qs[0])

            elif isinstance(op.gate, cirq.PhasedXPowGate):
                # PhasedXPowGate(3/4, 1/2) = Z^(3/4)·X^(1/2)·Z^(-3/4)
                #                          = P(3π/4) · SX · P(-3π/4)
                p_exp = op.gate.phase_exponent
                t_exp = op.gate.exponent
                qiskit_circ.p(-p_exp * np.pi, qs[0])
                qiskit_circ.sx(qs[0])   # X^(1/2) up to global phase
                qiskit_circ.p( p_exp * np.pi, qs[0])

            elif (isinstance(op.gate, cirq.ControlledGate)
                  and isinstance(op.gate.sub_gate, cirq.MatrixGate)):
                # ctrl = qs[0], tgt = qs[1]
                sub_mat = cirq.unitary(op.gate.sub_gate)
                if np.allclose(sub_mat, _sevendiag_mat):
                    # sevendiaggate = S·X  →  controlled: CX then CP(+π/2)
                    qiskit_circ.cx(qs[0], qs[1])
                    qiskit_circ.cp(np.pi / 2, qs[0], qs[1])
                elif np.allclose(sub_mat, _sodddiag_mat):
                    # sodddiaggate = S†·X  →  controlled: CX then CP(-π/2)
                    qiskit_circ.cx(qs[0], qs[1])
                    qiskit_circ.cp(-np.pi / 2, qs[0], qs[1])
                else:
                    raise ValueError(f"Unhandled controlled MatrixGate: {op}")

            elif isinstance(op.gate, cirq.ResetChannel):
                qiskit_circ.reset(qs[0])

            elif isinstance(op.gate, cirq.MeasurementGate):
                n = len(qs)
                qiskit_circ.measure(qs, list(range(creg_idx, creg_idx + n)))
                creg_idx += n

            else:
                raise ValueError(f"Unhandled gate: {op.gate!r} (type: {type(op.gate).__name__})")

        qiskit_circ.barrier()

    from qiskit import qasm2
    qasm_str = qasm2.dumps(qiskit_circ)
    with open("folded_cultivation_circ_qasm.txt", 'w') as f:
        f.write(qasm_str)
