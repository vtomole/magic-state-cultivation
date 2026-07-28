
import cultiv
import stimcirq
import cirq
import numpy as np
from collections import Counter

from cultiv._construction.full_clifford_sim.main_compiled_fxns_fault5 import full_circuit_f5
from cultiv._construction.full_clifford_sim.main_complied_fxns import full_circuit

import sys
def measurement_moment_count(circuit: cirq.Circuit):
	"""
	Returns a tuple (num_cirq_moments, num_stimcirq_moments) where:
	  - num_cirq_moments: number of moments containing cirq.MeasurementGate
	  - num_stimcirq_moments: number of moments containing stimcirq.MeasureAndOrResetGate
	"""
	num_cirq_moments = 0
	num_stimcirq_moments = 0
	for moment in circuit:
		has_cirq_meas = any(
			hasattr(op, 'gate') and isinstance(op.gate, cirq.MeasurementGate)
			for op in moment.operations
		)
		has_stim_meas = any(
			hasattr(op, 'gate') and isinstance(op.gate, stimcirq.MeasureAndOrResetGate)
			for op in moment.operations
		)
		if has_cirq_meas:
			num_cirq_moments += 1
		if has_stim_meas:
			num_stimcirq_moments += 1
	return num_cirq_moments, num_stimcirq_moments

def make_stim_circuit(code_distance: int, fault_distance: int):
	if fault_distance == 5:
		stim_circuit = full_circuit_f5(nm=.001, prep='hookinj', dfinal=code_distance, cultiv_only=False).without_noise()
	elif fault_distance == 3:
		stim_circuit = full_circuit(nm=.001, prep='hookinj', dfinal=code_distance, cultiv_only=False).without_noise()
	else:
		raise ValueError
	stim_rep = stim_circuit.with_inlined_feedback()
	return stim_rep

def make_cirq_circuit(code_distance: int, fault_distance: int):
	stim_rep = make_stim_circuit(code_distance=code_distance, fault_distance=fault_distance)
	cirq_rep = stimcirq.stim_circuit_to_cirq_circuit(stim_rep)
	d2 = code_distance
	#### process Clifford circuit into nonClifford Circ
	edited_cirq_rep = cirq.Circuit()
	sevendiaggate = cirq.MatrixGate(np.array([[0,1],[1j,0]]))
	sodddiaggate = cirq.MatrixGate(np.array([[0,1],[-1j,0]]))

	cy_moments = 0
	ghzh_counter = 0
	for midx, moment in enumerate(cirq_rep):
		new_moment = []
		this_moment_is_CY = False

		for op in moment:

			if op.gate == cirq.ops.SingleQubitCliffordGate(_clifford_tableau=cirq.CliffordTableau( #ONLY FOR UNITARY PREP
				1, rs=np.array([True, False]), xs=np.array([[True], [True]]), zs=np.array([[False], [True]]))):
				new_moment.append(cirq.PhasedXPowGate(phase_exponent=3/4, exponent=1/2).on(op.qubits[0]))

			elif op.gate == cirq.ControlledGate(cirq.Y) or op.gate == cirq.CY: #add a CH gate layer
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
	return edited_cirq_rep

def dirty_count(circuit: cirq.Circuit):
    serial = Counter()
    parallel = Counter()
    for moment in circuit:
        gates = [op.gate for op in moment.operations if op not in cirq.GateFamily(cirq.I)]
        if gates and gates[0] is not None:
            if gates[0] in cirq.GateFamily(cirq.MeasurementGate):
                gates = [cirq.MeasurementGate] * len(gates)
            if gates[0] in cirq.GateFamily(stimcirq.MeasureAndOrResetGate):
                gates = [cirq.MeasurementGate] * len(gates)
            if isinstance(gates[0], cirq.ControlledGate) and gates[0] not in cirq.GateFamily(cirq.CCZ):
                gates = [cirq.CNOT] * len(gates)
            if isinstance(gates[0], cirq.ResetChannel):
                gates = [cirq.ResetChannel] * len(gates)
            serial += Counter(gates)
            parallel += Counter(gates[:1])
    op_map = {
            cirq.H: cirq.PhasedXZGate,
            cirq.S: cirq.PhasedXZGate,
            cirq.CNOT: cirq.CZ,
            (cirq.T**-1): cirq.PhasedXZGate,
            }
    for key, val in op_map.items():
        serial[val] += serial[key]
        del serial[key]
        parallel[val] += parallel[key]
        del parallel[key]
    return {'serial': serial, 'parallel': parallel}

if __name__ == "__main__":
    print(sys.argv)
    cd = int(sys.argv[1])
    fd = int(sys.argv[2])
    circuit = make_cirq_circuit(code_distance=cd, fault_distance=fd)
    result = dirty_count(circuit)
    print(f"Serial:\n\t{result['serial']}")
    print(sum(result['serial'].values()), len([op for op in circuit.all_operations()]))
    print(f"Parallel:\n\t{result['parallel']}")
    print(sum(result['parallel'].values()), len(circuit))
