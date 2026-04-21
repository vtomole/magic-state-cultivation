import stim
from typing import List
from .ghz_fxns import GHZstate
from .full_circuit_fxns import FullCircuit
from .noise_model import insert_circuit_errs

def full_circuit_f5(nm: float,
                dfinal: int,
                 latter_rounds: int = 3,
                 ghz_size: int = 3,
                 component_array: List = [1,1,1,1,1,1],
                 cultiv_only: bool = False,
                 prep: str = "unitstab",
                 neutralatom: bool = False,
                 r1: int = 3,
                 r2: int = 0,
                 ) -> stim.Circuit:

    rsc = FullCircuit(dx=dfinal,
                    dy=dfinal ,
                    glen=ghz_size,
                    smallsc=False,
                    basis="Y")
    
    # --- Injection: hook inject into Rot(3), then grow to Reg(3) ---
    if prep == "hookinj":
        prep_circuit = rsc.cstage_circ.d3rot_hookinj()
        prep_circuit += rsc.cstage_circ.d3_rot_to_reg()
    else:
        raise NotImplementedError

    rsc.qcircuit += insert_circuit_errs(prep_circuit, nm,
        valid = component_array[0], convert_nac=neutralatom)

    # --- Cultivation: two HXY checks at Reg(3) using GHZ(3) ---
    ghzcirc_3 = GHZstate(dfinal + 1, dfinal + 1, 3)
    rsc.qcircuit += ghzcirc_3.layout_ghz_state()
    rsc.update_ghz_circ(ghzcirc_3, 3)
    for _ in range(2):
        ghz_prep = rsc.ghzcirc.prepare_ghz_state()
        rsc.qcircuit += insert_circuit_errs(ghz_prep, nm,
            valid = component_array[1], convert_nac=neutralatom)
        noisy_check = rsc.cbasis_check()
        rsc.qcircuit += insert_circuit_errs(noisy_check, nm+(nm*neutralatom),
            valid = component_array[2], convert_nac=neutralatom)
        ghz_meas = rsc.ghzcirc.measure_ghz_state()
        rsc.qcircuit += insert_circuit_errs(ghz_meas, nm,
            valid = component_array[1], convert_nac=neutralatom)

    # --- Unitary grow Reg(3) -> Rot(5) ---
    grow_d3d5 = rsc.cstage_circ.grow_3u5r()
    rsc.qcircuit += insert_circuit_errs(grow_d3d5, nm,
        valid = component_array[3], convert_nac=neutralatom)

    # --- r1 stabilizer measurement rounds at Rot(5) ---
    for i in range(r1):
        rot5_stab = rsc.sc_stab_round(d_rest=5)
        rot5_stab += rsc.sc_detectors(curr_only=(i == 0),
                                        d_rest=5,
                                        ps_round=True)
        rsc.qcircuit += insert_circuit_errs(rot5_stab, nm,
            valid = component_array[4], convert_nac=neutralatom)
        
    # --- grow Rot(5) to Reg    (5) ---
    rsc.cstage_circ.d5_rot_to_reg()
    
    # --- Two HXY checks at Rot(5) using GHZ(5) ---
    ghzcirc_5 = GHZstate(dfinal + 2, dfinal + 2, 5)
    rsc.qcircuit += ghzcirc_5.layout_ghz_state()
    rsc.update_ghz_circ(ghzcirc_5, 5)
    for _ in range(2):
        ghz_prep = rsc.ghzcirc.prepare_ghz_state()
        rsc.qcircuit += insert_circuit_errs(ghz_prep, nm,
            valid = component_array[1], convert_nac=neutralatom)
        noisy_check = rsc.cbasis_check()
        rsc.qcircuit += insert_circuit_errs(noisy_check, nm+(nm*neutralatom),
            valid = component_array[2], convert_nac=neutralatom)
        ghz_meas = rsc.ghzcirc.measure_ghz_state()
        rsc.qcircuit += insert_circuit_errs(ghz_meas, nm,
            valid = component_array[1], convert_nac=neutralatom)
    
    # TODO: Maybe add reg(3) stabalizer measurments for r2 meas (But superior performance for (2, 0) Ref Fig A.12)
    if cultiv_only:
        return rsc.qcircuit
    
    # --- Escape: reset qubits and grow to Rot(dfin) ---
    d_rest = 5
    grow_final = rsc.larger_code_reset(d_rest=d_rest)

    # In total 5 rounds of stab measurement rounds Ref Page 22 Section E3
    # do first stabilizer msmt round on larger code
    grow_final += rsc.sc_stab_round()
    grow_final += rsc.sc_detectors(curr_only=True,
                                    first_round=True,
                                    d_rest=d_rest)

    # do latter stab msmts on larger code
    for _ in range(latter_rounds - 1):
        grow_final += rsc.sc_stab_round()
        grow_final += rsc.sc_detectors()

    rsc.qcircuit += insert_circuit_errs(grow_final, nm,
        valid = component_array[5], convert_nac=neutralatom)

    # last perfect stab msmt round for decoding
    rsc.qcircuit += rsc.sc_stab_round()
    rsc.qcircuit += rsc.sc_detectors()

    #5. add logical msmt
    rsc.qcircuit += rsc.logYMeas()
    print(len(rsc.qcircuit.shortest_graphlike_error()))

    return rsc.qcircuit

# circuit = full_circuit(nm = 0.001, ghz_size=5, dfinal=7, prep="hookinj")
# print(circuit)
