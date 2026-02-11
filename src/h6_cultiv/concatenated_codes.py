from Code614 import Code6
from Code614 import get_init_circ as Code6_init
import stim 

class Code36:
    def __init__(self, dstart, dstart_aux, astart, astart_aux, p, m):
        self.dstart = dstart
        self.astart = astart
        self.dstart_aux = dstart_aux
        self.astart_aux = astart_aux
        self.num_data_blocks = 6
        self.num_aux_blocks = 2
        self.size_of_unit_data_block = 6
        self.size_of_unit_anc_block = 2
        self.p = p
        self.m = m
        self.data_blocks = [ 
            Code6(
                dstart + i* self.size_of_unit_data_block, 
                self.dstart_aux + i*self.size_of_unit_anc_block, self.p, self.m) 
            for i in range(self.num_data_blocks)]
        
        self.anc_blocks = [
            Code6(
                astart + i*self.size_of_unit_data_block,
                self.astart_aux + i*self.size_of_unit_anc_block, self.p, self.m)
            for i in range(self.num_aux_blocks)
        ]
        self.level = 2
        self.num_stabalizers = 2 **(self.level)
        self.num_observables = 2 **(self.level)
    
    def hadamard(self):
        c = stim.Circuit()
        # Hadamard applied only on data blocks in Code6
        for b in self.data_blocks:
            c += b.hadamard()
        return c
    
    def cnot(self, b2):
        assert isinstance(b2, Code36)
        c = stim.Circuit()
        for i in range(self.num_data_blocks):
            c += self.data_blocks[i].cnot(b2.data_blocks[i])
        return c
    
    def z(self, offset, flag):
        targets = []
        c = stim.Circuit()
        if flag == '00':
            for i in [0, 2, 4]:
                offset_curr = offset - i*self.size_of_unit_data_block
                targets += self.data_blocks[i].z(offset_curr, 0)
        elif flag == '01':
            for i in [0, 2, 4]:
                offset_curr = offset - i*self.size_of_unit_data_block
                targets += self.data_blocks[i].z(offset_curr, 1)
        elif flag == '10':
            for i in [1, 3, 5]:
                offset_curr = offset - i*self.size_of_unit_data_block
                targets += self.data_blocks[i].z(offset_curr, 0)
        elif flag == '11':
            for i in [1, 3, 5]:
                offset_curr = offset - i*self.size_of_unit_data_block
                targets += self.data_blocks[i].z(offset_curr, 1)
        return targets
    
    def measure(self, with_stab:bool=False):
        c = stim.Circuit()
        for b in self.data_blocks:
            c += b.measure()
            if with_stab:
                c += b.add_stabalizers()
        
        return c

    def add_stabalizers(self):
        c = stim.Circuit()
        targets = []
        total_meas = self.num_data_blocks * self.size_of_unit_data_block
        for i in [0, 1, 2, 3]:
                offset = total_meas - i*self.size_of_unit_data_block
                targets += self.data_blocks[i].z(offset, 0)
        c.append("DETECTOR", targets)

        targets = []
        for i in [0, 1, 2, 3]:
                offset = total_meas - i*self.size_of_unit_data_block
                targets += self.data_blocks[i].z(offset, 1)
        c.append("DETECTOR", targets)

        targets = []
        for i in [2, 3, 4, 5]:
                offset = total_meas - i*self.size_of_unit_data_block
                # based on z0z2z4
                targets += self.data_blocks[i].z(offset, 0)
        c.append("DETECTOR", targets)

        targets = []
        for i in [2, 3, 4, 5]:
                offset = total_meas - i*self.size_of_unit_data_block
                # based on z0z2z4
                targets += self.data_blocks[i].z( offset, 1)
        c.append("DETECTOR", targets)

        return c

    def add_logicals(self, log_offset:int= 0):
        c = stim.Circuit()
        total_meas = self.num_data_blocks * self.size_of_unit_data_block
        targets = []
        # z0z2z4
        for i in range(self.num_data_blocks):
            if i in [0, 2 , 4]:
                offset = total_meas - i*self.size_of_unit_data_block
                targets += self.data_blocks[i].z(offset, 0)
        c.append("OBSERVABLE_INCLUDE", targets, 0 + log_offset) 

        targets = []
        for i in range(self.num_data_blocks):
            if i in [0, 2 , 4]:
                offset = total_meas - i*self.size_of_unit_data_block
                targets += self.data_blocks[i].z(offset, 1)
        c.append("OBSERVABLE_INCLUDE", targets, 1+log_offset) 

        targets = []
        for i in range(self.num_data_blocks):
            if i in [1, 3, 5]:
                offset = total_meas - i*self.size_of_unit_data_block
                targets += self.data_blocks[i].z(offset, 0)
        c.append("OBSERVABLE_INCLUDE", targets, 2+log_offset)         

        targets = []
        for i in range(self.num_data_blocks):
            if i in [1, 3, 5]:
                offset = total_meas - i*self.size_of_unit_data_block
                targets += self.data_blocks[i].z(offset, 1)
        c.append("OBSERVABLE_INCLUDE", targets, 3+log_offset)

        return c
    
    def get_init_circ(self, with_aux_stab: bool = False):
        c = stim.Circuit()

        for i in range(self.num_data_blocks):
            c += self.data_blocks[i].get_init_circ(with_aux_stab=True)
        
        for i in range(self.num_aux_blocks):
            c += self.anc_blocks[i].get_init_circ(with_aux_stab=True)

        c += self.data_blocks[0].hadamard()
        c += self.data_blocks[2].hadamard()

        c += self.data_blocks[0].cnot(self.anc_blocks[0])
        c += self.data_blocks[2].cnot(self.anc_blocks[1])
        c += self.data_blocks[0].cnot(self.data_blocks[1])
        c += self.data_blocks[2].cnot(self.data_blocks[3])
        c += self.data_blocks[0].cnot(self.data_blocks[4])
        c += self.data_blocks[2].cnot(self.data_blocks[5])
        c += self.data_blocks[0].cnot(self.data_blocks[5])
        c += self.data_blocks[2].cnot(self.data_blocks[4])
        c += self.data_blocks[0].cnot(self.anc_blocks[0])
        c += self.data_blocks[2].cnot(self.anc_blocks[1])

        for anc in self.anc_blocks:
            c += anc.measure()
            # add detectors for ancilla
            if with_aux_stab:
                c += anc.add_stabalizers()
        
        return c
    
    def create_magic_state(
        self,
        with_aux_stab: bool = False,
        target_block=None,
        resource_dstart: int = 2000,
        resource_astart: int = 2500,
        bell_detectors: bool = False,
    ):
        c = stim.Circuit()
        c += self.data_blocks[0].get_dist_circ(with_aux_stab=True)
        c += self.data_blocks[1].get_dist_circ(with_aux_stab=True)

        for i in range(2, 6):
            c += self.data_blocks[i].get_init_circ(with_aux_stab=True)
        
        for i in range(0, 2):
            c+= self.anc_blocks[i].get_init_circ(with_aux_stab=True)
        
        resource_states = []
        resource_state_dstart = resource_dstart
        resource_state_astart = resource_astart
        for i in range(12):
            block = Code6(
                resource_state_dstart + self.size_of_unit_data_block*i,
                resource_state_astart + self.size_of_unit_anc_block*i, 
                self.p, self.m)
            resource_states.append(block)
            c += block.get_dist_circ(with_aux_stab=(target_block is None))
        
        c += self.data_blocks[2].hadamard()
        c += self.data_blocks[2].cnot(self.data_blocks[3])
        c += self.data_blocks[4].hadamard()
        c += self.data_blocks[4].cnot(self.data_blocks[5])
        
        # TODO: Check if this is needed in init circ
        # #Adding memory error to the two codeblocks not acted on by transversal gates in this depth 1
        # qss = [self.data_blocks[0].qubits, self.data_blocks[1].qubits]
        # c.append("Depolarize1", [q for qs in qss for q in qs], self.m)
        c += self.data_blocks[2].cnot(self.data_blocks[0])
        c += self.data_blocks[3].cnot(self.data_blocks[1])
        
        # qss = [self.data_blocks[4].qubits, self.data_blocks[5].qubits]
        # c.append("Depolarize1", [q for qs in qss for q in qs], self.m)
        c += self.data_blocks[0].cnot(self.data_blocks[4])
        c += self.data_blocks[1].cnot(self.data_blocks[5])
        
        # qss = [self.data_blocks[2].qubits, self.data_blocks[3].qubits]
        # c.append("Depolarize1", [q for qs in qss for q in qs], self.m)
        c += self.data_blocks[4].cnot(self.data_blocks[2])
        c += self.data_blocks[5].cnot(self.data_blocks[3])
        
        qss = [self.data_blocks[0].qubits, self.data_blocks[1].qubits]
        c.append("Depolarize1", [q for qs in qss for q in qs], self.m)

        #We rotate all logical qubits to |0> by performing Ry(-pi/2) rotations
        for i in range(0,6):
            c += self.data_blocks[i].rypi2(resource_states[i])
            # Skip resource stabilizers when used for Bell measurement (target_block),
            # since these detectors can be non-deterministic in that mode.
            if target_block is None:
                c += resource_states[i].add_stabalizers()

        #Now we can measure the Z value into the aux blocks
        c += self.anc_blocks[0].hadamard()
        c += self.anc_blocks[0].cnot(self.anc_blocks[1])

        c += self.anc_blocks[0].cz(self.data_blocks[0])
        c += self.anc_blocks[1].cz(self.data_blocks[1])
        # qss = [self.data_blocks[2].qubits, self.data_blocks[3].qubits, self.data_blocks[4].qubits, self.data_blocks[5].qubits]
        # c.append("Depolarize1", [q for qs in qss for q in qs], self.m)
    
        c += self.anc_blocks[0].cz(self.data_blocks[2])
        c += self.anc_blocks[1].cz(self.data_blocks[3])
        # qss = [self.data_blocks[0].qubits, self.data_blocks[1].qubits, self.data_blocks[4].qubits, self.data_blocks[5].qubits]
        # c.append("Depolarize1", [q for qs in qss for q in qs], self.m)
        
        c += self.anc_blocks[0].cz(self.data_blocks[4])
        c += self.anc_blocks[1].cz(self.data_blocks[5])
        # qss = [self.data_blocks[0].qubits, self.data_blocks[1].qubits, self.data_blocks[2].qubits, self.data_blocks[3].qubits]
        # c.append("Depolarize1", [q for qs in qss for q in qs], self.m)

        c += self.anc_blocks[0].cnot(self.anc_blocks[1])
        c += self.anc_blocks[0].hadamard()

        #We rotate all logical qubits back to |+> by performing Ry(+pi/2) rotations
        for i in range(0,6):
            c += self.data_blocks[i].rypi2(resource_states[i+6])
            # Skip resource stabilizers when used for Bell measurement (target_block)
            if target_block is None:
                c += resource_states[i+6].add_stabalizers()

        for anc in self.anc_blocks:
            c += anc.measure()
            if with_aux_stab:
                c += anc.add_stabalizers()

        # If target_block provided, CNOT before H+MR (makes the data measurement a Bell measurement)
        if target_block is not None:
            c += self.cnot(target_block)

        for i in range(0, 6):
            c += self.data_blocks[i].hadamard()
            c += self.data_blocks[i].measure()
            if target_block is None or bell_detectors:
                c += self.data_blocks[i].add_stabalizers()


        # Skip Code36-level detectors/observables in Bell mode unless explicitly requested.
        if target_block is not None and not bell_detectors:
            return c

        # --- Code36-level stabilizers & observables with resource corrections ---
        # Total measurements: 160
        #   Init: 40 | Ry(-) res: 36 | Ry(+) res: 36 | Anc: 12 | Data: 36
        # Ry(-) block j, qubit k: rec(-(120 - j*6 - k))
        # Ry(+) block j, qubit k: rec(-(84 - j*6 - k))
        total_meas = self.num_data_blocks * self.size_of_unit_data_block

        def _res_rec(block_j, qubit_indices):
            """Resource correction rec targets for both Ry(-) and Ry(+) rounds."""
            targets = []
            for k in qubit_indices:
                targets.append(stim.target_rec(-(120 - block_j*6 - k)))
                targets.append(stim.target_rec(-(84 - block_j*6 - k)))
            return targets

        # Code36 stabilizers (4 detectors)
        stab_config = [
            ([0, 1, 2, 3], [0, 2, 4]),
            ([0, 1, 2, 3], [1, 3, 5]),
            ([2, 3, 4, 5], [0, 2, 4]),
            ([2, 3, 4, 5], [1, 3, 5]),
        ]
        for block_indices, qubit_indices in stab_config:
            flag = 0 if qubit_indices == [0, 2, 4] else 1
            targets = []
            for i in block_indices:
                offset = total_meas - i * self.size_of_unit_data_block
                targets += self.data_blocks[i].z(offset, flag)
            for j in block_indices:
                targets += _res_rec(j, qubit_indices)
            c.append("DETECTOR", targets)

        # Ancilla logical detectors (CZ parity check)
        # aux0 measures Z parity of ALL data blocks (Bell pair funnels info to aux0)
        # Y correction from Ry(-) propagates through CZ → need Ry(-) resource corrections
        # Ry(+) happens after CZ, doesn't affect ancilla
        # aux0: rec(-48) to rec(-43), aux1: rec(-42) to rec(-37)
        def _ry_minus_all_blocks_rec(qubit_indices):
            targets = []
            for j in range(6):
                for k in qubit_indices:
                    targets.append(stim.target_rec(-(120 - j*6 - k)))
            return targets

        # aux0 Z_L0 + Ry(-) corrections for all 6 blocks
        targets = [stim.target_rec(-48), stim.target_rec(-46), stim.target_rec(-44)]
        targets += _ry_minus_all_blocks_rec([0, 2, 4])
        c.append("DETECTOR", targets)

        # aux0 Z_L1 + Ry(-) corrections for all 6 blocks
        targets = [stim.target_rec(-47), stim.target_rec(-45), stim.target_rec(-43)]
        targets += _ry_minus_all_blocks_rec([1, 3, 5])
        c.append("DETECTOR", targets)

        # aux1 Z_L0 (Bell pair reference arm, always |0>)
        c.append("DETECTOR", [stim.target_rec(-42), stim.target_rec(-40), stim.target_rec(-38)])

        # aux1 Z_L1 (Bell pair reference arm, always |0>)
        c.append("DETECTOR", [stim.target_rec(-41), stim.target_rec(-39), stim.target_rec(-37)])

        # Skip stabilizer/observable checks when used for Bell measurement (target_block)
        if target_block is not None:
            return c
        
        # Code36 observables (L0-L3)
        obs_config = [
            (0, [0, 2, 4], [0, 2, 4]),
            (1, [0, 2, 4], [1, 3, 5]),
            (2, [1, 3, 5], [0, 2, 4]),
            (3, [1, 3, 5], [1, 3, 5]),
        ]
        for obs_idx, block_indices, qubit_indices in obs_config:
            flag = 0 if qubit_indices == [0, 2, 4] else 1
            targets = []
            for i in block_indices:
                offset = total_meas - i * self.size_of_unit_data_block
                targets += self.data_blocks[i].z(offset, flag)
            for j in block_indices:
                targets += _res_rec(j, qubit_indices)
            c.append("OBSERVABLE_INCLUDE", targets, obs_idx)

        return c

    
class Code216():
    def __init__(self, dstart, astart, p2, m, shift:int = 0):
        self.dstart = dstart
        self.astart = astart
        self.dstart_dstart = dstart
        self.dstart_dstart_aux = 1500 + shift
        self.dstart_aux = 1000 + shift
        self.dstart_aux_aux = 2000 + shift
        self.astart_dstart = astart
        self.astart_dstart_aux = 3500 + shift
        self.astart_aux = 3000 + shift
        self.astart_aux_aux = 4000 + shift
        self.num_data_blocks = 6
        self.num_aux_blocks = 2
        self.size_of_unit_data_block = 36
        self.size_of_unit_anc_block = 12
        self.p = p2
        self.m = m
        self.data_blocks = [ 
            Code36(
                self.dstart_dstart + i* self.size_of_unit_data_block, 
                self.dstart_dstart_aux + i*self.size_of_unit_anc_block, 
                self.dstart_aux + i*self.size_of_unit_data_block,
                self.dstart_aux_aux + i*self.size_of_unit_anc_block,
                self.p, self.m) 
            for i in range(self.num_data_blocks)]
        
        self.anc_blocks = [
            Code36(
                self.astart_dstart + i* self.size_of_unit_data_block, 
                self.astart_dstart_aux + i*self.size_of_unit_anc_block, 
                self.astart_aux + i*self.size_of_unit_data_block,
                self.astart_aux_aux + i*self.size_of_unit_anc_block,
                self.p, self.m)
            for i in range(self.num_aux_blocks)
        ]
        self.level = 3
        self.num_stabalizers = 2 **(self.level)
        self.num_observables = 2 **(self.level)
    
    def hadamard(self):
        c = stim.Circuit()
        for b in self.data_blocks:
            c += b.hadamard()
        return c
    
    def cnot(self, b2):
        assert isinstance(b2, Code216)

        c = stim.Circuit()
        for i in range(self.num_data_blocks):
            c += self.data_blocks[i].cnot(b2.data_blocks[i])
        
        return c
    
    def measure(self, with_stab: bool = False):
        c = stim.Circuit()
        for b in self.data_blocks:
            c += b.measure(with_stab=with_stab)
            if with_stab:
                c += b.add_stabalizers()
        return c

    def add_stabalizers(self):
        c = stim.Circuit()
        targets = []
        total_meas = self.num_data_blocks * self.size_of_unit_data_block
        for i in [0, 1, 2, 3]:
            offset = total_meas - i*self.size_of_unit_data_block
            targets += self.data_blocks[i].z(offset, '00')
        c.append("DETECTOR", targets)

        targets = []
        for i in [0, 1, 2, 3]:
            offset = total_meas - i*self.size_of_unit_data_block
            targets += self.data_blocks[i].z(offset, '01')
        c.append("DETECTOR", targets)

        targets = []
        for i in [0, 1, 2, 3]:
            offset = total_meas - i*self.size_of_unit_data_block
            targets += self.data_blocks[i].z(offset, '10')
        c.append("DETECTOR", targets)

        targets = []
        for i in [0, 1, 2, 3]:
            offset = total_meas - i*self.size_of_unit_data_block
            targets += self.data_blocks[i].z(offset, '11')
        c.append("DETECTOR", targets)

        targets = []
        for i in [2, 3, 4, 5]:
            offset = total_meas - i*self.size_of_unit_data_block
            targets += self.data_blocks[i].z(offset, '00')
        c.append("DETECTOR", targets)

        targets = []
        for i in [2, 3, 4, 5]:
            offset = total_meas - i*self.size_of_unit_data_block
            targets += self.data_blocks[i].z(offset, '01')
        c.append("DETECTOR", targets)

        targets = []
        for i in [2, 3, 4, 5]:
            offset = total_meas - i*self.size_of_unit_data_block
            targets += self.data_blocks[i].z(offset, '10')
        c.append("DETECTOR", targets)

        targets = []
        for i in [2, 3, 4, 5]:
            offset = total_meas - i*self.size_of_unit_data_block
            targets += self.data_blocks[i].z(offset, '11')
        c.append("DETECTOR", targets)

        return c

    def add_logicals(self):
        c = stim.Circuit()
        total_meas = self.num_data_blocks * self.size_of_unit_data_block
        
        targets = []
        for i in [0, 2, 4]:
            offset = total_meas - i*self.size_of_unit_data_block
            targets += self.data_blocks[i].z(offset, '00') 
        c.append("OBSERVABLE_INCLUDE", targets, 0) 
        targets = []
        for i in [0, 2, 4]:
            offset = total_meas - i*self.size_of_unit_data_block
            targets += self.data_blocks[i].z(offset, '01') 
        c.append("OBSERVABLE_INCLUDE", targets, 1) 
        targets = []
        for i in [0, 2, 4]:
            offset = total_meas - i*self.size_of_unit_data_block
            targets += self.data_blocks[i].z(offset, '10') 
        c.append("OBSERVABLE_INCLUDE", targets, 2) 
        targets = []
        for i in [0, 2, 4]:
            offset = total_meas - i*self.size_of_unit_data_block
            targets += self.data_blocks[i].z(offset, '11') 
        c.append("OBSERVABLE_INCLUDE", targets, 3) 

        targets = []
        for i in [1, 3, 5]:
            offset = total_meas - i*self.size_of_unit_data_block
            targets += self.data_blocks[i].z(offset, '00') 
        c.append("OBSERVABLE_INCLUDE", targets, 4) 
        targets = []
        for i in [1, 3, 5]:
            offset = total_meas - i*self.size_of_unit_data_block
            targets += self.data_blocks[i].z(offset, '01') 
        c.append("OBSERVABLE_INCLUDE", targets, 5) 
        targets = []
        for i in [1, 3, 5]:
            offset = total_meas - i*self.size_of_unit_data_block
            targets += self.data_blocks[i].z(offset, '10') 
        c.append("OBSERVABLE_INCLUDE", targets, 6) 
        targets = []
        for i in [1, 3, 5]:
            offset = total_meas - i*self.size_of_unit_data_block
            targets += self.data_blocks[i].z(offset, '11') 
        c.append("OBSERVABLE_INCLUDE", targets, 7) 
        return c
        
    def get_init_circ(self, with_aux_stab: bool = False, measure_anc: bool = True):
        c = stim.Circuit()

        for i in range(self.num_data_blocks):
            c += self.data_blocks[i].get_init_circ(with_aux_stab=True)
        
        for i in range(self.num_aux_blocks):
            c += self.anc_blocks[i].get_init_circ(with_aux_stab=True)

        c += self.data_blocks[0].hadamard()
        c += self.data_blocks[2].hadamard()

        c += self.data_blocks[0].cnot(self.anc_blocks[0])
        c += self.data_blocks[2].cnot(self.anc_blocks[1])
        c += self.data_blocks[0].cnot(self.data_blocks[1])
        c += self.data_blocks[2].cnot(self.data_blocks[3])
        c += self.data_blocks[0].cnot(self.data_blocks[4])
        c += self.data_blocks[2].cnot(self.data_blocks[5])
        c += self.data_blocks[0].cnot(self.data_blocks[5])
        c += self.data_blocks[2].cnot(self.data_blocks[4])
        c += self.data_blocks[0].cnot(self.anc_blocks[0])
        c += self.data_blocks[2].cnot(self.anc_blocks[1])

        if measure_anc:
            for anc in self.anc_blocks:
                c += anc.measure(with_stab=True)
                # add detectors for ancilla
                if with_aux_stab:
                    c += anc.add_stabalizers()
        
        return c
    
    def decode_to_lower_level(self, with_stab: bool = True):
        # deocode block level: self.level-1 using decoding circuit
        c = stim.Circuit()

        c += self.data_blocks[2].cnot(self.anc_blocks[1])
        c += self.data_blocks[0].cnot(self.anc_blocks[0])
        c += self.data_blocks[2].cnot(self.data_blocks[4])
        c += self.data_blocks[0].cnot(self.data_blocks[5])
        c += self.data_blocks[2].cnot(self.data_blocks[5])
        c += self.data_blocks[0].cnot(self.data_blocks[4])
        c += self.data_blocks[2].cnot(self.data_blocks[3])
        c += self.data_blocks[0].cnot(self.data_blocks[1])
        c += self.data_blocks[2].cnot(self.anc_blocks[1])
        c += self.data_blocks[0].cnot(self.anc_blocks[0])

        c += self.data_blocks[0].hadamard()
        c += self.data_blocks[2].hadamard()

        for anc in self.anc_blocks:
            c += anc.measure(with_stab=True)
            # add detectors for ancilla
            # if with_aux_stab:
            if with_stab:
                c += anc.add_stabalizers()

        for b in self.data_blocks[1:]:
            c += b.measure(with_stab=with_stab)
            if with_stab:
                c += b.add_stabalizers()
        return c
        

        
