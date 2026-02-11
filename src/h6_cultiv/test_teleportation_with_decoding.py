from Code614 import Code6
from concatenated_codes import Code36, Code216
import stim

def create_Code216_logical(start, anc, p, m, shift:int = 0):
    b = Code216(start, anc, p, m, shift)

    c = stim.Circuit()
    c += b.get_init_circ(with_aux_stab=True, measure_anc=False)
    # c += b.measure(with_stab=True)
    # c += b.add_stabalizers()

    return c, b

def create_Code36_magic_state(start, anc, p, m, target_block=None):
    b = Code36(start, anc, 15000, 15500, p, m)
    c = stim.Circuit()
    # Use resource qubit range starting at 11000/11500 to avoid collisions with
    # b1 (0-4071) and b2 (5000-11023) qubit ranges
    c += b.create_magic_state(with_aux_stab=True, add_observables=False,
                               target_block=target_block,
                               resource_dstart=13000, resource_astart=13500)

    return c, b

def get_b1_b2_targets(b1, b2):
    total_meas = b1.num_data_blocks * b1.size_of_unit_data_block

    # Offsets for decode measurements of b1.data_blocks[1:].
    # Measurement order (from end): b2 (216), b1.block0 (36), then b1.blocks 5..1 (each 36).
    def _decode_offset(block_idx: int) -> int:
        return total_meas + 36 + (6 - block_idx) * b1.size_of_unit_data_block

    targets = []
    offset = _decode_offset(4)
    targets += b1.data_blocks[4].z(offset, '00')
    for i in [0, 2, 4]:
        offset = total_meas - i*b1.size_of_unit_data_block
        targets += b2.data_blocks[i].z(offset, '00')
    c.append("OBSERVABLE_INCLUDE", targets, 0)

    targets = []
    offset = _decode_offset(4)
    targets += b1.data_blocks[4].z(offset, '01') 
    for i in [0, 2, 4]:
        offset = total_meas - i*b1.size_of_unit_data_block
        targets += b2.data_blocks[i].z(offset, '01')
    c.append("OBSERVABLE_INCLUDE", targets, 1)

    targets = []
    offset = _decode_offset(4)
    targets += b1.data_blocks[4].z(offset, '10')  
    for i in [0, 2, 4]:
        offset = total_meas - i*b1.size_of_unit_data_block
        targets += b2.data_blocks[i].z(offset, '10')
    c.append("OBSERVABLE_INCLUDE", targets, 2)

    targets = []
    offset = _decode_offset(4)
    targets += b1.data_blocks[4].z(offset, '11') 
    for i in [0, 2, 4]:
        offset = total_meas - i*b1.size_of_unit_data_block
        targets += b2.data_blocks[i].z(offset, '11')
    c.append("OBSERVABLE_INCLUDE", targets, 3)

    targets = []
    for i in [1, 3, 5]:
        offset = _decode_offset(i)
        targets += b1.data_blocks[i].z(offset, '00')
    for i in [1, 3, 5]:
        offset = total_meas - i*b1.size_of_unit_data_block
        targets += b2.data_blocks[i].z(offset, '00')
    c.append("OBSERVABLE_INCLUDE", targets, 4)

    targets = []
    for i in [1, 3, 5]:
        offset = _decode_offset(i)
        targets += b1.data_blocks[i].z(offset, '01')
    for i in [1, 3, 5]:
        offset = total_meas - i*b1.size_of_unit_data_block
        targets += b2.data_blocks[i].z(offset, '01')
    c.append("OBSERVABLE_INCLUDE", targets, 5)

    targets = []
    for i in [1, 3, 5]:
        offset = _decode_offset(i)
        targets += b1.data_blocks[i].z(offset, '10')
    for i in [1, 3, 5]:
        offset = total_meas - i*b1.size_of_unit_data_block
        targets += b2.data_blocks[i].z(offset, '10')
    c.append("OBSERVABLE_INCLUDE", targets, 6)

    targets = []
    for i in [1, 3, 5]:
        offset = _decode_offset(i)
        targets += b1.data_blocks[i].z(offset, '11')
    for i in [1, 3, 5]:
        offset = total_meas - i*b1.size_of_unit_data_block
        targets += b2.data_blocks[i].z(offset, '11')
    c.append("OBSERVABLE_INCLUDE", targets, 7)

    return c

# create 2 logical 216 states
start = 0
anc = 400
p = 0.001
m = 0.000
c = stim.Circuit()
# create b1 NOTE: need to add meas and stab
c1, b1 = create_Code216_logical(start, anc, p, m)
c += c1
# create b2 NOTE: need to add meas and stab
c2, b2 = create_Code216_logical(start+5000, anc+5000, p, m, shift=6000)
c += c2

# entangle the two pairs
c += b1.hadamard()
c += b1.cnot(b2)

c += b1.decode_to_lower_level()

c += b1.data_blocks[0].measure(with_stab=True)
c += b1.data_blocks[0].add_stabalizers()

c += b2.measure(with_stab=True)
c += b2.add_stabalizers()

c += get_b1_b2_targets(b1, b2)

errors = c.search_for_undetectable_logical_errors(
        dont_explore_detection_event_sets_with_size_above=16,
        dont_explore_edges_with_degree_above=16,
        dont_explore_edges_increasing_symptom_degree=True,
        canonicalize_circuit_errors=True,
    )

distance = len(errors)
print("Distance 216 code:", distance)
assert distance == 8
