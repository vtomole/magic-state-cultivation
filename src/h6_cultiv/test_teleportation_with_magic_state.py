from Code614 import Code6
from concatenated_codes import Code36, Code216
import stim

# Set to a set like {0} or {3} to include only specific observables.
# Use None to include all.
OBS_ONLY = {0}

def create_Code216_logical(start, anc, p, m, shift:int = 0):
    b = Code216(start, anc, p, m, shift)

    c = stim.Circuit()
    c += b.get_init_circ(with_aux_stab=True, measure_anc=False)
    # c += b.measure(with_stab=True)
    # c += b.add_stabalizers()

    return c, b

def create_Code36_magic_state(start, anc, p, m, target_block=None, with_aux_stab: bool = True, bell_detectors: bool = False):
    b = Code36(start, anc, 15000, 15500, p, m)
    c = stim.Circuit()
    # Use resource qubit range starting at 11000/11500 to avoid collisions with
    # b1 (0-4071) and b2 (5000-11023) qubit ranges
    c += b.create_magic_state(
        with_aux_stab=with_aux_stab,
        target_block=target_block,
        resource_dstart=13000,
        resource_astart=13500,
        bell_detectors=bell_detectors,
    )

    return c, b

def bell_x_correction_targets(b1, b_ms, flag):
    """Bell X-correction targets for Z-logical observables/detectors."""
    total_meas = b1.num_data_blocks * b1.size_of_unit_data_block
    block_indices = [0, 2, 4] if flag in ("00", "01") else [1, 3, 5]
    qubit_indices = [0, 2, 4] if flag in ("00", "10") else [1, 3, 5]

    targets = []

    # b_ms data measurements (X-basis via H before MR)
    b_ms_offset = total_meas + 2 * b1.size_of_unit_data_block  # 216 + 72 = 288
    targets += b_ms.z(b_ms_offset, flag)

    # b_ms aux measurements from distillation (one pair per block)
    # aux0 offset base = 412, aux1 base = 411
    for i in block_indices:
        targets.append(stim.target_rec(-(412 - 2 * i)))
        targets.append(stim.target_rec(-(411 - 2 * i)))

    # Resource data corrections (Ry(-) and Ry(+))
    for i in block_indices:
        for k in qubit_indices:
            targets.append(stim.target_rec(-(372 - i * 6 - k)))  # Ry(-)
            targets.append(stim.target_rec(-(336 - i * 6 - k)))  # Ry(+)

    # Resource aux measurements (one pair per resource block, Ry(-) and Ry(+))
    for i in block_indices:
        targets.append(stim.target_rec(-(396 - 2 * i)))  # Ry(-) aux0
        targets.append(stim.target_rec(-(395 - 2 * i)))  # Ry(-) aux1
        targets.append(stim.target_rec(-(384 - 2 * i)))  # Ry(+) aux0 (block i+6)
        targets.append(stim.target_rec(-(383 - 2 * i)))  # Ry(+) aux1 (block i+6)

    return targets


def add_b2_stabalizers_bell(b1, b2, b_ms):
    """Custom b2 stabilizers in Bell mode with Bell X corrections."""
    c = stim.Circuit()
    total_meas = b2.num_data_blocks * b2.size_of_unit_data_block

    def _emit(blocks, flag):
        targets = []
        for i in blocks:
            offset = total_meas - i * b2.size_of_unit_data_block
            targets += b2.data_blocks[i].z(offset, flag)
        c.append("DETECTOR", targets)

    for flag in ["00", "01", "10", "11"]:
        _emit([0, 1, 2, 3], flag)
    for flag in ["00", "01", "10", "11"]:
        _emit([2, 3, 4, 5], flag)

    return c


def get_b1_b2_targets(b1, b2, b_ms):
    total_meas = b1.num_data_blocks * b1.size_of_unit_data_block
    # Offsets for decode measurements of b1.data_blocks[1:].
    # Verified (from end of circuit) for this magic-state circuit:
    # b1[1]=592, b1[2]=556, b1[3]=520, b1[4]=484, b1[5]=448.
    decode_offsets = {1: 592, 2: 556, 3: 520, 4: 484, 5: 448}

    def _b1_decode_correction_targets(flag):
        targets = []
        for k in [1, 3, 5]:
            targets += b1.data_blocks[k].z(decode_offsets[k], flag)
        return targets

    def _bell_x_correction_targets(flag):
        return bell_x_correction_targets(b1, b_ms, flag)

    include_only = OBS_ONLY

    def _emit(obs_idx, targets):
        if include_only is None or obs_idx in include_only:
            c.append("OBSERVABLE_INCLUDE", targets, obs_idx)

    targets = []
    for i in [0, 2, 4]:
        offset = total_meas - i*b1.size_of_unit_data_block
        targets += b2.data_blocks[i].z(offset, '00')
    # Decode correction from b1 block4 (from pauli_prop L0 propagation)
    targets += b1.data_blocks[4].z(decode_offsets[4], '00')
    targets += _bell_x_correction_targets('00')
    _emit(0, targets)

    targets = []
    for i in [0, 2, 4]:
        offset = total_meas - i*b1.size_of_unit_data_block
        targets += b2.data_blocks[i].z(offset, '01')
    targets += b1.data_blocks[4].z(decode_offsets[4], '01')
    targets += _bell_x_correction_targets('01')
    _emit(1, targets)

    targets = []
    for i in [0, 2, 4]:
        offset = total_meas - i*b1.size_of_unit_data_block
        targets += b2.data_blocks[i].z(offset, '10')
    targets += b1.data_blocks[4].z(decode_offsets[4], '10')
    targets += _bell_x_correction_targets('10')
    _emit(2, targets)

    targets = []
    for i in [0, 2, 4]:
        offset = total_meas - i*b1.size_of_unit_data_block
        targets += b2.data_blocks[i].z(offset, '11')
    targets += b1.data_blocks[4].z(decode_offsets[4], '11')
    targets += _bell_x_correction_targets('11')
    _emit(3, targets)

    targets = []
    for i in [1, 3, 5]:
        offset = total_meas - i*b1.size_of_unit_data_block
        targets += b2.data_blocks[i].z(offset, '00')
    # targets += _b1_decode_correction_targets('00')
    _emit(4, targets)

    targets = []
    for i in [1, 3, 5]:
        offset = total_meas - i*b1.size_of_unit_data_block
        targets += b2.data_blocks[i].z(offset, '01')
    # targets += _b1_decode_correction_targets('01')
    _emit(5, targets)

    targets = []
    for i in [1, 3, 5]:
        offset = total_meas - i*b1.size_of_unit_data_block
        targets += b2.data_blocks[i].z(offset, '10')
    # targets += _b1_decode_correction_targets('10')
    _emit(6, targets)

    targets = []
    for i in [1, 3, 5]:
        offset = total_meas - i*b1.size_of_unit_data_block
        targets += b2.data_blocks[i].z(offset, '11')
    # targets += _b1_decode_correction_targets('11')
    _emit(7, targets)

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

c += b1.decode_to_lower_level(with_stab=False)

c_ms, b_ms = create_Code36_magic_state(
    20000,
    20400,
    p,
    m,
    target_block=b1.data_blocks[0],
    with_aux_stab=False,
    bell_detectors=False,
)
c += c_ms

c += b1.data_blocks[0].measure(with_stab=False)
c += b1.data_blocks[0].add_stabalizers()

c += b2.measure(with_stab=False)
c += add_b2_stabalizers_bell(b1, b2, b_ms)

c += get_b1_b2_targets(b1, b2, b_ms)


errors = c.search_for_undetectable_logical_errors(
        dont_explore_detection_event_sets_with_size_above=16,
        dont_explore_edges_with_degree_above=16,
        dont_explore_edges_increasing_symptom_degree=True,
        canonicalize_circuit_errors=True,
    )

distance = len(errors)
print("Distance 216 code:", distance)
assert distance == 8
