from Code614 import Code6
from concatenated_codes import Code36, Code216
import stim 

def test_Code6_distance():
    start = 0
    anc = 200
    p = 0.001
    m = 0.001
    b = Code6(start, anc, p, m)

    c = stim.Circuit()
    c += b.get_init_circ(with_aux_stab=True)

    c += b.measure()
    c += b.add_stabalizers()
    c += b.add_logicals()

    distance = len(c.shortest_graphlike_error())

    assert distance == 2

def test_Code36_distance():
    start = 0
    anc = 400
    p = 0.001
    m = 0.001

    b = Code36(start, anc, 1000, 1500, p, m)

    c = stim.Circuit()
    c += b.get_init_circ(with_aux_stab=True)

    c += b.measure(with_stab=True)
    c += b.add_stabalizers()
    c += b.add_logicals()

    distance = len(c.shortest_graphlike_error())
    assert distance == 4

def test_Code36_magic_state():
    start = 0
    anc = 400
    p = 0.001
    m = 0.001
    b = Code36(start, anc, 1000, 1500, p, m)
    c = stim.Circuit()
    c += b.create_magic_state(with_aux_stab=True)
    
    distance = len(c.shortest_graphlike_error())
    assert distance == 4

def test_Code216_distance():
    start = 0
    anc = 400
    p = 0.001
    m = 0.001
    b = Code216(start, anc, p, m)

    c = stim.Circuit()
    c += b.get_init_circ(with_aux_stab=True)
    c += b.measure(with_stab=True)
    c += b.add_stabalizers()
    c += b.add_logicals()

    errors = c.search_for_undetectable_logical_errors(
        dont_explore_detection_event_sets_with_size_above=16,
        dont_explore_edges_with_degree_above=16,
        dont_explore_edges_increasing_symptom_degree=True,
        canonicalize_circuit_errors=True,
    )
    distance = len(errors)
    print("Distance 216 code:", distance)
    assert distance == 8

def test_Code216_distance_post_decoding():
    start = 0
    anc = 400
    p = 0.001
    m = 0.001
    b = Code216(start, anc, p, m)

    c = stim.Circuit()
    c += b.get_init_circ(with_aux_stab=True, measure_anc=False)
    c += b.decode_to_lower_level()
    # c += b.measure(with_stab=True)

    
    c += b.data_blocks[0].measure(with_stab=True)
    c += b.data_blocks[0].add_stabalizers()
    
    c += b.data_blocks[0].add_logicals()

    errors = c.search_for_undetectable_logical_errors(
        dont_explore_detection_event_sets_with_size_above=16,
        dont_explore_edges_with_degree_above=16,
        dont_explore_edges_increasing_symptom_degree=True,
        canonicalize_circuit_errors=True,
    )
    distance = len(errors)
    assert distance == 4

test_Code216_distance()
test_Code216_distance_post_decoding()
test_Code36_magic_state()