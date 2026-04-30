import cirq
try:
    import stimcirq
except ImportError:
    stimcirq = None

def measurement_moment_count(circuit):
    """
    Returns a tuple (num_cirq_moments, num_stimcirq_moments) where:
      - num_cirq_moments: number of moments containing cirq.MeasurementGate
      - num_stimcirq_moments: number of moments containing stimcirq.MeasureAndOrResetGate
    """
    num_cirq_moments = 0
    num_stimcirq_moments = 0
    for moment in circuit:
        has_cirq_meas = any(isinstance(op.gate, cirq.MeasurementGate) for op in moment.operations if hasattr(op, 'gate'))
        has_stim_meas = False
        if stimcirq is not None:
            has_stim_meas = any(isinstance(op.gate, stimcirq.MeasureAndOrResetGate) for op in moment.operations if hasattr(op, 'gate'))
        if has_cirq_meas:
            num_cirq_moments += 1
        if has_stim_meas:
            num_stimcirq_moments += 1
    return num_cirq_moments, num_stimcirq_moments
