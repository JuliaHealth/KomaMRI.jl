import math
import numpy as np
import pypulseq as pp


# Set system limits
system = pp.Opts(
    max_grad=28,
    grad_unit='mT/m',
    max_slew=150,
    slew_unit='T/m/s',
    rf_ringdown_time=20e-6,
    rf_dead_time=100e-6,
    adc_dead_time=10e-6,
)

seq = pp.Sequence(system)  # Create a new sequence object
seq.add_block(pp.make_label(label='REV', type='SET', value=0),pp.make_label(label='ECO', type='SET', value=0))

seq.add_block(pp.make_label(label='LIN', type='INC', value=1),pp.make_label(label='ECO', type='SET', value=0))

seq.add_block(pp.make_label(label='LIN', type='INC', value=1),pp.make_label(label='ECO', type='SET', value=2))

seq.add_block(pp.make_label(label='LIN', type='INC', value=1),pp.make_label(label='ECO', type='SET', value=1))

seq.add_block(pp.make_label(label='LIN', type='INC', value=1),pp.make_label(label='ECO', type='SET', value=2))

seq.add_block(pp.make_label(label='LIN', type='SET', value=0),pp.make_label(label='ECO', type='SET', value=1))

seq.write("label_test.seq")