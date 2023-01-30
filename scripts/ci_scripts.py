
import os

scripts = ['bevel_gears+gear.py',
           'reduction_gear_2_stages.py',
           'dynamic_positions/crank_rod.py',
           'unidimensional/double_ball_bearings.py'
            ]

for script_name in scripts:
    print('\n## Executing script {}'.format(script_name))

    exec(open(script_name).read())
