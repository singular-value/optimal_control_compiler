"""
- An example to run jobs on midway in batch
"""

import subprocess
import numpy as np
import os
import re

PATH = 'qasm/sqrt/'

angs = []
for filename in os.listdir(os.getcwd()+'/'+PATH):
    ang = re.compile('(\S+).qasm').search(filename)
    if(ang):
        angs.append(ang.group(1))

for k in angs:
    for i in np.linspace(9, 27, 6):
        content = """#!/bin/bash

#SBATCH --job-name=time%s-%s.batch
#SBATCH --output=time%s%s.out
#SBATCH --error=time%s%s.err
#SBATCH --partition broadwl
#SBATCH --time=31:10:00
#SBATCH --mem=55000

module load python_ucs4/2.7.12

python gen_pulse.py --circ %s.qasm --time %s
    """ % ("{0:3f}".format(i),k, "{0:3f}".format(i),k,  "{0:3f}".format(i), k, PATH+k, "{0:3f}".format(i))
        with open("new_coupling_time" + "{0:3f}".format(i)+k+".batch", "w") as f:
            f.write(content)

for k in angs:
    for i in np.linspace(9, 27, 6):
        subprocess.call(["sbatch", "new_coupling_time" + "{0:3f}".format(i)+k+".batch"])
