import numpy as np
import sys

reweight_card = """change mode NLO     
change helicity False 
change rwgt_dir rwgt
"""

c2v_values = np.linspace(-2, 2, 41)

for c2v in c2v_values:
    c2v = round(c2v, 2)
    reweight_card += f"""# C2V = {c2v}
launch --rewgt_name=scan_CV_1p0_C2V_{str(c2v).replace(".", "p").replace("-", "m")}_C3_1p0
set NEW 1 1.0
set NEW 2 {c2v}
set NEW 3 1.0
"""

with open(f"{sys.argv[1]}_4f_LO_reweight_card.dat", "w") as f:
    f.write(reweight_card)