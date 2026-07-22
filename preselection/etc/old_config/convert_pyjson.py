import xsec_ref
import json

import sys
keyRun=sys.argv[1]
suffix={"run2":"_13TeV_4f_LO",
        "run3":"_13p6TeV_4f_LO"}

#print(xsec_ref.xsec_dict.keys())

runX=xsec_ref.xsec_dict[keyRun]

newdict={}
#json_string = json.dumps(data, indent=4)

for key in runX["bkg"]:
    newdict[key.split("_Tune")[0]]=runX["bkg"][key]

for key in runX["sig"]:
    newdict[key.split("_Tune")[0].replace("c3","C3").replace("c2v", "C2V_")+suffix[keyRun]]=runX["sig"][key]

json_string = json.dumps(newdict, indent=4)
print(json_string)

