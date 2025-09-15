#!/usr/bin/env python3
import numpy as np
import onnxruntime as ort
import argparse
import sys
from pathlib import Path
import uproot

def load_model(model_path):
    session = ort.InferenceSession(model_path)
    return session

def load_one_event(input_file, event_idx=0):
    try:
        with uproot.open(input_file + ":Events") as tree:
            data = tree.arrays(library="ak", entry_start=event_idx, entry_stop=event_idx+1)
            
            ak4_pt = data["Jet_pt"][0].to_numpy()
            ak4_eta = data["Jet_eta"][0].to_numpy() 
            ak4_phi = data["Jet_phi"][0].to_numpy() 
            ak4_mass = data["Jet_mass"][0].to_numpy() 
            ak4_isTightBTag = data["Jet_isTightBTag"][0].to_numpy() 
            ak4_isMediumBTag = data["Jet_isMediumBTag"][0].to_numpy()
            ak4_isLooseBTag = data["Jet_isLooseBTag"][0].to_numpy()
            
            ak8_pt = data["FatJet_pt"][0].to_numpy() 
            ak8_eta = data["FatJet_eta"][0].to_numpy() 
            ak8_phi = data["FatJet_phi"][0].to_numpy() 
            ak8_mass = data["FatJet_mass"][0].to_numpy() 
            ak8_nConstituents = data["FatJet_nConstituents"][0].to_numpy() 
            ak8_XbbScore = data["FatJet_globalParT3_Xbb"][0].to_numpy() 
            ak8_XqqScore = data["FatJet_globalParT3_Xqq"][0].to_numpy() 
            ak8_XccScore = data["FatJet_globalParT3_Xcc"][0].to_numpy() 
            ak8_XcsScore = data["FatJet_globalParT3_Xcs"][0].to_numpy() 
            ak8_XqcdScore = data["FatJet_globalParT3_QCD"][0].to_numpy() 
            
            met_pt = float(data["PuppiMET_pt"][0])
            met_phi = float(data["PuppiMET_phi"][0])
            
            max_ak4 = min(10, len(ak4_pt))
            max_ak8 = min(3, len(ak8_pt))
            
            ak4_data = np.zeros((10, 8), dtype=np.float32)
            ak4_mask = np.zeros(10, dtype=bool)
            
            for i in range(max_ak4):
                ak4_mask[i] = True
                ak4_data[i, 0] = np.log(ak4_mass[i] + 1.0) 
                ak4_data[i, 1] = np.log(ak4_pt[i] + 1.0) 
                ak4_data[i, 2] = ak4_eta[i] 
                ak4_data[i, 3] = np.sin(ak4_phi[i])
                ak4_data[i, 4] = np.cos(ak4_phi[i])
                ak4_data[i, 5] = float(ak4_isTightBTag[i])
                ak4_data[i, 6] = float(ak4_isMediumBTag[i])
                ak4_data[i, 7] = float(ak4_isLooseBTag[i])
            
            ak8_data = np.zeros((3, 11), dtype=np.float32)
            ak8_mask = np.zeros(3, dtype=bool)
            
            for i in range(max_ak8):
                ak8_mask[i] = True
                ak8_data[i, 0] = np.log(ak8_mass[i] + 1.0)      
                ak8_data[i, 1] = np.log(ak8_pt[i] + 1.0)        
                ak8_data[i, 2] = ak8_eta[i]                     
                ak8_data[i, 3] = np.sin(ak8_phi[i])             
                ak8_data[i, 4] = np.cos(ak8_phi[i])             
                ak8_data[i, 5] = float(ak8_nConstituents[i])    
                ak8_data[i, 6] = ak8_XbbScore[i]                
                ak8_data[i, 7] = ak8_XqqScore[i]                
                ak8_data[i, 8] = ak8_XccScore[i]                
                ak8_data[i, 9] = ak8_XcsScore[i]                
                ak8_data[i, 10] = ak8_XqcdScore[i]              
            
            met_data = np.zeros((1, 3), dtype=np.float32)
            met_data[0, 0] = np.log(met_pt + 1.0)               
            met_data[0, 1] = np.sin(met_phi)                    
            met_data[0, 2] = np.cos(met_phi)                    
            
            inputs = {
                "AK4Jets_data": ak4_data[np.newaxis, :, :],
                "AK4Jets_mask": ak4_mask[np.newaxis, :],
                "AK8Jets_data": ak8_data[np.newaxis, :, :],
                "AK8Jets_mask": ak8_mask[np.newaxis, :],
                "MET_data": met_data[np.newaxis, :, :],
                "MET_mask": np.ones((1, 1), dtype=bool)
            }
            
            cpp_outputs = {
                "h_detection": float(data["spanet_h_detection"][0]),
                "bh_detection": float(data["spanet_bh_detection"][0]),
                "v1_detection": float(data["spanet_v1_detection"][0]),
                "v2_detection": float(data["spanet_v2_detection"][0]),
                "bv1_detection": float(data["spanet_bv1_detection"][0]),
                "bv2_detection": float(data["spanet_bv2_detection"][0]),
                "vbs_detection": float(data["spanet_vbs_detection"][0]),
            }

            return inputs, cpp_outputs
            
    except Exception as e:
        print(f"Error loading event {event_idx} from ROOT file: {e}")
        return None, None

def compare_outputs(cpp_outputs, onnx_outputs, tolerance=0.01):
    output_mapping = {
        "h_detection": 7,
        "bh_detection": 8,
        "v1_detection": 9,
        "v2_detection": 10,
        "bv1_detection": 11,
        "bv2_detection": 12,
        "vbs_detection": 13,
    }
    
    results = {}
    all_match = True
    
    for output_name, onnx_idx in output_mapping.items():
        cpp_value = cpp_outputs[output_name]
        
        onnx_value = float(onnx_outputs[onnx_idx][0])
        
        if abs(cpp_value) > 0:
            rel_diff = abs(cpp_value - onnx_value) / abs(cpp_value)
        else:
            rel_diff = abs(cpp_value - onnx_value)
        
        match = rel_diff <= tolerance
        all_match = all_match and match
        
        results[output_name] = {
            "cpp": cpp_value,
            "onnx": onnx_value,
            "rel_diff": rel_diff,
            "match": match
        }
    
    return results, all_match

def run_inference(session, inputs):
    try:
        outputs = session.run(None, inputs)
        return outputs
    except Exception as e:
        print(f"Error during inference: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Validate ONNX inference against C++ outputs")
    parser.add_argument("--model", "-m", required=True, type=str,
                        help="Path to ONNX model file")
    parser.add_argument("--input", "-i", required=True, type=str,
                        help="Path to input root file")
    parser.add_argument("--max-events", type=int, default=10,
                        help="Maximum number of events to process (default: 10)")
    parser.add_argument("--tolerance", type=float, default=0.01,
                        help="Relative tolerance for comparison (default: 0.01 = 1%)")
    
    args = parser.parse_args()
    
    model_path = Path(args.model)
    if not model_path.exists():
        print(f"Model file not found: {model_path}")
        sys.exit(1)
    
    print(f"Loading model: {model_path}")
    session = load_model(str(model_path))
    
    max_events = args.max_events
    
    print(f"Processing {max_events} events from {args.input}")
    print(f"Using tolerance of {args.tolerance*100:.1f}%")
    
    successful_events = 0
    failed_events = 0
    all_results = []
    output_names = ["h_detection", "bh_detection", "v1_detection", "v2_detection", 
                   "bv1_detection", "bv2_detection", "vbs_detection"]
    
    for event_idx in range(max_events):
        print(f"Processing event {event_idx}/{max_events}")
        
        inputs, cpp_outputs = load_one_event(args.input, event_idx)
        if inputs is None:
            failed_events += 1
            continue
            
        try:
            onnx_outputs = run_inference(session, inputs)
            
            results, all_match = compare_outputs(cpp_outputs, onnx_outputs, args.tolerance)
            all_results.append(results)
            
            if all_match:
                successful_events += 1
            else:
                failed_events += 1
                print(f"Event {event_idx}: Outputs don't match within tolerance!")
                for output_name in output_names:
                    result = results[output_name]
                    if not result["match"]:
                        print(f"  {output_name}: C++={result['cpp']:.6f}, ONNX={result['onnx']:.6f}, "
                              f"rel_diff={result['rel_diff']:.4f}")
                        
        except Exception as e:
            print(f"Error processing event {event_idx}: {e}")
            failed_events += 1
    
    print(f"\n{'='*60}")
    print(f"VALIDATION SUMMARY")
    print(f"{'='*60}")
    print(f"Total events processed: {successful_events + failed_events}")
    print(f"Successful matches: {successful_events}")
    print(f"Failed matches: {failed_events}")
    print(f"Success rate: {successful_events/(successful_events + failed_events)*100:.2f}%")
   
if __name__ == "__main__":
    main()