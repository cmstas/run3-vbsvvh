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

def load_batch_events(input_file, start_idx=0, batch_size=10):
    try:
        with uproot.open(input_file + ":Events") as tree:
            data = tree.arrays(library="ak", entry_start=start_idx, entry_stop=start_idx+batch_size)
            
            batch_inputs = []
            batch_cpp_outputs = []
            
            for i in range(len(data["Jet_pt"])):
                ak4_arrays = {
                    'pt': data["Jet_pt"][i].to_numpy(),
                    'eta': data["Jet_eta"][i].to_numpy(),
                    'phi': data["Jet_phi"][i].to_numpy(),
                    'mass': data["Jet_mass"][i].to_numpy(),
                    'tight_btag': data["Jet_isTightBTag"][i].to_numpy().astype(np.float32),
                    'medium_btag': data["Jet_isMediumBTag"][i].to_numpy().astype(np.float32),
                    'loose_btag': data["Jet_isLooseBTag"][i].to_numpy().astype(np.float32)
                }
                
                ak8_arrays = {
                    'pt': data["FatJet_pt"][i].to_numpy(),
                    'eta': data["FatJet_eta"][i].to_numpy(),
                    'phi': data["FatJet_phi"][i].to_numpy(),
                    'mass': data["FatJet_mass"][i].to_numpy(),
                    'nConstituents': data["FatJet_nConstituents"][i].to_numpy().astype(np.float32),
                    'xbb': data["FatJet_globalParT3_Xbb"][i].to_numpy(),
                    'xqq': data["FatJet_globalParT3_Xqq"][i].to_numpy(),
                    'xcc': data["FatJet_globalParT3_Xcc"][i].to_numpy(),
                    'xcs': data["FatJet_globalParT3_Xcs"][i].to_numpy(),
                    'xqcd': data["FatJet_globalParT3_QCD"][i].to_numpy()
                }
                
                met_pt = float(data["PuppiMET_pt"][i])
                met_phi = float(data["PuppiMET_phi"][i])

                lep_arrays = {
                    'pt': data["Lepton_pt"][i].to_numpy(),
                    'eta': data["Lepton_eta"][i].to_numpy(),
                    'phi': data["Lepton_phi"][i].to_numpy(),
                    'mass': data["Lepton_mass"][i].to_numpy(),
                    'charge': data["Lepton_charge"][i].to_numpy().astype(np.float32)
                }
                
                n_leptons = len(lep_arrays['pt'])
                lepton_data = np.zeros((2, 5))
                if n_leptons > 0:
                    lepton_data[0] = [lep_arrays['pt'][0], lep_arrays['eta'][0], lep_arrays['phi'][0], 
                                     lep_arrays['mass'][0], lep_arrays['charge'][0]]
                if n_leptons > 1:
                    lepton_data[1] = [lep_arrays['pt'][1], lep_arrays['eta'][1], lep_arrays['phi'][1],
                                     lep_arrays['mass'][1], lep_arrays['charge'][1]]

                max_ak4 = min(10, len(ak4_arrays['pt']))
                ak4_data = np.zeros((10, 8), dtype=np.float32)
                ak4_mask = np.zeros(10, dtype=bool)
                
                if max_ak4 > 0:
                    ak4_mask[:max_ak4] = True
                    ak4_data[:max_ak4, 0] = np.log(ak4_arrays['mass'][:max_ak4] + 1.0)
                    ak4_data[:max_ak4, 1] = np.log(ak4_arrays['pt'][:max_ak4] + 1.0)
                    ak4_data[:max_ak4, 2] = ak4_arrays['eta'][:max_ak4]
                    ak4_data[:max_ak4, 3] = np.sin(ak4_arrays['phi'][:max_ak4])
                    ak4_data[:max_ak4, 4] = np.cos(ak4_arrays['phi'][:max_ak4])
                    ak4_data[:max_ak4, 5] = ak4_arrays['tight_btag'][:max_ak4]
                    ak4_data[:max_ak4, 6] = ak4_arrays['medium_btag'][:max_ak4]
                    ak4_data[:max_ak4, 7] = ak4_arrays['loose_btag'][:max_ak4]
                
                max_ak8 = min(3, len(ak8_arrays['pt']))
                ak8_data = np.zeros((3, 11), dtype=np.float32)
                ak8_mask = np.zeros(3, dtype=bool)
                
                if max_ak8 > 0:
                    ak8_mask[:max_ak8] = True
                    ak8_data[:max_ak8, 0] = np.log(ak8_arrays['mass'][:max_ak8] + 1.0)
                    ak8_data[:max_ak8, 1] = np.log(ak8_arrays['pt'][:max_ak8] + 1.0)
                    ak8_data[:max_ak8, 2] = ak8_arrays['eta'][:max_ak8]
                    ak8_data[:max_ak8, 3] = np.sin(ak8_arrays['phi'][:max_ak8])
                    ak8_data[:max_ak8, 4] = np.cos(ak8_arrays['phi'][:max_ak8])
                    ak8_data[:max_ak8, 5] = ak8_arrays['nConstituents'][:max_ak8]
                    ak8_data[:max_ak8, 6] = ak8_arrays['xbb'][:max_ak8]
                    ak8_data[:max_ak8, 7] = ak8_arrays['xqq'][:max_ak8]
                    ak8_data[:max_ak8, 8] = ak8_arrays['xcc'][:max_ak8]
                    ak8_data[:max_ak8, 9] = ak8_arrays['xcs'][:max_ak8]
                    ak8_data[:max_ak8, 10] = ak8_arrays['xqcd'][:max_ak8]
                
                met_data = np.array([[np.log(met_pt + 1.0), np.sin(met_phi), np.cos(met_phi)]], dtype=np.float32)
                
                lep1_data = np.zeros((1, 6), dtype=np.float32)
                lep1_mask = np.array([True], dtype=bool)
                lep1_data[0] = [np.log(lepton_data[0, 0] + 1.0), lepton_data[0, 1], 
                               np.sin(lepton_data[0, 2]), np.cos(lepton_data[0, 2]),
                               np.log(lepton_data[0, 3] + 1.0), lepton_data[0, 4]]

                lep2_data = np.zeros((1, 6), dtype=np.float32)
                lep2_mask = np.array([True], dtype=bool)
                lep2_data[0] = [np.log(lepton_data[1, 0] + 1.0), lepton_data[1, 1],
                               np.sin(lepton_data[1, 2]), np.cos(lepton_data[1, 2]),
                               np.log(lepton_data[1, 3] + 1.0), lepton_data[1, 4]]         

                inputs = {
                    "AK4Jets_data": ak4_data[np.newaxis, :, :],
                    "AK4Jets_mask": ak4_mask[np.newaxis, :],
                    "AK8Jets_data": ak8_data[np.newaxis, :, :],
                    "AK8Jets_mask": ak8_mask[np.newaxis, :],
                    "MET_data": met_data[np.newaxis, :, :],
                    "MET_mask": np.ones((1, 1), dtype=bool),
                    "Lepton1_data": lep1_data[np.newaxis, :, :],
                    "Lepton1_mask": lep1_mask[np.newaxis, :],
                    "Lepton2_data": lep2_data[np.newaxis, :, :],
                    "Lepton2_mask": lep2_mask[np.newaxis, :],
                }
                
                cpp_outputs = {
                    "h_detection": float(data["spanet_h_detection"][i]),
                    "bh_detection": float(data["spanet_bh_detection"][i]),
                    "v1_detection": float(data["spanet_v1_detection"][i]),
                    "v2_detection": float(data["spanet_v2_detection"][i]),
                    "bv1_detection": float(data["spanet_bv1_detection"][i]),
                    "bv2_detection": float(data["spanet_bv2_detection"][i]),
                    "vbs_detection": float(data["spanet_vbs_detection"][i]),
                }
                
                batch_inputs.append(inputs)
                batch_cpp_outputs.append(cpp_outputs)
            
            return batch_inputs, batch_cpp_outputs
            
    except Exception as e:
        print(f"Error loading batch starting at event {start_idx}: {e}")
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
    parser.add_argument("--max-events", type=int, default=-1,
                        help="Maximum number of events to process (default: -1 = all events)")
    parser.add_argument("--batch-size", type=int, default=10,
                        help="Batch size for processing events (default: 10)")
    parser.add_argument("--tolerance", type=float, default=0.01,
                        help="Relative tolerance for comparison (default: 0.01 = 1%)")
    
    args = parser.parse_args()
    
    model_path = Path(args.model)
    if not model_path.exists():
        print(f"Model file not found: {model_path}")
        sys.exit(1)
    
    print(f"Loading model: {model_path}")
    session = load_model(str(model_path))
    
    max_events = args.max_events if args.max_events > 0 else float('inf')
    batch_size = min(args.batch_size, max_events)
    
    print(f"Processing {max_events} events from {args.input} in batches of {batch_size}")
    print(f"Using tolerance of {args.tolerance*100:.1f}%")
    
    successful_events = 0
    failed_events = 0
    all_results = []
    output_names = ["h_detection", "bh_detection", "v1_detection", "v2_detection", 
                   "bv1_detection", "bv2_detection", "vbs_detection"]
    
    for batch_start in range(0, max_events, batch_size):
        batch_end = min(batch_start + batch_size, max_events)
        current_batch_size = batch_end - batch_start
        
        print(f"Processing batch {batch_start//batch_size + 1}: events {batch_start}-{batch_end-1}")
        
        batch_inputs, batch_cpp_outputs = load_batch_events(args.input, batch_start, current_batch_size)
        if batch_inputs is None:
            failed_events += current_batch_size
            continue
        
        for i, (inputs, cpp_outputs) in enumerate(zip(batch_inputs, batch_cpp_outputs)):
            event_idx = batch_start + i
            
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