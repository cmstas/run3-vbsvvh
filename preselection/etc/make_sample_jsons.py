import os
import xsec_ref


# Dictionary of the skim sets and the analysis channels they serve
# NOTE: 2lepSS has no skim?
SKIM_INFO_DICT = {
    "0lep_0FJ" : {
        "ana_chans": ["0lep_0FJ_6j"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_0Lep0FJ_11Feb2026_v4",
    },
    "0lep_1FJ" : {
        "ana_chans": ["0lep_1FJ_met", "0lep_1FJ_4j"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_0Lep1FJ_11Feb2026_v4",
    },
    "0lep_2FJ" : {
        "ana_chans": ["0lep_2FJ_met", "0lep_2FJ_2j"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_0Lep2FJ_11Feb2026_v4",
    },
    "0lep_3FJ" : {
        "ana_chans": ["0lep_3FJ"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_0Lep3FJ_11Feb2026_v4",
    },
    "1lep_1FJ" : {
        "ana_chans": ["1lep_1FJ","1lep_2FJ"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_1Lep1FJ_11Feb2026_v4",
    },
    "2lep_1FJ" : {
        "ana_chans": ["2lepOSOF_1FJ", "2lepOSSF_1FJ"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_2Lep1FJ_11Feb2026_v4",
    },
    "2lep_2FJ" : {
        "ana_chans": ["2lepOSSF_2FJ"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_2Lep2FJ_11Feb2026_v4",
    },
    "3lep"    : {
        "ana_chans": ["3lep"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_3Lep_11Feb2026_v4",
    },
    "4lep"    : {
        "ana_chans": ["4lep"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_4Lep_11Feb2026_v4",
    },
}

SAMPLE_NAMES_LST = [
    "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v2_NANOAODSIM",
    "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v2_NANOAODSIM",
    "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v5_NANOAODSIM",
    "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v2_NANOAODSIM",
    "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v2_NANOAODSIM",
    "EWKWMinus2Jets_WToLNu_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "EWKWMinus2Jets_WToLNu_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "EWKWMinus2Jets_WToLNu_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "EWKWMinus2Jets_WToLNu_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "EWKWPlus2Jets_WToLNu_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "EWKWPlus2Jets_WToLNu_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "EWKWPlus2Jets_WToLNu_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "EWKWPlus2Jets_WToLNu_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "EWKWminus2Jets_WToQQ_dipoleRecoilOn_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "EWKWminus2Jets_WToQQ_dipoleRecoilOn_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "EWKWminus2Jets_WToQQ_dipoleRecoilOn_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "EWKWminus2Jets_WToQQ_dipoleRecoilOn_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "EWKWplus2Jets_WToQQ_dipoleRecoilOn_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "EWKWplus2Jets_WToQQ_dipoleRecoilOn_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "EWKWplus2Jets_WToQQ_dipoleRecoilOn_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "EWKWplus2Jets_WToQQ_dipoleRecoilOn_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "EWKZ2Jets_ZToLL_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "EWKZ2Jets_ZToLL_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "EWKZ2Jets_ZToLL_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "EWKZ2Jets_ZToLL_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "EWKZ2Jets_ZToNuNu_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "EWKZ2Jets_ZToNuNu_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "EWKZ2Jets_ZToNuNu_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "EWKZ2Jets_ZToNuNu_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "EWKZ2Jets_ZToQQ_dipoleRecoilOn_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v2_NANOAODSIM",
    "EWKZ2Jets_ZToQQ_dipoleRecoilOn_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "EWKZ2Jets_ZToQQ_dipoleRecoilOn_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "EWKZ2Jets_ZToQQ_dipoleRecoilOn_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v2_NANOAODSIM",
    "GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "GluGluZH_HToWWTo2L2Nu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "GluGluZH_HToWWTo2L2Nu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "HZJ_HToWWTo2L2Nu_ZTo2L_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "HZJ_HToWWTo2L2Nu_ZTo2L_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "SSWW_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "SSWW_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "SSWW_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "SSWW_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v2_NANOAODSIM",
    "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v2_NANOAODSIM",
    "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v2_NANOAODSIM",
    "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v2_NANOAODSIM",
    "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v2_NANOAODSIM",
    "TTToHadronic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v2_NANOAODSIM",
    "TTToHadronic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v2_NANOAODSIM",
    "TTToHadronic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v2_NANOAODSIM",
    "TTToHadronic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v2_NANOAODSIM",
    "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v2_NANOAODSIM",
    "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v2_NANOAODSIM",
    "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v2_NANOAODSIM",
    "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v2_NANOAODSIM",
    "TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "TTWW_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "TTWW_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "TTWW_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "TTWW_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "TTWZ_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "TTWZ_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "TTWZ_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "TTWZ_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "TWZToLL_tlept_Wlept_5f_DR_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "TWZToLL_tlept_Wlept_5f_DR_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "TWZToLL_tlept_Wlept_5f_DR_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "TWZToLL_tlept_Wlept_5f_DR_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "VBFWH_HToBB_WToLNu_M-125_dipoleRecoilOn_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "VBFWH_HToBB_WToLNu_M-125_dipoleRecoilOn_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "VBFWH_HToBB_WToLNu_M-125_dipoleRecoilOn_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "VBFWH_HToBB_WToLNu_M-125_dipoleRecoilOn_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "VHToNonbb_M125_TuneCP5_13TeV-amcatnloFXFX_madspin_pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "VHToNonbb_M125_TuneCP5_13TeV-amcatnloFXFX_madspin_pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "VHToNonbb_M125_TuneCP5_13TeV-amcatnloFXFX_madspin_pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "VHToNonbb_M125_TuneCP5_13TeV-amcatnloFXFX_madspin_pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1_ext1-v1_NANOAODSIM",
    "WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1_ext1-v1_NANOAODSIM",
    "WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1_ext1-v1_NANOAODSIM",
    "WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1_ext1-v1_NANOAODSIM",
    "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1_ext1-v1_NANOAODSIM",
    "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1_ext1-v1_NANOAODSIM",
    "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1_ext1-v1_NANOAODSIM",
    "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1_ext1-v1_NANOAODSIM",
    "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1_ext2-v1_NANOAODSIM",
    "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1_ext2-v1_NANOAODSIM",
    "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1_ext2-v1_NANOAODSIM",
    "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1_ext2-v1_NANOAODSIM",
    "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1_ext1-v1_NANOAODSIM",
    "WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1_ext1-v1_NANOAODSIM",
    "WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1_ext1-v1_NANOAODSIM",
    "WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1_ext1-v1_NANOAODSIM",
    "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WWJJToLNuLNu_EWK_noTop_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v2_NANOAODSIM",
    "WWJJToLNuLNu_EWK_noTop_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WWJJToLNuLNu_EWK_noTop_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WWJJToLNuLNu_EWK_noTop_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WWTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WWTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WWTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WWTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WWTo4Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WWTo4Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WWTo4Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WWTo4Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1_ext1-v1_NANOAODSIM",
    "WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1_ext1-v1_NANOAODSIM",
    "WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1_ext1-v1_NANOAODSIM",
    "WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1_ext1-v1_NANOAODSIM",
    "WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1_ext1-v1_NANOAODSIM",
    "WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1_ext1-v1_NANOAODSIM",
    "WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1_ext1-v1_NANOAODSIM",
    "WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1_ext1-v1_NANOAODSIM",
    "WW_TuneCP5_13TeV-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WW_TuneCP5_13TeV-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WW_TuneCP5_13TeV-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WW_TuneCP5_13TeV-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WZJJ_EWK_InclusivePolarization_TuneCP5_13TeV_madgraph-madspin-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WZJJ_EWK_InclusivePolarization_TuneCP5_13TeV_madgraph-madspin-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WZJJ_EWK_InclusivePolarization_TuneCP5_13TeV_madgraph-madspin-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WZJJ_EWK_InclusivePolarization_TuneCP5_13TeV_madgraph-madspin-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WZTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WZTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WZTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WZTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WZTo1L3Nu_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WZTo1L3Nu_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WZTo1L3Nu_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WZTo1L3Nu_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1_ext1-v1_NANOAODSIM",
    "WZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1_ext1-v1_NANOAODSIM",
    "WZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1_ext1-v1_NANOAODSIM",
    "WZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1_ext1-v1_NANOAODSIM",
    "WZ_TuneCP5_13TeV-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WZ_TuneCP5_13TeV-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WZ_TuneCP5_13TeV-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WZ_TuneCP5_13TeV-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WminusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WminusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WminusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WminusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WminusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v2_NANOAODSIM",
    "WminusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WminusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WminusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WplusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WplusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WplusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WplusH_HToBB_WToLNu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "WplusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "WplusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "WplusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "WplusH_HToBB_WToQQ_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v2_NANOAODSIM",
    "ZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v2_NANOAODSIM",
    "ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ZZTo2Nu2Q_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ZZTo2Nu2Q_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ZZTo2Nu2Q_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ZZTo2Nu2Q_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ZZTo4L_M-1toInf_TuneCP5_13TeV_powheg_pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ZZTo4L_M-1toInf_TuneCP5_13TeV_powheg_pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ZZTo4L_M-1toInf_TuneCP5_13TeV_powheg_pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ZZTo4L_M-1toInf_TuneCP5_13TeV_powheg_pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ZZTo4Q_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ZZTo4Q_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ZZTo4Q_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ZZTo4Q_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1_ext1-v1_NANOAODSIM",
    "ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1_ext1-v1_NANOAODSIM",
    "ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1_ext1-v1_NANOAODSIM",
    "ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1_ext1-v1_NANOAODSIM",
    "ZZ_TuneCP5_13TeV-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ZZ_TuneCP5_13TeV-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ZZ_TuneCP5_13TeV-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ZZ_TuneCP5_13TeV-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ggZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ggZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ggZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ggZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ggZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ggZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ggZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ggZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ggZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ggZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ggZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ggZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ggZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v2_NANOAODSIM",
    "ggZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ggZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ggZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv15-150X_mcRun2_asymptotic_preVFP_v1-v1_NANOAODSIM",
    "ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
    "ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8_RunIISummer20UL16NanoAODv15-150X_mcRun2_asymptotic_v1-v1_NANOAODSIM",
    "ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8_RunIISummer20UL17NanoAODv15-150X_mc2017_realistic_v1-v1_NANOAODSIM",
    "ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8_RunIISummer20UL18NanoAODv15-150X_mc2018_realistic_v1-v1_NANOAODSIM",
]

# Get the list of root files in a given dir
#     - Input should be a full path to the given dir
#     - Return the list of root files joined with the full path
def get_root_file_lst(path_to_dir):
    out_lst = []
    file_name_lst = os.listdir(path_to_dir)
    for fname in file_name_lst:
        if not fname.endswith(".root"): continue
        fname_fullpath = os.path.join(path_to_dir,fname)
        out_lst.append(fname_fullpath)
    return out_lst


# Match a sample name to an xsec in the given list
# Return the name of the xsec and the value
def match_xsec(dataset_name,xsec_lst):
    n_matches = 0
    for xsec_name in xsec_lst:
        if dataset_name.startswith(xsec_name):
                dataset_name_short = xsec_name
                dataset_xsec = xsec_dict[xsec_name]
                match_bool = True
                n_matches += n_matches

    # Throw and error if we did not find a match
    if n_matches == 0:
        raise Exception(f"Failed to find xsec match for sample: \"{dataset_name}\"")
    if n_matches > 1:
        raise Exception(f"Found more than one matching xsec for sample: \"{dataset_name}\"")

    return [dataset_name_short,dataset_xsec]


# Create a json for a given sample
#     - Input is the full path to the sample in question
def make_json_for_dataset(path_to_dataset,xsec_lst):

        # Grab dataset_name from /full/path/to/dataset_name
        dataset_name = path_to_dataset.split("/")[-1]

        # Find the xsec for this dataset
        dataset_name_short, xsec_val = match_xsec(dataset_name,xsec_lst)
        print(dataset_name_short,xsec_val)

        # Get full paths to all root files for this dataset
        file_fullpath_lst = get_root_file_lst(path_to_dataset)
        #print(file_fullpath_lst)


def main():

    skim_fulpath_dict = {}
    xsec_lst = xsec_ref.xsec_dict["bkg_run2"]
    print(xsec_lst)
    exit()

    # Loop over the skim sets
    for skim_set_name in SKIM_INFO_DICT:
        path = SKIM_INFO_DICT[skim_set_name]["path"]
        skim_fulpath_dict[skim_set_name] = {}

        # Get the modified tag at the end of the skim name, remove this when skim formatting is updated
        # NOTE: Remove this when skim formatting is updated in later versions
        tag_tmp = path.split("/")[-1] 
        tag_tmp = tag_tmp[:-1] + "2"

        # Loop over the set of datasets
        for dataset_name in SAMPLE_NAMES_LST:
            skim_fulpath_dict[skim_set_name][dataset_name] = []

            # Modify the dataset name with the second tag
            # NOTE Remove this when skim formatting is updated in later versions
            dataset_name_tmp = f"{dataset_name}_{tag_tmp}"

            # Get the full path to this dataset, and create the json
            path_full = os.path.join(path,dataset_name_tmp)
            make_json_for_dataset(path_full,xsec_lst)



    #import json
    #with open('tmp.json', 'w') as fp:
        #json.dump(skim_fulpath_dict, fp, indent=4)



main()
