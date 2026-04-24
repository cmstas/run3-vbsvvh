# Note the keys of each dict should correspond to a key in the SKIM_PATH_DICT

xsec_dict = {

    "data" : {
    },

    "sig_c2v1p0_c3_1p0" : {
        "VBSWWH_OS_c2v1p0_c3_1p0": 0.000643,
        "VBSWWH_SS_c2v1p0_c3_1p0": 0.000316,
        "VBSWZH_c2v1p0_c3_1p0"   : 0.000418,
        "VBSZZH_c2v1p0_c3_1p0"   : 0.000108
    },
    "sig_c2v1p5_c3_1p0" : {
        "VBSWWH_OS_c2v1p5_c3_1p0": 0.001091,
        "VBSWWH_SS_c2v1p5_c3_1p0": 0.000620,
        "VBSWZH_c2v1p5_c3_1p0"   : 0.000707,
        "VBSZZH_c2v1p5_c3_1p0"   : 0.000368 
    },

    "sig_c2v1p0_c3_10p0" : {
        "VBSWWH_OS_c2v1p0_c3_10p0": 0.001587,
        "VBSWWH_SS_c2v1p0_c3_10p0": 0.000445,
        "VBSWZH_c2v1p0_c3_10p0"   : 0.000624,
        "VBSZZH_c2v1p0_c3_10p0"   : 0.000416
    },

    "bkg": {
        #AN v9 2024 (https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=HIG-24-003)
        "QCD_HT100to200_TuneCP5_PSWeights_13TeV"                             : 27849880.0,
        "QCD_HT200to300_TuneCP5_PSWeights_13TeV"                             : 1716997.0,
        "QCD_HT300to500_TuneCP5_PSWeights_13TeV"                             : 351302.0,
        "QCD_HT500to700_TuneCP5_PSWeights_13TeV"                             : 31630.0,
        "QCD_HT700to1000_TuneCP5_PSWeights_13TeV"                            : 6802.0,
        "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV"                           : 1206.0,
        "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV"                           : 98.71,
        "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV"                            : 20.2,
        "TTToSemiLeptonic_TuneCP5_13TeV"                           : 365.34,
        "TTToHadronic_TuneCP5_13TeV"                               : 377.96,
        "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV"           : 19.559,
        "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV"               : 19.559,
        "WJetsToQQ_HT-200to400_TuneCP5_13TeV"                      : 2549.0,
        "WJetsToQQ_HT-400to600_TuneCP5_13TeV"                      : 276.5,
        "WJetsToQQ_HT-600to800_TuneCP5_13TeV"                      : 59.25,
        "WJetsToQQ_HT-800toInf_TuneCP5_13TeV"                      : 28.75,
        "EWKWminus2Jets_WToQQ_dipoleRecoilOn_TuneCP5_13TeV"        : 10.67,
        "EWKWplus2Jets_WToQQ_dipoleRecoilOn_TuneCP5_13TeV"         : 10.67,
        "EWKZ2Jets_ZToLL_M-50_TuneCP5_withDipoleRecoil_13TeV"      : 6.22,
        "EWKZ2Jets_ZToNuNu_M-50_TuneCP5_withDipoleRecoil_13TeV"    : 10.72,
        "EWKZ2Jets_ZToQQ_dipoleRecoilOn_TuneCP5_13TeV"             : 10.67,
        "TTWJetsToQQ_TuneCP5_13TeV"                                : 0.4377,
        "TTWW_TuneCP5_13TeV"                                       : 0.0115,
        "TTWZ_TuneCP5_13TeV"                                       : 0.003884,
        "ttHToNonbb_M125_TuneCP5_13TeV"                            : 0.215,
        "ttHTobb_M125_TuneCP5_13TeV"                               : 0.1279,
        "VHToNonbb_M125_TuneCP5_13TeV"                             : 2.207,
        "WWW_4F_TuneCP5_13TeV"                                     : 0.2086,
        "WWZ_4F_TuneCP5_13TeV"                                     : 0.1651,
        "WZJJ_EWK_InclusivePolarization_TuneCP5_13TeV"             : 0.01701,
        "WZTo1L1Nu2Q_4f_TuneCP5_13TeV"                             : 49.997,
        "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV"                         : 5.6,
        "WZZ_TuneCP5_13TeV"                                        : 0.05565,
        "ZZTo4Q_5f_TuneCP5_13TeV"                                  : 3.451,
        "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV"                         : 3.28,
        "ZZZ_TuneCP5_13TeV"                                        : 0.01398,
        "ZJetsToQQ_HT-200to400_TuneCP5_13TeV"                      : 1012.0,
        "ZJetsToQQ_HT-400to600_TuneCP5_13TeV"                      : 114.2,
        "ZJetsToQQ_HT-600to800_TuneCP5_13TeV"                      : 25.34,
        "ZJetsToQQ_HT-800toInf_TuneCP5_13TeV"                      : 12.99,
        "ZZTo2Nu2Q_5f_TuneCP5_13TeV"                               : 4.58725,
        #-----------------------------------------------------------------------------
        #AN v15 2023 (https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=HIG-24-003)
        "TTTo2L2Nu_TuneCP5_13TeV"                                  : 88.29,
        "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV"    : 80.95,
        "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV"        : 136.02,
        "DYJetsToLL_M-10to50_TuneCP5_13TeV"                        : 20657.0,
        "DYJetsToLL_M-50_TuneCP5_13TeV"                            : 6198.0,
        "EWKWMinus2Jets_WToLNu_M-50_TuneCP5_withDipoleRecoil_13TeV": 32.26,
        "EWKWPlus2Jets_WToLNu_M-50_TuneCP5_withDipoleRecoil_13TeV" : 39.33,
        "TTWJetsToLNu_TuneCP5_13TeV"                               : 0.2043,
        "TTZToLLNuNu_M-10_TuneCP5_13TeV"                           : 0.2529,
        "TTbb_4f_TTTo2L2Nu"                                        : 0.04,
        "TTbb_4f_TTToSemiLeptonic"                                 : 0.62,
        "VBFWH_HToBB_WToLNu_M-125_dipoleRecoilOn_TuneCP5_13TeV"    : 0.02656,
        "WWJJToLNuLNu_EWK_noTop_TuneCP5_13TeV"                     : 0.284,
        "WWTo1L1Nu2Q_4f_TuneCP5_13TeV"                             : 49.997,
        "WWTo2L2Nu_TuneCP5_13TeV"                                  : 12.178,
        "WZTo1L3Nu_4f_TuneCP5_13TeV"                               : 3.05402,
        "WminusH_HToBB_WToLNu_M-125_TuneCP5_13TeV"                 : 0.0490124,
        "WplusH_HToBB_WToLNu_M-125_TuneCP5_13TeV"                  : 0.084876,
        "ZH_HToBB_ZToLL_M-125_TuneCP5_13TeV"                       : 0.0262749,
        "ZZJJTo4L_TuneCP5_13TeV"                                   : 0.00884,
        "ZZTo2L2Nu_TuneCP5_13TeV"                                  : 0.564,
        "ZZTo4L_M-1toInf_TuneCP5_13TeV"                            : 1.256,
        "ggZH_HToBB_ZToLL_M-125_TuneCP5_13TeV"                     : 0.0024614,
        "ZZJJTo4L_EWKnotop_TuneCP5_13TeV"                          : 0.00884, 
        # for numbers below we use averages
        "WJetsToLNu_HT-70To100_TuneCP5_13TeV"                      : 1308.9025,
        "WJetsToLNu_HT-100To200_TuneCP5_13TeV"                     : 1324.85,
        "WJetsToLNu_HT-200To400_TuneCP5_13TeV"                     : 347.935075,
        "WJetsToLNu_HT-400To600_TuneCP5_13TeV"                     : 46.62084375,
        "WJetsToLNu_HT-600To800_TuneCP5_13TeV"                     : 11.23292175,
        "WJetsToLNu_HT-800To1200_TuneCP5_13TeV"                    : 5.07420335,
        "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV"                   : 1.177367725,
        "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV"                    : 0.02426248275,
        #-----------------------------------------------------------------------------
        
        #AN v12 2022 (https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=HIG-24-003)
        "QCD_HT50to100_TuneCP5_PSWeights_13TeV"                    : 187700000.0,
        "ST_s-channel_4f_leptonDecays_TuneCP5_13TeV"               : 3.74,
        "TTbb_4f_TTToHadronic"                                     : 5.5,
        "ttWJets_TuneCP5_13TeV"                                    : 0.4611,
        "ttZJets_TuneCP5_13TeV"                                    : 0.5407,
        "WminusH_HToBB_WToQQ_M-125_TuneCP5_13TeV"                  : 0.3675,
        "WplusH_HToBB_WToQQ_M-125_TuneCP5_13TeV"                   : 0.589,
        "ZH_HToBB_ZToBB_M-125_TuneCP5_13TeV"                       : 0.5612,
        "ZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV"                     : 0.1573,
        "ZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV"                       : 0.5612,
        "ggZH_HToBB_ZToBB_M-125_TuneCP5_13TeV"                     : 0.04319,
        "ggZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV"                   : 0.01222,
        "ggZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV"                     : 0.04319,
        "WZ_TuneCP5_13TeV"                                         : 27.59,
        "WW_TuneCP5_13TeV"                                         : 75.95,
        "ZZ_TuneCP5_13TeV"                                         : 12.17,
        #-----------------------------------------------------------------------------

        #AN-19-004(https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2019/004) (k factor of 1.7 used)
        "GluGluToContinToZZTo2e2mu_TuneCP5_13TeV" : 0.005423,
        "GluGluToContinToZZTo2e2tau_TuneCP5_13TeV" : 0.005423,
        "GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV" : 0.005423,
        "GluGluToContinToZZTo4e_TuneCP5_13TeV" : 0.002703,
        "GluGluToContinToZZTo4mu_TuneCP5_13TeV" : 0.002703,
        "GluGluToContinToZZTo4tau_TuneCP5_13TeV" : 0.002703,
        #-----------------------------------------------------------------------------

        # from https://github.com/TopEFT/topcoffea/blob/main/topcoffea/params/xsec.yml
        "GluGluZH_HToWWTo2L2Nu_TuneCP5_13TeV": 0.00282,
        "GluGluZH_HToWWTo2L2Nu_M-125_TuneCP5_13TeV": 0.00282,
        "HZJ_HToWWTo2L2Nu_ZTo2L_M-125_TuneCP5_13TeV": 0.00177,
        #--------------------------------------------------------------------------
        
        #----------------------------------------------------------------------------
        #AN-19-004 (https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2019/156)
        "tZq_ll_4f_ckm_NLO_TuneCP5_13TeV": 0.0758,
        #----------------------------------------------------------------------------
        
        #from https://github.com/cmstas/VVVNanoLooper/blob/master/weights/xsec.txt
        "SSWW": 0.02794,
        #----------------------------------------------------------------------------
        
        #from TOP-22-008
        "TWZToLL_tlept_Wlept_5f_DR_TuneCP5_13TeV": 0.0015,
        #----------------------------------------------------------------------------
        
        #-----------XSDB--------------------------------------------------------------
        "WJetsToLNu_TuneCP5_13TeV"                                 : 66680.0,
        "WWTo4Q_4f_TuneCP5_13TeV"                                  : 51.03,
        #-----------------------------------------------------------------------------

        # From AN-19-144 (https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2019/144)
        "GluGluHToZZTo4L": 0.011814

    }
}
