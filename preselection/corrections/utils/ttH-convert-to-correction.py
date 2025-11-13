# %%
import uproot
import correctionlib.schemav2 as cs
import correctionlib.convert

# %%
def write_correction(paths, corrname, outname):
    tth = []
    for path in paths:
        with uproot.open(path) as h:
            tth.append(correctionlib.convert.from_histogram(h.to_hist()))
            tth.append(correctionlib.convert.from_histogram(h.to_hist()))
            tth[-1].data.content = list(h.values().flatten() + h.errors().flatten())
            tth.append(correctionlib.convert.from_histogram(h.to_hist()))
            tth[-1].data.content = list(h.values().flatten() - h.errors().flatten())

    for corr_ in tth:
        corr_.data.flow = "clamp"

    corr = cs.Correction(
        name=outname,
        description=corrname,
        version=1,
        inputs=[
            cs.Variable(
                name="variation",
                type="string",
                description="variation"
            ),
            cs.Variable(
                name="year",
                type="string",
                description="year",
            ),
            cs.Variable(
                name="xaxis",
                type="real",
                description="eta",
            ),
            cs.Variable(
                name="yaxis",
                type="real",
                description="pt",
            ),
        ],
        output=cs.Variable(name="weight", type="real", description="weight"),
        data=cs.Category(
            nodetype="category",
            input="year",
            content=[
                cs.CategoryItem(
                    key="2016preVFP",
                    value=cs.Category( 
                        nodetype="category",
                        input="variation",
                        content=[
                            cs.CategoryItem(key="nominal", value=tth[0].data),
                            cs.CategoryItem(key="up", value=tth[1].data),
                            cs.CategoryItem(key="down", value=tth[2].data),
                        ],
                    )
                ),
                cs.CategoryItem(
                    key="2016postVFP",
                    value=cs.Category( 
                        nodetype="category",
                        input="variation",
                        content=[
                            cs.CategoryItem(key="nominal", value=tth[3].data),
                            cs.CategoryItem(key="up", value=tth[4].data),
                            cs.CategoryItem(key="down", value=tth[5].data),
                        ],
                    )
                ),
                cs.CategoryItem(
                    key="2017",
                    value=cs.Category( 
                        nodetype="category",
                        input="variation",
                        content=[
                            cs.CategoryItem(key="nominal", value=tth[6].data),
                            cs.CategoryItem(key="up", value=tth[7].data),
                            cs.CategoryItem(key="down", value=tth[8].data),
                        ],
                    )
                ),
                cs.CategoryItem(
                    key="2018",
                    value=cs.Category( 
                        nodetype="category",
                        input="variation",
                        content=[
                            cs.CategoryItem(key="nominal", value=tth[9].data),
                            cs.CategoryItem(key="up", value=tth[10].data),
                            cs.CategoryItem(key="down", value=tth[11].data),
                        ],
                    )
                ),
            ]
        )
    )

    cset2 = cs.CorrectionSet(
        schema_version=2,
        description=corrname,
        corrections=[
            corr,
        ],
    )
    
    with open(f"{outname}.json", "w") as fout:
        fout.write(cset2.model_dump_json(indent=4, exclude_unset=True))

# %%
if __name__ == "__main__":
    # #muon
    # corr_hists = ["/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/muon/egammaEffi2016APV_EGM2D.root:EGamma_SF2D", "/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/muon/egammaEffi2016_EGM2D.root:EGamma_SF2D", "/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/muon/egammaEffi2017_EGM2D.root:EGamma_SF2D", "/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/muon/egammaEffi2018_EGM2D.root:EGamma_SF2D"]
    # write_correction(corr_hists, "ttH SFs", "ttH_MuonID_SF")

    # corr_hists = ["/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/muon/egammaEffi2016APV_iso_EGM2D.root:EGamma_SF2D", "/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/muon/egammaEffi2016_iso_EGM2D.root:EGamma_SF2D", "/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/muon/egammaEffi2017_iso_EGM2D.root:EGamma_SF2D", "/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/muon/egammaEffi2018_iso_EGM2D.root:EGamma_SF2D"]
    # write_correction(corr_hists, "ttH SFs", "ttH_MuonISO_SF")

    # # electron
    # corr_hists = ["/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/elec/egammaEffi2016APV_2lss_EGM2D.root:EGamma_SF2D", "/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/elec/egammaEffi2016_2lss_EGM2D.root:EGamma_SF2D", "/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/elec/egammaEffi2017_2lss_EGM2D.root:EGamma_SF2D", "/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/elec/egammaEffi2018_2lss_EGM2D.root:EGamma_SF2D"]
    # write_correction(corr_hists, "ttH SFs", "ttH_ElectronID_SF")

    # corr_hists = ["/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/elec/egammaEffi2016APV_iso_EGM2D.root:EGamma_SF2D", "/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/elec/egammaEffi2016_iso_EGM2D.root:EGamma_SF2D", "/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/elec/egammaEffi2017_iso_EGM2D.root:EGamma_SF2D", "/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/elec/egammaEffi2018_iso_EGM2D.root:EGamma_SF2D"]
    # write_correction(corr_hists, "ttH SFs", "ttH_ElectronISO_SF")

    corr_hists = ["/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/elec/UL2016_preVFP.root:EGamma_SF2D", "/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/elec/UL2016_postVFP.root:EGamma_SF2D", "/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/elec/UL2017.root:EGamma_SF2D", "/home/users/aaarora/phys/vbs/vbs-1lep-new/junk/root_sf/elec/UL2018.root:EGamma_SF2D"]
    write_correction(corr_hists, "Electron Trigger SFs", "trigger")
# %%
