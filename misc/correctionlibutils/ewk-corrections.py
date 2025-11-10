# %%
import uproot
import correctionlib.schemav2 as cs
import correctionlib.convert

# %%
def write_correction(paths, corrname, outname):
    hists = []
    for path in paths:
        with uproot.open(path) as h:
            hists.append(h.to_hist())

    ttH = [correctionlib.convert.from_histogram(hist) for hist in hists]
    
    for corr_ in ttH:
        corr_.data.flow = "clamp"

    corr = cs.Correction(
        name=outname,
        description=corrname,
        version=1,
        inputs=[
            cs.Variable(
                name="xaxis",
                type="real",
                description="pt",
            ),
        ],
        output=cs.Variable(name="weight", type="real", description="weight"),
        data=ttH[0].data
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
    corr_hists_2016preVFP = ["/home/users/aaarora/phys/analysis/vbs-1lep/junk/root_sf/ewk_fix.root:Wgt__pdgid5_quarks_pt_varbin"]
    
    write_correction(corr_hists_2016preVFP, "EWK_Reweighting", "EWK")
    
# %%
