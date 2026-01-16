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
                name="flavor",
                type="string",
                description="B/C/L",
            ),
            cs.Variable(
                name="WP",
                type="string",
                description="T/L/N (tight, loose, notag)",
            ),
            cs.Variable(
                name="xaxis",
                type="real",
                description="pt",
            ),
            cs.Variable(
                name="yaxis",
                type="real",
                description="eta",
            ),
        ],
        output=cs.Variable(name="weight", type="real", description="weight"),
        data=cs.Category(
            nodetype="category",
            input="flavor",
            content=[
                cs.CategoryItem(key="B", value=cs.Category(
                    nodetype="category",
                    input="WP",
                    content=[
                        cs.CategoryItem(key="T", value=ttH[0].data),
                        cs.CategoryItem(key="L", value=ttH[1].data),
                        cs.CategoryItem(key="N", value=ttH[2].data),
                    ],
                )),
                cs.CategoryItem(key="C", value=cs.Category(
                    nodetype="category",
                    input="WP",
                    content=[
                        cs.CategoryItem(key="T", value=ttH[3].data),
                        cs.CategoryItem(key="L", value=ttH[4].data),
                        cs.CategoryItem(key="N", value=ttH[5].data),
                    ],
                )),
                cs.CategoryItem(key="L", value=cs.Category(
                    nodetype="category",
                    input="WP",
                    content=[
                        cs.CategoryItem(key="T", value=ttH[6].data),
                        cs.CategoryItem(key="L", value=ttH[7].data),
                        cs.CategoryItem(key="N", value=ttH[8].data),
                    ],
                )),
            ],
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
    corr_hists_2016preVFP = ["hist_beff.root:hist_btight_2016preVFP", 
                             "hist_beff.root:hist_bloose_2016preVFP", 
                             "hist_beff.root:hist_bnotag_2016preVFP",
                             "hist_beff.root:hist_ctight_2016preVFP", 
                             "hist_beff.root:hist_cloose_2016preVFP",
                             "hist_beff.root:hist_cnotag_2016preVFP", 
                             "hist_beff.root:hist_ltight_2016preVFP", 
                             "hist_beff.root:hist_lloose_2016preVFP",
                             "hist_beff.root:hist_lnotag_2016preVFP"]
    
    write_correction(corr_hists_2016preVFP, "Btag Efficiency", "btag_2016preVFP")
    
    corr_hists_2016postVFP = ["hist_beff.root:hist_btight_2016postVFP",
                                "hist_beff.root:hist_bloose_2016postVFP",
                                "hist_beff.root:hist_bnotag_2016postVFP",
                                "hist_beff.root:hist_ctight_2016postVFP",
                                "hist_beff.root:hist_cloose_2016postVFP",
                                "hist_beff.root:hist_cnotag_2016postVFP",
                                "hist_beff.root:hist_ltight_2016postVFP",
                                "hist_beff.root:hist_lloose_2016postVFP",
                                "hist_beff.root:hist_lnotag_2016postVFP"]
    
    write_correction(corr_hists_2016postVFP, "Btag Efficiency", "btag_2016postVFP")

    corr_hists_2017 = ["hist_beff.root:hist_btight_2017",
                        "hist_beff.root:hist_bloose_2017",
                        "hist_beff.root:hist_bnotag_2017",
                        "hist_beff.root:hist_ctight_2017",
                        "hist_beff.root:hist_cloose_2017",
                        "hist_beff.root:hist_cnotag_2017",
                        "hist_beff.root:hist_ltight_2017",
                        "hist_beff.root:hist_lloose_2017",
                        "hist_beff.root:hist_lnotag_2017"]
    
    write_correction(corr_hists_2017, "Btag Efficiency", "btag_2017")

    corr_hists_2018 = ["hist_beff.root:hist_btight_2018",
                        "hist_beff.root:hist_bloose_2018",
                        "hist_beff.root:hist_bnotag_2018",
                        "hist_beff.root:hist_ctight_2018",
                        "hist_beff.root:hist_cloose_2018",
                        "hist_beff.root:hist_cnotag_2018",
                        "hist_beff.root:hist_ltight_2018",
                        "hist_beff.root:hist_lloose_2018",
                        "hist_beff.root:hist_lnotag_2018"]
    
    write_correction(corr_hists_2018, "Btag Efficiency", "btag_2018")