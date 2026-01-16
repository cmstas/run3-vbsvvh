# %%
import correctionlib.schemav2 as cs

# %%
# ParticleNet Corrections

# %%
import correctionlib.schemav2 as cs

# %%
# ParticleNet Corrections
corr_W = cs.Correction(
    name="PNET_W",
    version=1,
    inputs=[
        cs.Variable(name="pt", type="real", description="pt"),
        cs.Variable(name="year", type="string", description="year"),
        cs.Variable(name="sf", type="string", description="systematic variation"),
    ],
    output=cs.Variable(name="weight", type="real", description="Multiplicative event weight"),
    data=cs.Category(
        nodetype="category",
        input="year",
        content=[
            cs.CategoryItem(
                key="2016preVFP",
                value=cs.Category(
                    nodetype="category",
                    input="sf",
                    content=[
                        cs.CategoryItem(
                            key="nominal", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[200, 300, 400, 800],
                                content=[0.92, 0.93, 0.94],
                                flow="clamp",
                            ),
                        ),
                        cs.CategoryItem(
                            key="up", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[200, 300, 400, 800],
                                content=[0.96, 0.97, 1.0],
                                flow=1.06,
                            ),
                        ),
                        cs.CategoryItem(
                            key="down", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[200, 300, 400, 800],
                                content=[0.88, 0.89, 0.88],
                                flow=0.82,
                            ),
                        ),
                    ]
                ),
            ),
            cs.CategoryItem(
                key="2016postVFP",
                value=cs.Category(
                    nodetype="category",
                    input="sf",
                    content=[
                        cs.CategoryItem(
                            key="nominal", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[200, 300, 400, 800],
                                content=[0.96, 0.89, 0.85],
                                flow="clamp",
                            ),
                        ),
                        cs.CategoryItem(
                            key="up", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[200, 300, 400, 800],
                                content=[1.01, 0.93, 0.92],
                                flow=0.99,
                            ),
                        ),
                        cs.CategoryItem(
                            key="down", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[200, 300, 400, 800],
                                content=[0.91, 0.85, 0.79],
                                flow=0.73,
                            ),
                        ),
                    ]
                ),
            ),
            cs.CategoryItem(
                key="2017",
                value=cs.Category(
                    nodetype="category",
                    input="sf",
                    content=[
                        cs.CategoryItem(
                            key="nominal", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[200, 300, 400, 800],
                                content=[0.96, 0.94, 0.93],
                                flow="clamp",
                            ),
                        ),
                        cs.CategoryItem(
                            key="up", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[200, 300, 400, 800],
                                content=[0.98, 0.96, 0.97],
                                flow=1.01,
                            ),
                        ),
                        cs.CategoryItem(
                            key="down", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[200, 300, 400, 800],
                                content=[0.94, 0.92, 0.89],
                                flow=0.85,
                            ),
                        ),
                    ]
                ),
            ),
            cs.CategoryItem(
                key="2018",
                value=cs.Category(
                    nodetype="category",
                    input="sf",
                    content=[
                        cs.CategoryItem(
                            key="nominal", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[200, 300, 400, 800],
                                content=[0.91, 0.90, 0.85],
                                flow="clamp",
                            ),
                        ),
                        cs.CategoryItem(
                            key="up", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[200, 300, 400, 800],
                                content=[0.93, 0.92, 0.89],
                                flow=0.93,
                            ),
                        ),
                        cs.CategoryItem(
                            key="down", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[200, 300, 400, 800],
                                content=[0.89, 0.88, 0.81],
                                flow=0.77,
                            ),
                        ),
                    ]
                ),
            ),
        ]
    )
)

corr_H = cs.Correction(
    name="PNET_H",
    version=1,
    inputs=[
        cs.Variable(name="pt", type="real", description="pt"),
        cs.Variable(name="year", type="string", description="year"),
        cs.Variable(name="sf", type="string", description="systematic variation"),
    ],
    output=cs.Variable(name="weight", type="real", description="Multiplicative event weight"),
    data=cs.Category(
        nodetype="category",
        input="year",
        content=[
            cs.CategoryItem(
                key="2016preVFP",
                value=cs.Category(
                    nodetype="category",
                    input="sf",
                    content=[
                        cs.CategoryItem(
                            key="nominal", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[250, 350, 450, 550, 650, 10000],
                                content=[1.131, 1.356, 1.374, 1.376, 1.360],
                                flow="clamp",
                            ),
                        ),
                        cs.CategoryItem(
                            key="up", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[250, 350, 450, 550, 650, 10000],
                                content=[round(i,3) for i in [1.131+0.160, 1.356+0.192, 1.374+0.223, 1.376+0.242, 1.360+0.299]],
                                flow=round(1.360 + 2 * 0.299, 3),
                            ),
                        ),
                        cs.CategoryItem(
                            key="down", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[250, 350, 450, 550, 650, 10000],
                                content=[round(i,3) for i in [1.131-0.138, 1.356-0.148, 1.374-0.206, 1.376-0.216, 1.360-0.257]],
                                flow=round(1.360 - 2 * 0.257, 3),
                            ),
                        ),
                    ]
                ),
            ),
            cs.CategoryItem(
                key="2016postVFP",
                value=cs.Category(
                    nodetype="category",
                    input="sf",
                    content=[
                        cs.CategoryItem(
                            key="nominal", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[250, 350, 450, 550, 650, 10000],
                                content=[1.071, 1.104, 1.163, 1.377, 1.338],
                                flow="clamp",
                            ),
                        ),
                        cs.CategoryItem(
                            key="up", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[250, 350, 450, 550, 650, 10000],
                                content=[round(i,3) for i in [1.071+0.211, 1.104+0.180, 1.163+0.155, 1.377+0.255, 1.338+0.232]],
                                flow=round(1.338 + 2 * 0.232, 3),
                            ),
                        ),
                        cs.CategoryItem(
                            key="down", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[250, 350, 450, 550, 650, 10000],
                                content=[round(i,3) for i in [1.071-0.174, 1.104-0.165, 1.163-0.142, 1.377-0.240, 1.338-0.181]],
                                flow=round(1.338 - 2 * 0.181, 3)
                            ),
                        ),
                    ]
                ),
            ),
            cs.CategoryItem(
                key="2017",
                value=cs.Category(
                    nodetype="category",
                    input="sf",
                    content=[
                        cs.CategoryItem(
                            key="nominal", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[250, 350, 450, 550, 650, 10000],
                                content=[1.213, 1.100, 1.216, 1.192, 1.081],
                                flow="clamp",
                            ),
                        ),
                        cs.CategoryItem(
                            key="up", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[250, 350, 450, 550, 650, 10000],
                                content=[round(i,3) for i in [1.213+0.212, 1.100+0.275, 1.216+0.192, 1.192+0.201, 1.081+0.242]],
                                flow=round(1.081 + 2 * 0.242, 3)
                            ),
                        ),
                        cs.CategoryItem(
                            key="down", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[250, 350, 450, 550, 650, 10000],
                                content=[round(i,3) for i in [1.213-0.201, 1.100-0.271, 1.216-0.178, 1.192-0.188, 1.081-0.237]],
                                flow=round(1.081 - 2 * 0.237, 3)
                            ),
                        ),
                    ]
                ),
            ),
            cs.CategoryItem(
                key="2018",
                value=cs.Category(
                    nodetype="category",
                    input="sf",
                    content=[
                        cs.CategoryItem(
                            key="nominal", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[250, 350, 450, 550, 650, 10000],
                                content=[1.003, 1.172, 1.15, 0.977, 0.998],
                                flow="clamp",
                            ),
                        ),
                        cs.CategoryItem(
                            key="up", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[250, 350, 450, 550, 650, 10000],
                                content=[round(i,3) for i in [1.003+0.067, 1.172+0.180, 1.15+0.121, 0.977+0.117, 0.998+0.168]],
                                flow=round(0.998 + 2 * 0.168, 3),
                            ),
                        ),
                        cs.CategoryItem(
                            key="down", 
                            value=cs.Binning(
                                nodetype="binning",
                                input="pt",
                                edges=[250, 350, 450, 550, 650, 10000],
                                content=[round(i,3) for i in [1.003-0.083, 1.172-0.185, 1.15-0.095, 0.977-0.116, 0.998-0.162]],
                                flow=round(0.998 - 2 * 0.162, 3),
                            ),
                        ),
                    ]
                ),
            ),
        ]
    )
)

cset = cs.CorrectionSet(
        schema_version=2,
        description="PNET corrections MD",
        corrections=[
            corr_W,
            corr_H
        ],
    )

#%%   
with open("pnet.json", "w") as fout:
    fout.write(cset.model_dump_json(indent=4, exclude_unset=True))
