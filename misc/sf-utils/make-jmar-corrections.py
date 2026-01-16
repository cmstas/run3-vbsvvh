# %%
import correctionlib.schemav2 as cs

# %%
# JMAR Corrections
corr_JMS = cs.Correction(
    name="JMS",
    version=1,
    inputs=[
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
                        cs.CategoryItem(key="nominal", value=-0.3),
                        cs.CategoryItem(key="up", value=0),
                        cs.CategoryItem(key="down", value=-0.6),
                    ]
                ),
            ),
            cs.CategoryItem(
                key="2016postVFP",
                value=cs.Category(
                    nodetype="category",
                    input="sf",
                    content=[
                        cs.CategoryItem(key="nominal", value=-0.24),
                        cs.CategoryItem(key="up", value=0),
                        cs.CategoryItem(key="down", value=-0.48),
                    ]
                ),
            ),
            cs.CategoryItem(
                key="2017",
                value=cs.Category(
                    nodetype="category",
                    input="sf",
                    content=[
                        cs.CategoryItem(key="nominal", value=-0.11),
                        cs.CategoryItem(key="up", value=0),
                        cs.CategoryItem(key="down", value=-0.22),
                    ]
                ),
            ),
            cs.CategoryItem(
                key="2018",
                value=cs.Category(
                    nodetype="category",
                    input="sf",
                    content=[
                        cs.CategoryItem(key="nominal", value=-0.04),
                        cs.CategoryItem(key="up", value=0.06),
                        cs.CategoryItem(key="down", value=-0.14),
                    ]
                ),
            ),
        ]
    )
)

corr_JMR = cs.Correction(
    name="JMR",
    version=1,
    inputs=[
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
                        cs.CategoryItem(key="nominal", value=0.13),
                        cs.CategoryItem(key="up", value=0.43),
                        cs.CategoryItem(key="down", value=0),
                    ]
                ),
            ),
            cs.CategoryItem(
                key="2016postVFP",
                value=cs.Category(
                    nodetype="category",
                    input="sf",
                    content=[
                        cs.CategoryItem(key="nominal", value=0.15),
                        cs.CategoryItem(key="up", value=0.45),
                        cs.CategoryItem(key="down", value=0),
                    ]
                ),
            ),
            cs.CategoryItem(
                key="2017",
                value=cs.Category(
                    nodetype="category",
                    input="sf",
                    content=[
                        cs.CategoryItem(key="nominal", value=0.09),
                        cs.CategoryItem(key="up", value=0.39),
                        cs.CategoryItem(key="down", value=0),
                    ]
                ),
            ),
            cs.CategoryItem(
                key="2018",
                value=cs.Category(
                    nodetype="category",
                    input="sf",
                    content=[
                        cs.CategoryItem(key="nominal", value=0.10),
                        cs.CategoryItem(key="up", value=0.40),
                        cs.CategoryItem(key="down", value=0),
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
            corr_JMS,
            corr_JMR
        ],
    )




#%%   
with open("jmar.json", "w") as fout:
    fout.write(cset.model_dump_json(indent=4, exclude_unset=True))
