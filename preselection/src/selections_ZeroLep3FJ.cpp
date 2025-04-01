#include "commonSelections.h"
#include "selections_ZeroLep3FJ.h"

namespace ZeroLep3FJ {
    RNode TriggerSelections(RNode df_) {
        return df_.Define("passesTriggers", "HLT_PFHT1050");
    }

    RNode ObjectSelections(RNode df_) {
        return df_;
    }

    RNode EventSelections(RNode df_) {
        return df_.Define("passes_filters_triggers", "passesTriggers")
            .Define("passes_lepton_selection", "passes_filters_triggers && "
                "(nLooseMuons == 0 && nLooseElectrons == 0)");
    }

    RNode runPreselection(RNode df_) {
        auto df = CommonSelections(df_);
        df = ObjectSelections(df);
        df = EventSelections(df);
        return df;
    }
} // OneLep2FJ