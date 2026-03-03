#include "cutflow.h"
#include "tabulate.hpp"

#include <iomanip>
#include <sstream>

std::string &Cutflow::weightCol() {
    static std::string s = "1";
    return s;
}

std::vector<Cutflow::Entry> &Cutflow::entries() {
    static std::vector<Entry> s;
    return s;
}

void Cutflow::SetWeightCol(const std::string &col) {
    weightCol() = col;
}

void Cutflow::Add(RNode df, const std::string &label) {
    Entry e;
    e.label = label;
    if (weightCol() == "1" || weightCol().empty()) {
        e.cResult = df.Count();
        e.isCount = true;
    } else {
        e.wResult = df.Sum<double>(weightCol());
        e.isCount = false;
    }
    entries().push_back(std::move(e));
}

void Cutflow::Clear() {
    entries().clear();
}

std::size_t Cutflow::Size() {
    return entries().size();
}

void Cutflow::Print(std::ostream &out) {
    auto &ents = entries();
    if (ents.empty()) {
        out << "[Cutflow] No entries booked.\n";
        return;
    }

    std::vector<double> vals;
    vals.reserve(ents.size());
    for (auto &e : ents)
        vals.push_back(e.value());

    const double first = vals.front();

    auto fmt = [](double v, int prec = 3) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(prec) << v;
        return oss.str();
    };

    tabulate::Table table;

    table.add_row({"#", "Cut name", "Sum(w)", "rel. eff.", "abs. eff."});

    for (std::size_t i = 0; i < ents.size(); ++i) {
        const double prev   = (i == 0) ? vals[0] : vals[i - 1];
        const double relEff = (prev  > 0.) ? vals[i] / prev  : 0.;
        const double absEff = (first > 0.) ? vals[i] / first : 0.;

        table.add_row({
            std::to_string(i),
            ents[i].label,
            fmt(vals[i]),
            fmt(relEff),
            fmt(absEff)
        });
    }

    table[0].format()
        .font_style({tabulate::FontStyle::bold})
        .font_align(tabulate::FontAlign::center);

    const std::size_t nrows = ents.size() + 1;
    for (std::size_t r = 1; r < nrows; ++r) {
        for (std::size_t c : {0, 2, 3, 4})
            table[r][c].format().font_align(tabulate::FontAlign::right);
    }

    out << "\n[Cutflow] weight column: \"" << weightCol() << "\"\n";
    table.print(out);
    out << "\n";
}
