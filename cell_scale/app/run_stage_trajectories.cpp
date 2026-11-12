#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <cxxopts.hpp>

#ifdef BOOTSTRAPPING_USE_OMP
#include <omp.h>
#endif

#include "../src/application/simulation/Simulation.h"
#include "../src/application/bootstrapping/BootstrappingMetrics.h"
#include "../src/domain/adapters/NullLogger.h"
#include "../src/domain/adapters/RandomNoise.h"
#include "../src/domain/cell/CellFactory.h"
#include "../src/domain/gene/GenomeFactory.h"

using application::bootstrapping::makeCellSeed;

class TrackingSimulation final : public application::Simulation {
public:
    explicit TrackingSimulation(int max_t): application::Simulation(max_t) {}
    std::vector<domain::ICell*> liveCells() { return tissue_.getLiveCells(); }
    bool isSaturated() const { return tissue_.isSaturated(); }
};

int main(int argc, char* argv[]) {
    cxxopts::Options options("run_stage_trajectories", "Simulate cell stages over time");
    options.add_options()
        ("n_tissues", "Number of virtual tissues", cxxopts::value<int>()->default_value("100"))
        ("n_cells", "Initial cells per tissue", cxxopts::value<int>()->default_value("1000"))
        ("max_cells", "Carrying capacity (homeostasis)", cxxopts::value<int>()->default_value("1000"))
        ("max_t", "Max years", cxxopts::value<int>()->default_value("80"))
        ("p,division_rate", "Basal division rate", cxxopts::value<double>()->default_value("0.0518"))
        ("beta,brca1_rate", "BRCA1 loss rate", cxxopts::value<double>()->default_value("0.0870"))
        ("tp53_rate", "TP53 loss rate", cxxopts::value<double>()->default_value("0.003"))
        ("delta_low", "Low delta", cxxopts::value<double>()->default_value("0.0787"))
        ("delta_high", "High delta", cxxopts::value<double>()->default_value("0.3549"))
        ("nu,neoplastic_div_rate", "Neoplastic division rate", cxxopts::value<double>()->default_value("0.518"))
        ("kappa_brca1", "Kappa BRCA1", cxxopts::value<double>()->default_value("2.0"))
        ("age_drift", "Age drift", cxxopts::value<double>()->default_value("0.00001"))
        ("d1_threshold", "D1 primer threshold", cxxopts::value<double>()->default_value("2.0"))
        ("d2_threshold", "D2 threshold", cxxopts::value<double>()->default_value("5.0"))
        ("tumor_threshold", "Tumor threshold fraction", cxxopts::value<double>()->default_value("0.05"))
        ("output", "Output CSV path", cxxopts::value<std::string>()->default_value("stage_trajectories.csv"))
        ("seed", "Random seed base", cxxopts::value<int>()->default_value("42"));

    auto result = options.parse(argc, argv);
    int n_tissues = result["n_tissues"].as<int>();
    int n_cells = result["n_cells"].as<int>();
    int max_cells = result["max_cells"].as<int>();
    int max_t = result["max_t"].as<int>();
    double p = result["division_rate"].as<double>();
    double beta = result["brca1_rate"].as<double>();
    double tp53_rate = result["tp53_rate"].as<double>();
    double dl = result["delta_low"].as<double>();
    double dh = result["delta_high"].as<double>();
    double nu = result["neoplastic_div_rate"].as<double>();
    double kappa_brca1 = result["kappa_brca1"].as<double>();
    double age_drift = result["age_drift"].as<double>();
    double d1_th = result["d1_threshold"].as<double>();
    double d2_th = result["d2_threshold"].as<double>();
    std::string out_path = result["output"].as<std::string>();
    int seed_base = result["seed"].as<int>();

    std::ofstream out(out_path);
    out << "tissue_id,year,total_alive,baseline,unstable,unprotected,primer,tumoral,brca1_pm,brca1_mm,tp53_pp,tp53_pm,tp53_mm\n";

    for (int t_idx = 0; t_idx < n_tissues; ++t_idx) {
        int seed = seed_base + t_idx * 1000;
        auto logger = std::make_shared<domain::adapters::NullLogger>();
        TrackingSimulation sim(max_t);
        if (max_cells > 0) {
            sim.setMaxHealthyCells(static_cast<std::size_t>(max_cells));
        }

        // Initialize tissue cells
        for (int i = 0; i < n_cells; ++i) {
            auto genome = domain::genome_factory::makeDefaultGenome(
                {{"BRCA1", beta}, {"TP53", tp53_rate}},
                {{"BRCA1", 0.0}, {"TP53", 0.0}}, logger);
            auto noise = std::make_unique<domain::adapters::RandomNoise>(makeCellSeed(seed, i));
            auto cell = domain::CellFactory::createCustomCell(
                std::move(noise), genome,
                dl, dh,
                p, nu, true, // big_bang = true (neoplastic divides at nu = 10*p)
                d1_th, d2_th, logger,
                0.0, kappa_brca1, age_drift);
            sim.addCell(std::move(cell));
        }

        int last_alive = n_cells;
        int last_baseline = n_cells;
        int last_unstable = 0;
        int last_unprotected = 0;
        int last_primer = 0;
        int last_tumoral = 0;
        int last_b_pm = n_cells, last_b_mm = 0;
        int last_tp_pp = n_cells, last_tp_pm = 0, last_tp_mm = 0;

        for (int year = 1; year <= max_t; ++year) {
            if (!sim.isSaturated()) {
                sim.executeCellCycle();
                auto cells = sim.liveCells();
                last_alive = static_cast<int>(cells.size());

                last_baseline = 0;
                last_unstable = 0;
                last_unprotected = 0;
                last_primer = 0;
                last_tumoral = 0;

                last_b_pm = 0;
                last_b_mm = 0;
                last_tp_pp = 0;
                last_tp_pm = 0;
                last_tp_mm = 0;

                for (auto* c : cells) {
                    auto stage = c->getCurrentCellLifeStage();
                    switch (stage) {
                        case domain::CellLifeStage::BASELINE: last_baseline++; break;
                        case domain::CellLifeStage::UNSTABLE: last_unstable++; break;
                        case domain::CellLifeStage::UNPROTECTED: last_unprotected++; break;
                        case domain::CellLifeStage::PRIMER: last_primer++; break;
                        case domain::CellLifeStage::TUMORAL: last_tumoral++; break;
                        default: break;
                    }

                    auto b = c->getBRCA1Status();
                    if (b == "+/-") last_b_pm++;
                    else if (b == "-/-") last_b_mm++;

                    auto tp = c->getTP53Status();
                    if (tp == "+/+") last_tp_pp++;
                    else if (tp == "+/-") last_tp_pm++;
                    else if (tp == "-/-") last_tp_mm++;
                }
            }
            // Record state for this year (if saturated, carries forward the final expanded tumor state)
            out << t_idx << "," << year << "," << last_alive << ","
                << last_baseline << "," << last_unstable << "," << last_unprotected << "," << last_primer << "," << last_tumoral << ","
                << last_b_pm << "," << last_b_mm << ","
                << last_tp_pp << "," << last_tp_pm << "," << last_tp_mm << "\n";
        }
    }
    out.close();
    std::cout << "Done generating trajectories -> " << out_path << std::endl;
    return 0;
}
