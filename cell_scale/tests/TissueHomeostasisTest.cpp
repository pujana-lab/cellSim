#include <gtest/gtest.h>
#include "domain/tissue/Tissue.h"
#include "domain/cell/CellFactory.h"
#include "domain/gene/GenomeFactory.h"
#include "domain/adapters/FixedNoise.h"
#include "domain/adapters/NullLogger.h"

using namespace domain;

class TissueHomeostasisTest : public ::testing::Test {
protected:
    std::shared_ptr<adapters::NullLogger> logger_ = std::make_shared<adapters::NullLogger>();
};

TEST_F(TissueHomeostasisTest, HealthyCellsDoNotExceedCarryingCapacityWithHighDivisionRate) {
    Tissue tissue(logger_);
    tissue.setMaxHealthyCells(10);

    // Add 10 initial healthy cells with 100% division rate (division_rate = 1.0)
    // and no mutations
    for (int i = 0; i < 10; ++i) {
        Genome genome = genome_factory::makeDefaultGenome(
            {{"BRCA1", 0.0}, {"TP53", 0.0}},
            {{"BRCA1", 0.0}, {"TP53", 0.0}},
            logger_
        );
        auto noise = std::make_unique<adapters::FixedNoise>(CellNoise{0.0}); // always triggers division
        auto cell = CellFactory::createCustomCell(
            std::move(noise), genome,
            0.0, 0.0,
            1.0,   // division_rate = 1.0
            0.0,   // neoplastic_division_rate
            false, // big_bang
            2.0, 5.0,
            logger_
        );
        tissue.addCell(std::move(cell));
    }

    ASSERT_EQ(tissue.size(), 10u);

    // Run multiple cycles with 100% division rate
    for (int cycle = 0; cycle < 5; ++cycle) {
        tissue.live();
        // Since capacity is 10, healthy cell count must not exceed 10
        EXPECT_EQ(tissue.size(), 10u);
        EXPECT_EQ(tissue.lastBirthCount(), 0);
    }
}

TEST_F(TissueHomeostasisTest, TumoralCellsCanDivideAndExpandBeyondCapacity) {
    Tissue tissue(logger_);
    tissue.setMaxHealthyCells(5);

    // Add 5 initial healthy cells with division_rate = 0.0
    for (int i = 0; i < 5; ++i) {
        Genome genome = genome_factory::makeDefaultGenome(
            {{"BRCA1", 0.0}, {"TP53", 0.0}},
            {{"BRCA1", 0.0}, {"TP53", 0.0}},
            logger_
        );
        auto noise = std::make_unique<adapters::FixedNoise>(CellNoise{0.5});
        auto cell = CellFactory::createCustomCell(
            std::move(noise), genome,
            0.0, 0.0,
            0.0, 0.0, false,
            2.0, 5.0,
            logger_
        );
        tissue.addCell(std::move(cell));
    }

    // Add 1 cell that transforms into neoplastic via high delta instability (D1>2, D2>5)
    {
        auto genome = genome_factory::makeDefaultGenome(
            {{"BRCA1", 0.0}, {"TP53", 0.0}},
            {{"BRCA1", 0.0}, {"TP53", 0.0}},
            logger_
        );
        genome.mutate("TP53");
        genome.mutate("TP53"); // TP53 -/-

        auto noise = std::make_unique<adapters::FixedNoise>(CellNoise{0.0}); // triggers division
        auto neo_cell = CellFactory::createCustomCell(
            std::move(noise), genome,
            5.0, 5.0, // jump D1 > 2.0 and D2 > 5.0 in 1 cycle
            0.0,   // normal division rate
            1.0,   // neoplastic division rate = 1.0
            true,  // big bang mode enabled
            2.0, 5.0,
            logger_
        );
        tissue.addCell(std::move(neo_cell));
    }

    EXPECT_EQ(tissue.size(), 6u); // 5 healthy + 1 candidate neoplastic

    // Cycle 1: neo_cell accumulates D1>2 and D2>5, becomes neoplastic
    tissue.live();

    // Cycle 2: neoplastic cell divides with 100% neoplastic rate despite size >= 5
    tissue.live();

    // Tumoral cells are not blocked by max_healthy_cells
    EXPECT_GT(tissue.size(), 6u);
}
