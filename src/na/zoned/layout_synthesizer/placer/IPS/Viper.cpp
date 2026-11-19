
#include "na/zoned/layout_synthesizer/placer/IPS/Viper.hpp"

#include "na/zoned/Architecture.hpp"
#include "na/zoned/layout_synthesizer/placer/IPS/Common.hpp"

#include <array>
#include <format>
#include <iostream>

namespace na::zoned::ips::viper {

namespace {

struct Phase;

/**
 * @brief Represents a single or paired horizontal site entry within a column
 * pair.
 *
 * Encapsulates 1 site (for an odd/unpaired center column) or 2 sites in the
 * same row across paired columns, allowing two interacting qubits to occupy
 * adjacent columns.
 */
struct RowEntry {
  std::array<Site, 2> sites;
  std::uint8_t count; // 1 for single column, 2 for pair

  RowEntry(Site s1, Site s2) : sites{std::move(s1), std::move(s2)}, count{2} {}

  // Constructor for a single column (1 site)
  RowEntry(const Site& s1) : sites{s1, s1}, count{1} {}
};

/**
 * @brief Represents a row slot that has been occupied by physical qubits.
 */
struct TakenRow {
  RowEntry traps;
  std::array<qc::Qubit, 2> qubits;
};

/**
 * @brief Represents an interaction block of two-qubit gate pairs from a circuit
 * phase.
 *
 * Qubits interacting in the same phase are grouped into blocks to be placed
 * compactly in adjacent columns, minimizing shuttle routing overhead.
 */
struct Block {
  Block* parentBlock{nullptr};
  std::size_t id{0};
  std::size_t phaseIndex{0};
  double weight{0};

  // before being assigned
  std::vector<std::pair<qc::Qubit, qc::Qubit>> qubitPairs;
  // after being assigned
  std::vector<TakenRow> takenRows;

  [[nodiscard]] auto isSubblock() const noexcept -> bool {
    return parentBlock != nullptr;
  }
  [[nodiscard]] auto isAssigned() const noexcept -> bool {
    return !takenRows.empty();
  }

  void
  updateWeight(const std::unordered_map<qc::Qubit, double>& qubitWeighMap) {
    weight = 0;

    for (const auto& [q0, q1] : qubitPairs) {
      weight += (qubitWeighMap.at(q0) + qubitWeighMap.at(q1)) / 2.0;
    }
  }

  auto operator<(const Block& other) const noexcept -> bool {

    if (weight != other.weight) {
      return weight > other.weight;
    }
    if (phaseIndex != other.phaseIndex) {
      return phaseIndex < other.phaseIndex;
    } else {
      return id < other.id;
    }
  }
};

/**
 * @brief Represents a paired physical column (or single center column) in the
 * storage zone.
 *
 * Tracks available free row slots, occupied taken rows, and an aggregated
 * heuristic distance score representing proximity to the entanglement zones.
 */
struct ColumnPair {

  double score{0};
  std::size_t id{0};
  std::vector<RowEntry> freeRows;
  std::vector<TakenRow> takenRows;
  std::vector<RowEntry> prunedRows;

  [[nodiscard]] auto hasCapacity() const noexcept -> bool {
    return !freeRows.empty();
  }

  void updateScore(const std::unordered_map<Site, double>& siteScoreMap) {
    if (freeRows.empty()) {
      score = static_cast<double>(INTMAX_MAX);
      return;
    }

    score = 0;
    for (const auto& rowEntry : freeRows) {
      for (std::size_t i = 0; i < rowEntry.count; i++) {
        score += siteScoreMap.at(rowEntry.sites[i]);
      }
    }
    score /= static_cast<double>(freeRows.size()) * freeRows.front().count;
  }

  auto operator<(const ColumnPair& other) const noexcept -> bool {
    if (score != other.score) {
      return score < other.score;
    }
    return id < other.id;
  }
};

/**
 * @brief Encapsulates the overall state of the storage zone during VIPER
 * placement.
 *
 * Holds available column pairs, prioritized interaction blocks, height
 * thresholds, and site/qubit heuristic score mappings.
 */
struct PlacementZone {
  std::size_t nFreeColumnPairs;
  std::size_t heightLimit;
  std::unordered_map<std::size_t, std::size_t> nextId;

  std::unordered_map<Site, double> siteScoreMap;
  std::unordered_map<qc::Qubit, double> qubitWeightMap;

  std::set<ColumnPair> columnPairs;
  std::set<Block> blocks;
};

/**
 * @brief Initializes the placement zone from architecture specs and storage
 * layout.
 *
 * Groups storage columns into column pairs (or single center column if odd),
 * populates the site score map, and computes the vertical height limit.
 *
 * @param nQubits Total number of qubits to place.
 * @param architecture Target architecture description.
 * @param reverseInitialPlacement Whether to flip row orientation relative to
 * entanglement zones.
 * @return Initialized PlacementZone ready for block assignment.
 */
[[nodiscard]] auto initPlacementZone(const size_t nQubits,
                                     const Architecture& architecture,
                                     const bool reverseInitialPlacement)
    -> PlacementZone {

  const auto zonesSize = common::getStorageZoneSizes(architecture);
  const auto zoneSlots =
      common::getStorageZoneSlots(architecture, reverseInitialPlacement);
  const auto nCols = zonesSize.second;

  PlacementZone placementZone;
  std::vector<ColumnPair> temp;
  placementZone.siteScoreMap.insert(zoneSlots.begin(), zoneSlots.end());
  temp.resize((nCols + 1) / 2);

  if (nCols % 2 == 0) {
    for (std::size_t i = 0; i < zoneSlots.size(); i += 2) {

      auto& columnPair = temp[i % nCols / 2];
      RowEntry row_entry{zoneSlots[i].first, zoneSlots[i + 1].first};

      columnPair.freeRows.emplace_back(row_entry);
    }
  } else {
    const std::size_t midCol = nCols / 2;
    const std::size_t singleCol = midCol % 2 == 0 ? midCol : midCol + 1;

    for (std::size_t i = 0; i < zoneSlots.size();) {
      const std::size_t col = i % nCols;
      const std::size_t colIdx =
          col <= singleCol ? col / 2
                           : (singleCol / 2 + 1 + (col - singleCol - 1) / 2);
      auto& columnPair = temp[colIdx];

      if (col == singleCol) {
        RowEntry row_entry{zoneSlots[i].first};
        columnPair.freeRows.emplace_back(std::move(row_entry));
        columnPair.score += zoneSlots[i].second;
        i += 1;
      } else {
        RowEntry row_entry{zoneSlots[i].first, zoneSlots[i + 1].first};
        columnPair.freeRows.emplace_back(std::move(row_entry));
        i += 2;
      }
    }
  }

  for (std::size_t i = 0; i < temp.size(); ++i) {
    temp[i].id = i;
    temp[i].updateScore(placementZone.siteScoreMap);
    placementZone.columnPairs.insert(std::move(temp[i]));
  }

  placementZone.nFreeColumnPairs = temp.size();
  placementZone.heightLimit = common::getHeightLimit(nQubits, architecture);

  return placementZone;
}

/**
 * @brief Partitions two-qubit gate interactions into prioritized blocks.
 *
 * Scans gate layers within the lookahead window, groups interacting qubit pairs
 * into blocks, weights each block, and inserts them into the placement zone.
 *
 * @param placementZone Target placement zone to populate with blocks.
 * @param twoQubitGateLayers Circuit gate layers containing two-qubit
 * operations.
 * @param lookaheadWindow Number of initial gate layers to consider.
 */
[[nodiscard]] auto getQubitBlockAssignment(
    PlacementZone& placementZone,
    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers,
    std::size_t lookaheadWindow) {

  std::map<qc::Qubit, double> qubitMap =
      common::getQubitWeightMap(twoQubitGateLayers);

  placementZone.qubitWeightMap = {qubitMap.begin(), qubitMap.end()};
  std::size_t idx = 0;
  for (const auto& layer : twoQubitGateLayers) {
    Block block;
    std::vector<std::pair<std::pair<qc::Qubit, qc::Qubit>, double>>
        qubitPairList;
    block.qubitPairs.reserve(layer.size());
    qubitPairList.reserve(layer.size());

    block.phaseIndex = idx;
    placementZone.nextId[idx] = 0;

    for (const auto& pair : layer) {
      if (qubitMap.contains(pair[0]) && qubitMap.contains(pair[1])) {
        auto weight = (qubitMap[pair[0]] + qubitMap[pair[1]]) / 2.0;
        qubitPairList.emplace_back(pair, weight);
        qubitMap.erase(pair[0]);
        qubitMap.erase(pair[1]);
      }
    }

    std::ranges::sort(qubitPairList, [](const auto& a, const auto& b) {
      return a.second > b.second;
    });

    for (auto& pair : qubitPairList | std::views::keys) {
      block.qubitPairs.push_back(std::move(pair));
    }

    block.updateWeight(placementZone.qubitWeightMap);
    placementZone.blocks.insert(std::move(block));

    if (idx == lookaheadWindow) {
      break;
    }
    idx++;
  }
}

/**
 * @brief Assigns qubit pairs from a block into consecutive free rows of a
 * column pair.
 *
 * Extracts the block and column pair from the placement sets, maps qubit pairs
 * to physical row traps, updates scores and weights, and reinserts them into
 * the active sets.
 *
 * @param placementZone The active placement zone.
 * @param block The block whose qubit pairs are being placed.
 * @param columnPair The target column pair receiving the qubits.
 */
[[nodiscard]] auto colToBlockInsertHelper(PlacementZone& placementZone,
                                          const Block& block,
                                          const ColumnPair& columnPair) {

  auto columnPairNode = placementZone.columnPairs.extract(columnPair);
  auto blockNode = placementZone.blocks.extract(block);

  if (columnPairNode.empty() || blockNode.empty()) {
    placementZone.columnPairs.insert(std::move(columnPairNode));
    return;
  }

  auto& columnPairValue = columnPairNode.value();
  auto& blockValue = blockNode.value();

  for (const auto& pair : blockValue.qubitPairs) {
    RowEntry openRow = columnPairValue.freeRows.front();
    columnPairValue.freeRows.erase(columnPairValue.freeRows.begin());
    const auto [q0, q1] = pair;
    columnPairValue.takenRows.emplace_back(std::move(openRow),
                                           std::array{q0, q1});
    blockValue.takenRows.emplace_back(columnPairValue.takenRows.back());
  }

  blockValue.qubitPairs.clear();
  blockValue.qubitPairs.shrink_to_fit();
  columnPairValue.freeRows.shrink_to_fit();

  columnPairValue.updateScore(placementZone.siteScoreMap);
  blockValue.updateWeight(placementZone.qubitWeightMap);
  placementZone.columnPairs.insert(std::move(columnPairNode));
  placementZone.blocks.insert(std::move(blockNode));
}

/**
 * @brief Splits an oversized block into smaller subblocks to fit column height
 * capacity.
 *
 * Distributes qubit pairs across multiple child blocks so each fits within
 * available vertical limits, and reinserts them into the block queue.
 *
 * @note It does not assure that these blocks will fit into their respective
 * column pair when necessary, it might be needed to split them again in the
 * future.
 *
 * @param placementZone The active placement zone.
 * @param block The oversized block to split.
 * @param bestColumnPair Iterator pointing to the candidate column pair.
 */
[[nodiscard]] auto
colToBlockSplitHelper(PlacementZone& placementZone, const Block& block,
                      const std::set<ColumnPair>::iterator& bestColumnPair) {

  auto blockValue = placementZone.blocks.extract(block).value();

  auto& pairs = blockValue.qubitPairs;
  std::size_t coef = 0;
  if (placementZone.heightLimit == bestColumnPair->takenRows.size()) {
    coef = 2;
  } else {
    coef = (pairs.size() - 1) /
               (placementZone.heightLimit - bestColumnPair->takenRows.size()) +
           1;
  }
  std::vector<Block> newBlocks;
  newBlocks.resize(coef);

  for (std::size_t i = 0; i < pairs.size(); i++) {
    newBlocks[i % coef].qubitPairs.emplace_back(std::move(pairs[i]));
  }

  std::size_t idx = placementZone.nextId.at(block.phaseIndex);
  for (auto& newBlock : newBlocks) {
    newBlock.phaseIndex = block.phaseIndex;
    newBlock.id = idx++;
    newBlock.updateWeight(placementZone.qubitWeightMap);
    placementZone.blocks.insert(std::move(newBlock));
  }
  placementZone.nextId[block.phaseIndex] = idx;
}

/**
 * @brief Iteratively maps interaction blocks to column pairs.
 *
 * Places blocks into candidate column pairs respecting the height limit. If a
 * block exceeds available capacity, it is either placed within an allowed
 * overflow tolerance or split into subblocks.
 *
 * @param placementZone The placement zone state.
 * @param heightCoef Overflow threshold percentage for relaxing height limit
 * constraints.
 */
[[nodiscard]] auto assignBlockToCol(PlacementZone& placementZone,
                                    std::size_t heightCoef) {

  while (!placementZone.blocks.empty()) {
    if (placementZone.blocks.begin()->weight <= 0.0) {
      break;
    }

    auto& block = *placementZone.blocks.begin();

    if (const auto& bestColumnPair = placementZone.columnPairs.begin();
        block.qubitPairs.size() > bestColumnPair->freeRows.size()) {
      // printf("Split1\n");
      colToBlockSplitHelper(placementZone, block, bestColumnPair);
      continue;
    }

    bool abort = false;
    std::vector<std::pair<const ColumnPair*, std::size_t>> heightDelta;
    heightDelta.reserve(placementZone.columnPairs.size());
    for (const auto& columnPair : placementZone.columnPairs) {
      if (columnPair.freeRows.begin()->count == 1) {
        continue;
      }
      if (columnPair.takenRows.size() + block.qubitPairs.size() <=
          placementZone.heightLimit) {
        // printf("Insert1\n");
        colToBlockInsertHelper(placementZone, block, columnPair);
        abort = true;
        break;
      } else {
        heightDelta.emplace_back(&columnPair, block.qubitPairs.size() +
                                                  columnPair.takenRows.size() -
                                                  placementZone.heightLimit);
      }
    }

    if (abort) {
      continue;
    }

    // TODO: if above doesnt work i.e. no column can completely fill it up hunt
    // for best solution (split)
    for (const auto& [columnPair, delta] : heightDelta) {
      if (delta < 1 + static_cast<std::size_t>(
                          (static_cast<double>(heightCoef) / 100.0) *
                          static_cast<double>(block.qubitPairs.size()))) {
        if (block.qubitPairs.size() > columnPair->freeRows.size()) {
          continue;
        }
        // printf("Insert2\n");
        colToBlockInsertHelper(placementZone, block, *columnPair);
        abort = true;
        break;
      }
    }

    if (abort) {
      continue;
    }
    printf("Split2\n");
    colToBlockSplitHelper(placementZone, block,
                          placementZone.columnPairs.begin());
  }
}

/**
 * @brief Debug utility to print the occupancy, scores, and free rows of each
 * column pair.
 *
 * @param placementZone The placement zone to inspect.
 */
void printPlacementZoneColumns(const PlacementZone& placementZone) {
  std::cout << "\n==================== COLUMN PAIRS ("
            << placementZone.columnPairs.size()
            << " total) ====================\n";

  for (const auto& cp : placementZone.columnPairs) {
    std::cout << std::format(
        "[ColumnPair ID: {}] Score: {:.4f} | Free: {} rows | Taken: {} rows\n",
        cp.id, cp.score, cp.freeRows.size(), cp.takenRows.size());

    // Print taken rows and which qubits occupy them
    if (!cp.takenRows.empty()) {
      std::cout << "  Taken Rows:\n";
      for (std::size_t i = 0; i < cp.takenRows.size(); ++i) {
        const auto& tr = cp.takenRows[i];
        std::cout << std::format("    [{:2d}] Qubits: [{}, {}] -> Sites: ", i,
                                 tr.qubits[0], tr.qubits[1]);
        for (std::size_t s = 0; s < tr.traps.count; ++s) {
          const auto& [slm, r, c] = tr.traps.sites[s];
          std::cout << std::format("(SLM {}, r:{}, c:{}) ", slm.get().id, r, c);
        }
        std::cout << "\n";
      }
    }

    // Print remaining free rows
    if (!cp.freeRows.empty()) {
      std::cout << std::format("  Free Rows (showing up to 3 of {}):\n",
                               cp.freeRows.size());
      for (std::size_t i = 0; i < std::min<std::size_t>(3, cp.freeRows.size());
           ++i) {
        const auto& fr = cp.freeRows[i];
        std::cout << std::format("    [{:2d}] Sites: ", i);
        for (std::size_t s = 0; s < fr.count; ++s) {
          const auto& [slm, r, c] = fr.sites[s];
          std::cout << std::format("(SLM {}, r:{}, c:{}) ", slm.get().id, r, c);
        }
        std::cout << "\n";
      }
    } else {
      std::cout
          << "  ** DEPLETED (freeRows is empty, score is INTMAX_MAX) **\n";
    }
    std::cout << "--------------------------------------------------------\n";
  }
}

/**
 * @brief Finalizes the mapping of all qubits to physical sites.
 *
 * Maps all qubits placed via block-to-column assignment, then greedily assigns
 * any remaining unplaced or single-qubit nodes into open row slots according to
 * weight.
 *
 * @param nQubits Total number of qubits in the circuit.
 * @param placementZone Placement zone containing populated column pairs and
 * blocks.
 * @param heightCoef Height coefficient used for block-to-column assignment.
 * @return Complete mapping from logical qubits to physical sites.
 */
[[nodiscard]] auto getAssignmentsColBlock(size_t nQubits,
                                          PlacementZone& placementZone,
                                          std::size_t heightCoef)
    -> std::unordered_map<qc::Qubit, Site> {

  assignBlockToCol(placementZone, heightCoef);
  std::unordered_map<qc::Qubit, Site> result;
  result.reserve(nQubits);

  printf("Initialising Qubit to Site placement...");
  for (const auto& block : placementZone.blocks) {
    if (!block.qubitPairs.empty()) {
      throw std::runtime_error(
          "Not all qubits were assigned, initial placement error");
    }

    for (const auto& row : block.takenRows) {
      for (std::size_t i = 0; i < row.qubits.size(); i++) {
        result.insert_or_assign(row.qubits[i], row.traps.sites[i]);
      }
    }
  }

  std::vector<std::pair<qc::Qubit, double>> qubitWeightList = {
      placementZone.qubitWeightMap.begin(), placementZone.qubitWeightMap.end()};
  std::ranges::sort(qubitWeightList, [](const auto& a, const auto& b) {
    return a.second > b.second;
  });
  std::optional<Site> freeSite;
  qc::Qubit* pendingQubitSlot = nullptr;

  for (std::size_t i = 0; i < nQubits; i++) {
    if (!placementZone.qubitWeightMap.contains(static_cast<qc::Qubit>(i))) {
      qubitWeightList.emplace_back(i, 0);
    }
  }

  for (const auto& qubit : qubitWeightList | std::views::keys) {
    if (!result.contains(qubit)) {
      if (freeSite) {
        result.insert_or_assign(qubit, *freeSite);
        *pendingQubitSlot = qubit;
        freeSite.reset();
        pendingQubitSlot = nullptr;
        continue;
      }

      auto it = std::ranges::find_if(
          placementZone.columnPairs, [&](const ColumnPair& cp) {
            return cp.takenRows.size() < placementZone.heightLimit;
          });
      if (it == placementZone.columnPairs.end()) {
        it = placementZone.columnPairs.begin();
      }

      auto columnPairNode = placementZone.columnPairs.extract(it);
      auto& columnPairValue = columnPairNode.value();

      if (columnPairValue.freeRows.empty()) {
        printf("%s\n", std::format("Score: {}", columnPairValue.score).c_str());
        placementZone.columnPairs.insert(std::move(columnPairNode));

        std::unordered_map<qc::Qubit, Site> map;
        map.reserve(nQubits);
        for (const auto& [qubit2, site] : result) {
          printf(" %zu,", static_cast<std::size_t>(qubit2));
          if (map.contains(qubit2)) {
            const auto& [slm1, row1, col1] = site;
            const auto& [slm2, row2, col2] = map.at(qubit2);
            throw std::runtime_error(
                std::format("Qubit found twice, ({}, {}) and ({}, {})", row1,
                            col1, row2, col2));
          } else {
            map.insert_or_assign(qubit2, site);
          }
        }
        printf("nQubits: %zu", nQubits);
        throw std::runtime_error("Too many qubits for architecture");
      }

      RowEntry openRow = std::move(columnPairValue.freeRows.front());
      columnPairValue.freeRows.erase(columnPairValue.freeRows.begin());

      result.insert_or_assign(qubit, openRow.sites[0]);

      if (openRow.count == 2) {
        freeSite.emplace(openRow.sites[1]);
      }

      columnPairValue.takenRows.push_back(
          {std::move(openRow),
           std::array{qubit, static_cast<qc::Qubit>(INTMAX_MAX)}});

      if (openRow.count == 2) {
        pendingQubitSlot = &columnPairValue.takenRows.back().qubits[1];
      }

      columnPairValue.updateScore(placementZone.siteScoreMap);
      placementZone.columnPairs.insert(std::move(columnPairNode));
    }
  }
  return result;
}

} // namespace

/**
 * @brief Computes the initial placement of qubits onto storage sites using
 * VIPER.
 *
 * Coordinates placement zone initialization, block extraction from two-qubit
 * gate layers, column-block mapping, and final site assignment.
 *
 * @param nQubits Total number of logical qubits to place.
 * @param twoQubitGateLayers Layers of 2-qubit gates representing circuit
 * interactions.
 * @param architecture Target zoned neutral-atom hardware architecture.
 * @param reverseInitialPlacement If true, reverses initial row orientation
 * relative to entanglement zones.
 * @return Initial mapping from logical qubits to physical sites (SLM traps).
 */
auto getAssignments(const std::size_t nQubits,
                    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers,
                    const Architecture& architecture,
                    const bool reverseInitialPlacement)
    -> std::unordered_map<qc::Qubit, Site> {
  printf("Initialising Placement Zone...\n");
  PlacementZone placement_zone =
      initPlacementZone(nQubits, architecture, reverseInitialPlacement);
  printf("Assigning Qubits to Blocks...\n");
  getQubitBlockAssignment(placement_zone, twoQubitGateLayers,
                          twoQubitGateLayers.size());
  printf("Initialising Assigning Blocks to Columns...\n");
  return getAssignmentsColBlock(nQubits, placement_zone, 20);
}

} // namespace na::zoned::ips::viper
