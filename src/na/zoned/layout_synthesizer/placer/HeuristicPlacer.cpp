/*
 * Copyright (c) 2023 - 2026 Chair for Design Automation, TUM
 * Copyright (c) 2025 - 2026 Munich Quantum Software Company GmbH
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

#include "na/zoned/layout_synthesizer/placer/HeuristicPlacer.hpp"

#include "ir/Definitions.hpp"
#include "na/zoned/Architecture.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <spdlog/spdlog.h>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na::zoned {

auto HeuristicPlacer::isGoal(const size_t nGates, const GateNode& node)
    -> bool {
  return node.level == nGates;
}
auto HeuristicPlacer::isGoal(const size_t nAtoms, const AtomNode& node)
    -> bool {
  return node.level == nAtoms;
}
auto HeuristicPlacer::discretizePlacementOfAtoms(
    const Placement& placement, const std::vector<qc::Qubit>& atoms) const
    -> std::pair<RowColumnMap<uint8_t>, RowColumnMap<uint8_t>> {
  std::map<size_t, RowColumnSet> rows;
  std::map<size_t, RowColumnSet> columns;
  for (const auto atom : atoms) {
    const auto& [slm, r, c] = placement[atom];
    const auto& [x, y] = architecture_.get().exactSLMLocation(slm, r, c);
    rows.try_emplace(y).first->second.emplace(std::cref(slm), r);
    columns.try_emplace(x).first->second.emplace(std::cref(slm), c);
  }
  RowColumnMap<uint8_t> rowIndices;
  uint8_t rowIndex = 0;
  for (const auto& sites : rows | std::views::values) {
    for (const auto& site : sites) {
      rowIndices.emplace(site, rowIndex);
    }
    ++rowIndex;
  }
  RowColumnMap<uint8_t> columnIndices;
  uint8_t columnIndex = 0;
  for (const auto& sites : columns | std::views::values) {
    for (const auto& site : sites) {
      columnIndices.emplace(site, columnIndex);
    }
    ++columnIndex;
  }
  return std::pair{rowIndices, columnIndices};
}

auto HeuristicPlacer::discretizeNonOccupiedStorageSites(
    const SiteSet& occupiedSites) const
    -> std::pair<RowColumnMap<uint8_t>, RowColumnMap<uint8_t>> {
  std::map<size_t, std::pair<std::reference_wrapper<const SLM>, size_t>> rows;
  std::map<size_t, std::pair<std::reference_wrapper<const SLM>, size_t>>
      columns;
  for (const auto& slm : architecture_.get().storageZones) {
    // find rows with free sites
    for (size_t r = 0; r < slm->nRows; ++r) {
      for (size_t c = 0; c < slm->nCols; ++c) {
        if (!occupiedSites.contains(std::tie(*slm, r, c))) {
          // free site in row r found at column c
          rows.emplace(slm->location.second + (slm->siteSeparation.second * r),
                       std::pair{std::cref(*slm), r});
          break;
        }
      }
    }
    // find columns with free sites
    for (size_t c = 0; c < slm->nCols; ++c) {
      for (size_t r = 0; r < slm->nRows; ++r) {
        if (!occupiedSites.contains(std::tie(*slm, r, c))) {
          // free site in column c found at row r
          columns.emplace(slm->location.first + (slm->siteSeparation.first * c),
                          std::pair{std::cref(*slm), c});
          break;
        }
      }
    }
  }
  RowColumnMap<uint8_t> rowIndices;
  uint8_t rowIndex = 0;
  for (const auto& site : rows | std::views::values) {
    rowIndices.emplace(site, rowIndex++);
  }
  RowColumnMap<uint8_t> columnIndices;
  uint8_t columnIndex = 0;
  for (const auto& site : columns | std::views::values) {
    columnIndices.emplace(site, columnIndex++);
  }
  return std::pair{rowIndices, columnIndices};
}

auto HeuristicPlacer::discretizeNonOccupiedEntanglementSites(
    const SiteSet& occupiedSites) const
    -> std::pair<RowColumnMap<uint8_t>, RowColumnMap<uint8_t>> {
  std::map<size_t, RowColumnSet> rows;
  std::map<size_t, RowColumnSet> columns;
  for (const auto& zone : architecture_.get().entanglementZones) {
    for (const auto& slm : *zone) {
      // find rows with free sites
      for (size_t r = 0; r < slm.nRows; ++r) {
        for (size_t c = 0; c < slm.nCols; ++c) {
          if (!occupiedSites.contains(std::tie(slm, r, c))) {
            // free site in row r found at column c
            rows.try_emplace(slm.location.second +
                             (slm.siteSeparation.second * r))
                .first->second.emplace(std::cref(slm), r);
            break;
          }
        }
      }
      // find columns with free sites
      for (size_t c = 0; c < slm.nCols; ++c) {
        for (size_t r = 0; r < slm.nRows; ++r) {
          if (!occupiedSites.contains(std::tie(slm, r, c))) {
            // free site in column c found at row r
            columns
                .try_emplace(slm.location.first +
                             (slm.siteSeparation.first * c))
                .first->second.emplace(std::cref(slm), c);
            break;
          }
        }
      }
    }
  }
  RowColumnMap<uint8_t> rowIndices;
  uint8_t rowIndex = 0;
  for (const auto& sites : rows | std::views::values) {
    for (const auto& site : sites) {
      rowIndices.emplace(site, rowIndex);
    }
    ++rowIndex;
  }
  RowColumnMap<uint8_t> columnIndices;
  uint8_t columnIndex = 0;
  for (const auto& sites : columns | std::views::values) {
    for (const auto& site : sites) {
      columnIndices.emplace(site, columnIndex);
    }
    ++columnIndex;
  }
  return std::pair{rowIndices, columnIndices};
}

auto HeuristicPlacer::getStorageZoneSizes() const
    -> std::pair<std::size_t, std::size_t> {

  auto& storageZones = architecture_.get().storageZones;
  std::size_t nRowsTotal = 0;
  std::size_t nColsTotal = 0;
  for (const auto& zonePointer : storageZones) {
    const auto& storageZone = *zonePointer;

    nRowsTotal += storageZone.nRows;
    nColsTotal += storageZone.nCols;
  }

  return std::pair {nRowsTotal, nColsTotal};
}

//TODO: Doesn't yet work for multiple storage zones (change rows and cols assigned to Site)
auto HeuristicPlacer::getStorageZoneSlots() const
    -> std::vector<std::pair<Site, double>> {

  std::vector<std::pair<Site, double>> zoneSlots;
  auto& storageZones = architecture_.get().storageZones;
  auto& entanglementZones = architecture_.get().entanglementZones;

  for (const auto& slmPointer : storageZones) {
    const auto& slm = *slmPointer;

    auto processRow = [&](std::size_t row) {
      for (std::size_t col = 0; col < slm.nCols; ++col) {

        double averageProximity = 0;
        std::size_t count = 0;

        for (const auto& slmPair : entanglementZones) {
          const auto& slm1 = (*slmPair)[0];
          const auto& slm2 = (*slmPair)[1];

          for (std::size_t rowEnt = 0; rowEnt < slm1.nRows; ++rowEnt) {
            for (std::size_t colEnt = 0; colEnt < slm1.nCols; ++colEnt) {

              const double distance1 = architecture_.get().distance(
                  slm, row, col, slm1, rowEnt, colEnt);
              const double distance2 = architecture_.get().distance(
                  slm, row, col, slm2, rowEnt, colEnt);

              averageProximity += distance1 + distance2;
            }
          }
          count += slm1.nRows * slm1.nCols + slm2.nRows * slm2.nCols;
        }

        if (count > 0) {
          averageProximity /= static_cast<double>(count);
        }

        Site site(std::cref(slm), row, col);
        zoneSlots.emplace_back(site, averageProximity);
      }
    };

    if (reverseInitialPlacement_) {
      for (std::size_t row = slm.nRows - 1; ; row--) {
        processRow(row);

        if (row == 0) {
          break;
        }
      }
    } else {
      for (std::size_t row = 0; row < slm.nRows; ++row) {
        processRow(row);
      }
    }
  }

  return zoneSlots;
}

auto HeuristicPlacer::getHeightLimit(std::size_t nQubits) const
    -> std::size_t {

  auto nColsTotal = getStorageZoneSizes().second;
  return 1 + (nQubits - 1) / nColsTotal;
}

auto HeuristicPlacer::initPlacementZone(const size_t nQubits) const
    -> PlacementZone {

  const auto zonesSize = getStorageZoneSizes();
  const auto zoneSlots = getStorageZoneSlots();
  const auto nCols = zonesSize.second;

  PlacementZone placementZone;
  std::vector<ColumnPair> temp;
  placementZone.siteScoreMap.insert(zoneSlots.begin(), zoneSlots.end());
  temp.resize((nCols + 1) / 2);


  if (nCols % 2 == 0) {
    for (std::size_t i = 0; i < zoneSlots.size(); i += 2) {

      auto& columnPair = temp[i % nCols / 2];
      RowEntry row_entry{zoneSlots[i].first, zoneSlots[i+1].first};

      columnPair.freeRows.emplace_back(row_entry);
    }
  } else {
    const std::size_t midCol = nCols / 2;
    const std::size_t singleCol = midCol % 2 == 0 ? midCol : midCol + 1;

    for (std::size_t i = 0; i < zoneSlots.size(); ) {
      const std::size_t col = i % nCols;
      const std::size_t colIdx = col <= singleCol ? col / 2 : (singleCol/2 + 1 + (col - singleCol - 1) / 2);
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
  placementZone.heightLimit = getHeightLimit(nQubits);

  return placementZone;
}

auto HeuristicPlacer::getQubitWeightMap(
    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers)
    -> std::map<qc::Qubit, double> {

  std::map<qc::Qubit, double> qubitMap;
  std::size_t idx = 0;

  for (const auto& layer : twoQubitGateLayers) {
    double weight = 1.0 /(static_cast<double>(idx + 1));

    for (const auto& pair : layer) {
      qubitMap[pair[0]] += weight;
      qubitMap[pair[1]] += weight;
    }
    idx++;
  }
  return qubitMap;
}

auto HeuristicPlacer::getQubitWeightList(
    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers)
    -> std::vector<std::pair<qc::Qubit, double>> {

  std::map<qc::Qubit, double> qubitMap = getQubitWeightMap(twoQubitGateLayers);
  return {qubitMap.begin(), qubitMap.end()};
}

auto HeuristicPlacer::getQubitPairWeightList(
    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers ,
    const std::size_t lookaheadWindow)
    -> std::vector<std::pair<std::pair<qc::Qubit, qc::Qubit>, double>> {

  std::map<qc::Qubit, double> qubitMap = getQubitWeightMap(twoQubitGateLayers);
  std::map<std::pair<qc::Qubit, qc::Qubit>, double> result;
  std::size_t idx = 0;

  idx = 0;
  for (const auto& layer : twoQubitGateLayers) {

    for (const auto& pair : layer) {

      result[pair] += (qubitMap[pair[0]] + qubitMap[pair[1]]) / 2;
    }

    idx++;

    if (idx == lookaheadWindow) {
      break;
    }
  }

  return {result.begin(), result.end()};
}

auto HeuristicPlacer::getQubitBlockAssignment(
    PlacementZone& placementZone,
    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers,
    const std::size_t lookaheadWindow) const {

  std::map<qc::Qubit, double> qubitMap = getQubitWeightMap(twoQubitGateLayers);

  placementZone.qubitWeightMap = {qubitMap.begin(), qubitMap.end()};
  std::size_t idx = 0;
  for (const auto& layer : twoQubitGateLayers) {
    Block block;
    std::vector<std::pair<std::pair<qc::Qubit, qc::Qubit>, double>> qubitPairList;
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

    std::ranges::sort(qubitPairList, [](const auto& a, const auto& b){return a.second > b.second;});

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

auto HeuristicPlacer::getFirstStages(const std::size_t nQubits,
    const std::size_t lookaheadWindow,
    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers) const
    -> std::unordered_map<qc::Qubit, std::size_t> {

  std::unordered_map<qc::Qubit, std::size_t> qubitMap;
  std::size_t idx = 0;

  for (const auto& layer : twoQubitGateLayers) {
    for (const auto& pair : layer) {

      if (!qubitMap.contains(pair[0])) {
        qubitMap[pair[0]] = idx;
      }
      if (!qubitMap.contains(pair[1])) {
        qubitMap[pair[1]] = idx;
      }
    }

    if (idx == lookaheadWindow || qubitMap.size() == nQubits) {
      break;
    }
    idx++;
  }
  return qubitMap;
}

auto HeuristicPlacer::getFirstStagesPair(const std::size_t lookaheadWindow,
    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers) const
    -> std::unordered_map<std::pair<qc::Qubit, qc::Qubit>, std::size_t> {

  std::unordered_map<std::pair<qc::Qubit, qc::Qubit>, std::size_t> qubitMap;
  std::size_t idx = 0;

  for (const auto& layer : twoQubitGateLayers) {
    for (const auto& pair : layer) {

      if (!qubitMap.contains(pair)) {
        qubitMap[pair] = idx;
      }
    }

    if (idx == lookaheadWindow) {
      break;
    }
    idx++;
  }
  return qubitMap;
}

auto HeuristicPlacer::hasConflict1D(std::size_t k_cand, std::size_t v_cand, std::size_t k, std::size_t v) -> bool {
  if (k_cand == k) {
    return v_cand != v;
  }
  if (k_cand < k) {
    return v_cand >= v;
  }
  // k_cand > k
  return v_cand <= v;
}

auto HeuristicPlacer::countConflictsWithGroup(
    const Movement& move_cand,
    const std::vector<Movement>& group) const
    -> size_t {
  size_t conflicts = 0;

  for (const auto& movement : group) {

    if (hasConflict1D(move_cand.srcRow, move_cand.dstRow, movement.srcRow, movement.dstRow)
        || hasConflict1D(move_cand.srcCol, move_cand.dstCol, movement.srcCol, movement.dstCol)) {
      conflicts++;
    }
  }

  return conflicts;
}

//TODO for initial smart split
/*
auto HeuristicPlacer::blockToColStats( const Block& block,
    const std::vector<ColumnPair*>& columnPairs,
    const std::size_t maxHeight) const
    -> std::pair<bool, std::size_t> {

  if (block.qubitPairs.size() > maxHeight * columnPairs.size()) {
    return {false, 0};
  }

  for (const auto& pair : block.qubitPairs) {

  }
}*/

auto HeuristicPlacer::getAssignments(const std::vector<std::pair<Site, double>>& zoneSlots,
    const std::vector<std::pair<qc::Qubit, double>>& qubitWeightList,
    const std::unordered_map<qc::Qubit, std::size_t>& firstStages,
    const std::size_t nQubits) const
    -> std::unordered_map<qc::Qubit, Site> {

  std::unordered_map zoneSlotsMap(zoneSlots.begin(), zoneSlots.end());
  std::vector<Movement> group;
  std::unordered_map<qc::Qubit, Site> results;

  for (const auto& qc : qubitWeightList | std::views::keys) {

    std::pair<std::optional<std::pair<Site, Site>>, double> minScore = {std::nullopt, 1e308};
    auto minScoreMove = Movement();

    for (const auto& site : zoneSlots | std::views::keys) {
      constexpr double lamda = 10.0;

      if (!zoneSlotsMap.contains(site)) {
        continue;
      }

      auto& slm = std::get<0>(site);
      auto row = std::get<1>(site);
      auto col = std::get<2>(site);

      const auto& entangleSite =
            architecture_.get().nearestEntanglementSite(slm, row, col, slm, row, col);
      const auto& [entangleSLM, entangleRow, entangleCol] = entangleSite;

      auto& [x0,y0] = slm.get().location;
      auto rowCoord = y0 + slm.get().siteSeparation.second * row;
      auto colCoord = x0 + slm.get().siteSeparation.first * col;

      auto& [entx0,enty0] = entangleSLM.get().location;
      auto entRowCoord = enty0 + entangleSLM.get().siteSeparation.second * entangleRow;
      auto entColCoord = entx0 + entangleSLM.get().siteSeparation.first * entangleCol;
      auto movement = Movement(rowCoord, entRowCoord, colCoord, entColCoord);

      std::size_t count = 0;
      if (firstStages.contains(qc)) {
        count = countConflictsWithGroup(movement , group);
      }
      double distance = architecture_.get()
        .nearestEntanglementSiteDistance(slm, row, col, slm, row, col);

      double score = lamda * static_cast<double>(count) + distance;

      if (score < minScore.second) {
        minScore.first = {site, entangleSite};
        minScore.second = score;
        minScoreMove = movement;
      }
    }

    if (!minScore.first.has_value()) {
      throw std::runtime_error("Error: No available storage slot found for qubit assignment!");
    }

    const auto& site = minScore.first->first;
    zoneSlotsMap.erase(site);
    if (firstStages.contains(qc)) {
      group.emplace_back(minScoreMove);
    }
    results.emplace(qc, site);
  }

  for (std::size_t i = 0; i < nQubits; ++i) {
    const auto qc = static_cast<qc::Qubit>(i);

    if (!results.contains(qc)) {

      std::pair<std::optional<Site>, double> minScore = {std::nullopt, 1e308};
      for (const auto& site : zoneSlots | std::views::keys) {
        if (!zoneSlotsMap.contains(site)) {
          continue;
        }

        auto& slm = std::get<0>(site);
        auto row = std::get<1>(site);
        auto col = std::get<2>(site);


        double distance = architecture_.get()
        .nearestEntanglementSiteDistance(slm, row, col, slm, row, col);

        if (distance < minScore.second) {
          minScore.first = site;
          minScore.second = distance;
        }
      }

      const auto& site = minScore.first;

      if (!site.has_value()) {
        throw std::runtime_error("Error: No available storage slot found for qubit assignment!");
      }

      zoneSlotsMap.erase(*site);
      results.emplace(qc, *site);
    }
  }

  return results;
}

auto HeuristicPlacer::getAssignmentsBasicPair( //please only use with lookaheadWindow = 1
    const std::vector<std::pair<Site, double>>& zoneSlots,
    const std::vector<std::pair<std::pair<qc::Qubit, qc::Qubit>, double>>& qubitWeightPairList,
    const std::vector<std::pair<qc::Qubit, double>>& qubitWeightList) const
    -> std::unordered_map<qc::Qubit, Site>{

  std::unordered_map<qc::Qubit, Site> result;
  std::size_t idx = 0;

  for (const auto& pair : qubitWeightPairList | std::views::keys) {
    result.emplace(pair.first, zoneSlots[idx++].first);
    result.emplace(pair.second, zoneSlots[idx++].first);
  }

  for (const auto& qubit : qubitWeightList | std::views::keys) {

    if (!result.contains(qubit)) {
      result.emplace(qubit, zoneSlots[idx++].first);
    }
  }

  return result;
}

auto HeuristicPlacer::colToBlockSplitHelper(
    PlacementZone& placementZone,
    const Block& block,
    const std::set<ColumnPair>::iterator& bestColumnPair) const {

  auto blockValue = placementZone.blocks.extract(block).value();

  auto& pairs = blockValue.qubitPairs;
  std::size_t coef = 0;
  if (placementZone.heightLimit == bestColumnPair->takenRows.size()) {
    coef = 2;
  } else {
    coef = coef = (pairs.size() - 1) / (placementZone.heightLimit - bestColumnPair->takenRows.size()) + 1;
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

auto HeuristicPlacer::colToBlockInsertHelper(
    PlacementZone& placementZone,
    const Block& block,
    const ColumnPair& columnPair) const {

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
    const auto [q0,q1] = pair;
    columnPairValue.takenRows.emplace_back(std::move(openRow), std::array{q0, q1});
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

auto HeuristicPlacer::assignBlockToCol(
    PlacementZone& placementZone,
    const std::size_t heightCoef) const {

  while (!placementZone.blocks.empty())  {
    if (placementZone.blocks.begin()->weight <= 0.0) {
      break;
    }

    auto& block = *placementZone.blocks.begin();

    if (const auto& bestColumnPair = placementZone.columnPairs.begin();
        block.qubitPairs.size() > placementZone.heightLimit) {
      printf("Split1\n");
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
      if (columnPair.takenRows.size() + block.qubitPairs.size() <= placementZone.heightLimit) {
        printf("Insert1\n");
        colToBlockInsertHelper(placementZone, block, columnPair);
        abort = true;
        break;
      } else {
        heightDelta.emplace_back(&columnPair,block.qubitPairs.size() + columnPair.takenRows.size() - placementZone.heightLimit);
      }
    }

    if (abort) {
      continue;
    }

    //TODO: if above doesnt work i.e. no column can completely fill it up hunt for best solution (split)
    for (const auto& [columnPair, delta] : heightDelta) {
      if (bool condition =
              delta < static_cast<std::size_t>(
                          (static_cast<double>(heightCoef) / 100.0) *
                          static_cast<double>(block.qubitPairs.size()))) {
        if (block.qubitPairs.size() > columnPair->freeRows.size()) {
          continue;
        }
        printf("Insert2\n");
        colToBlockInsertHelper(placementZone, block, *columnPair);
        abort = true;
        break;
      }
    }

    if (abort) {
      continue;
    }
    printf("Split2\n");
    colToBlockSplitHelper(placementZone, block, placementZone.columnPairs.begin());
  }
}

void printPlacementZoneColumns(const PlacementZone& placementZone) {
  std::cout << "\n==================== COLUMN PAIRS ("
            << placementZone.columnPairs.size() << " total) ====================\n";

  for (const auto& cp : placementZone.columnPairs) {
    std::cout << std::format("[ColumnPair ID: {}] Score: {:.4f} | Free: {} rows | Taken: {} rows\n",
                             cp.id, cp.score, cp.freeRows.size(), cp.takenRows.size());

    // Print taken rows and which qubits occupy them
    if (!cp.takenRows.empty()) {
      std::cout << "  Taken Rows:\n";
      for (std::size_t i = 0; i < cp.takenRows.size(); ++i) {
        const auto& tr = cp.takenRows[i];
        std::cout << std::format("    [{:2d}] Qubits: [{}, {}] -> Sites: ",
                                 i, tr.qubits[0], tr.qubits[1]);
        for (std::size_t s = 0; s < tr.traps.count; ++s) {
          const auto& [slm, r, c] = tr.traps.sites[s];
          std::cout << std::format("(SLM {}, r:{}, c:{}) ", slm.get().id, r, c);
        }
        std::cout << "\n";
      }
    }

    // Print remaining free rows
    if (!cp.freeRows.empty()) {
      std::cout << std::format("  Free Rows (showing up to 3 of {}):\n", cp.freeRows.size());
      for (std::size_t i = 0; i < std::min<std::size_t>(3, cp.freeRows.size()); ++i) {
        const auto& fr = cp.freeRows[i];
        std::cout << std::format("    [{:2d}] Sites: ", i);
        for (std::size_t s = 0; s < fr.count; ++s) {
          const auto& [slm, r, c] = fr.sites[s];
          std::cout << std::format("(SLM {}, r:{}, c:{}) ", slm.get().id, r, c);
        }
        std::cout << "\n";
      }
    } else {
      std::cout << "  ** DEPLETED (freeRows is empty, score is INTMAX_MAX) **\n";
    }
    std::cout << "--------------------------------------------------------\n";
  }
}

auto HeuristicPlacer::getAssignmentsColBlock(const size_t nQubits,
    PlacementZone& placementZone,
    const std::size_t heightCoef) const
    -> std::unordered_map<qc::Qubit, Site> {

  assignBlockToCol(placementZone, heightCoef);
  std::unordered_map<qc::Qubit, Site> result;
  result.reserve(nQubits);

  printf("Initialising Qubit to Site placement...");
  for (const auto& block : placementZone.blocks) {
    if (!block.qubitPairs.empty()) {
      throw std::runtime_error("Not all qubits were assigned, initial placement error");
    }

    for (const auto& row : block.takenRows) {
      for (std::size_t i = 0; i < row.qubits.size(); i++) {
        result.insert_or_assign(row.qubits[i], row.traps.sites[i]);
      }
    }
  }


  std::vector<std::pair<qc::Qubit, double>> qubitWeightList = {placementZone.qubitWeightMap.begin(), placementZone.qubitWeightMap.end()};
  std::ranges::sort(qubitWeightList, [](const auto& a, const auto& b){return a.second > b.second;});
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

      auto it = std::ranges::find_if(placementZone.columnPairs, [&](const ColumnPair& cp) {
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
        printPlacementZoneColumns(placementZone);

        std::unordered_map<qc::Qubit, Site> map;
        map.reserve(nQubits);
        for (const auto& [qubit, site] : result) {
          printf(" %zu,", static_cast<std::size_t>(qubit));
          if (map.contains(qubit)) {
            const auto& [slm1, row1, col1] = site;
            const auto& [slm2, row2, col2] = map.at(qubit);
            throw std::runtime_error(std::format("Qubit found twice, ({}, {}) and ({}, {})", row1, col1, row2, col2));
          } else {
            map.insert_or_assign(qubit, site);
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

      columnPairValue.takenRows.push_back({
          std::move(openRow),
          std::array{qubit, static_cast<qc::Qubit>(INTMAX_MAX)}
      });

      if (openRow.count == 2) {
        pendingQubitSlot = &columnPairValue.takenRows.back().qubits[1];
      }

      columnPairValue.updateScore(placementZone.siteScoreMap);
      placementZone.columnPairs.insert(std::move(columnPairNode));
    }
  }
  return result;
}

auto HeuristicPlacer::makeInitialPlacement(const size_t nQubits) const
    -> Placement {
  auto slmIt = architecture_.get().storageZones.cbegin();
  std::size_t c = 0;
  std::int64_t r = reverseInitialPlacement_
                       ? static_cast<std::int64_t>((*slmIt)->nRows) - 1
                       : 0;
  const std::int64_t step = reverseInitialPlacement_ ? -1 : 1;
  Placement initialPlacement;
  initialPlacement.reserve(nQubits);
  for (qc::Qubit qubit = 0; qubit < nQubits; ++qubit) {
    initialPlacement.emplace_back(**slmIt, r, c++);
    if (c == (*slmIt)->nCols) {
      // the end of the row reached, go to the next row
      r += step;
      c = 0;
      if (step == 1 ? r == static_cast<std::int64_t>((*slmIt)->nRows)
                    : r == -1) {
        // the end of the slm reached, go to the next slm
        ++slmIt;
        r = step == 1 ? static_cast<std::int64_t>((*slmIt)->nRows) - 1 : 0;
                    }
    }
  }
  return initialPlacement;
}

auto HeuristicPlacer::testOutput(
    const std::unordered_map<qc::Qubit, Site>& assignments, std::size_t nQubits) const {

  assert(assignments.size() == nQubits && "Not all qubits were assigned");

  std::unordered_map<Site, std::size_t> tmp;
  tmp.reserve(nQubits);
  for (std::size_t i = 0; i < nQubits; i++) {
    const auto& site = assignments.at(static_cast<qc::Qubit>(i));
    const auto [it, inserted] = tmp.try_emplace(site, i);
    assert(inserted || !std::format("Mapping qubit {} to same site as qubit {}", i, it->second).empty());
  }
}

auto HeuristicPlacer::makeInitialPlacementMulti(const size_t nQubits,
    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers,
    Config::InitialPlacementStrategy ips) const
    -> Placement {

  Placement initialPlacement;
  initialPlacement.reserve(nQubits);

  std::vector<std::pair<Site, double>> zoneSlots = getStorageZoneSlots();

  std::printf("Slots: %zu, Qubits: %zu\n", zoneSlots.size(), nQubits);
  if (zoneSlots.size() < nQubits) {
    throw std::invalid_argument("Too many qubits for this architecture.");
  }

  std::unordered_map<qc::Qubit, Site> assignments;
  switch (ips) {

    case InitialPlacementStrategy::ZAP: {
      std::ranges::sort(zoneSlots, [](const auto& a, const auto& b){return a.second < b.second;});
      if (zoneSlots.size() > nQubits) {
        zoneSlots.erase(zoneSlots.begin() + static_cast<std::ptrdiff_t>(nQubits), zoneSlots.end());
      }
      std::vector<std::pair<qc::Qubit, double>> qubitWeightList = getQubitWeightList(twoQubitGateLayers);
      constexpr std::size_t lookaheadWindow = 10;  //Can be changed but defined as 2 in ZAP
      std::unordered_map<qc::Qubit, std::size_t> firstStages = getFirstStages(nQubits, lookaheadWindow, twoQubitGateLayers);
      std::ranges::sort(qubitWeightList, [](const auto& a, const auto& b){return a.second > b.second;});

      assignments = getAssignments(zoneSlots, qubitWeightList, firstStages, nQubits);
      break;
    }

    case InitialPlacementStrategy::BasicPair: {
      zoneSlots.erase(zoneSlots.begin() + static_cast<std::ptrdiff_t>(nQubits), zoneSlots.end());
      constexpr std::size_t lookaheadWindow = 1; //must be 1 for this IPS
      std::vector<std::pair<qc::Qubit, double>> qubitWeightList = getQubitWeightList(twoQubitGateLayers);
      auto qubitPairWeightList = getQubitPairWeightList(twoQubitGateLayers, lookaheadWindow);
      const auto firstStages = getFirstStagesPair(lookaheadWindow, twoQubitGateLayers);

      std::ranges::sort(qubitWeightList, [](const auto& a, const auto& b){return a.second > b.second;});
      std::ranges::sort(qubitPairWeightList, [](const auto& a, const auto& b){return a.second > b.second;});

      assignments = getAssignmentsBasicPair(zoneSlots, qubitPairWeightList, qubitWeightList);
      break;
    }

    case InitialPlacementStrategy::VIPER: {
      printf("Initialising Placement Zone...\n");
      PlacementZone placement_zone = initPlacementZone(nQubits);
      printf("Assigning Qubits to Blocks...\n");
      getQubitBlockAssignment(placement_zone, twoQubitGateLayers, 5);
      printf("Initialising Assigning Blocks to Columns...\n");
      assignments = getAssignmentsColBlock(nQubits, placement_zone, 20);
      printPlacementZoneColumns(placement_zone);
      break;
    }

    case InitialPlacementStrategy::COBRA: {
      zoneSlots.erase(zoneSlots.begin() + static_cast<std::ptrdiff_t>(nQubits), zoneSlots.end());
      std::vector<std::pair<qc::Qubit, double>> qubitWeightList = getQubitWeightList(twoQubitGateLayers);
      auto qubitPairWeightList = getQubitPairWeightList(twoQubitGateLayers, twoQubitGateLayers.size());
      const auto& [maxWeightQubit, qubitWeight] = std::ranges::max(qubitWeightList, [](const auto& a, const auto& b){return a.second > b.second;});
    }

    default: {break;}
  }

  printf("Testing Placement...");
  testOutput(assignments, nQubits);
  for (std::size_t i = 0; i < nQubits; ++i) {
    const auto qc = static_cast<qc::Qubit>(i);

    initialPlacement.push_back(assignments.at(qc));
  }

  return initialPlacement;
}

auto HeuristicPlacer::makeIntermediatePlacement(
    const Placement& previousPlacement,
    const std::unordered_set<qc::Qubit>& previousReuseQubits,
    const std::unordered_set<qc::Qubit>& reuseQubits,
    const TwoQubitGateLayer& twoQubitGates,
    const TwoQubitGateLayer& nextTwoQubitGates) const
    -> std::pair<Placement, Placement> {
  const auto& gatePlacement = placeGatesInEntanglementZone(
      previousPlacement, previousReuseQubits, twoQubitGates, reuseQubits,
      nextTwoQubitGates);
  SPDLOG_TRACE("Placed gates");
  return {gatePlacement,
          placeAtomsInStorageZone(gatePlacement, reuseQubits, twoQubitGates,
                                  nextTwoQubitGates)};
}

auto HeuristicPlacer::addGateOption(
    const RowColumnMap<uint8_t>& discreteTargetRows,
    const RowColumnMap<uint8_t>& discreteTargetColumns, const SLM& leftSLM,
    const size_t leftRow, const size_t leftCol, const SLM& rightSLM,
    const size_t rightRow, const size_t rightCol, const SLM& nearestSLM,
    const size_t r, const size_t c, GateJob& job) const -> void {
  //                  other
  //         ┌─┐       ┌─┐ <-- Entanglement sites
  //         └┬┘       └┬┘
  //          │╲dis2   ╱│
  //     dis1 │  ╲   ╱  │
  //          │    ╳    │
  //          │  ╱   ╲  │ dis4
  //          │╱dis3   ╲│
  //         ┌┴┐       ┌┴┐ <-- Storage sites
  //         └─┘       └─┘
  //          ^         ^
  //        atom1     atom2
  const auto& [otherSLM, otherRow, otherCol] =
      architecture_.get().otherEntanglementSite(nearestSLM, r, c);
  const auto dis1 = static_cast<float>(architecture_.get().distance(
      leftSLM, leftRow, leftCol, nearestSLM, r, c));
  const auto dis2 = static_cast<float>(architecture_.get().distance(
      rightSLM, rightRow, rightCol, nearestSLM, r, c));
  const auto dis3 = static_cast<float>(architecture_.get().distance(
      leftSLM, leftRow, leftCol, otherSLM, otherRow, otherCol));
  const auto dis4 = static_cast<float>(architecture_.get().distance(
      rightSLM, rightRow, rightCol, otherSLM, otherRow, otherCol));
  if (dis1 + dis4 <= dis2 + dis3) {
    job.options.emplace_back(GateJob::Option{
        std::array{
            std::array{
                discreteTargetRows.at(std::pair{std::cref(nearestSLM), r}),
                discreteTargetColumns.at(std::pair{std::cref(nearestSLM), c})},
            std::array{
                discreteTargetRows.at(std::pair{std::cref(otherSLM), otherRow}),
                discreteTargetColumns.at(
                    std::pair{std::cref(otherSLM), otherCol})}},
        std::array{dis1, dis4}});
  } else {
    job.options.emplace_back(GateJob::Option{
        std::array{
            std::array{
                discreteTargetRows.at(std::pair{std::cref(otherSLM), otherRow}),
                discreteTargetColumns.at(
                    std::pair{std::cref(otherSLM), otherCol})},
            std::array{
                discreteTargetRows.at(std::pair{std::cref(nearestSLM), r}),
                discreteTargetColumns.at(std::pair{std::cref(nearestSLM), c})}},
        std::array{dis2, dis3}});
  }
}

auto HeuristicPlacer::placeGatesInEntanglementZone(
    const Placement& previousPlacement,
    const std::unordered_set<qc::Qubit>& reuseQubits,
    const TwoQubitGateLayer& twoQubitGates,
    const std::unordered_set<qc::Qubit>& nextReuseQubits,
    const TwoQubitGateLayer& nextTwoQubitGates) const -> Placement {
  // Duplicate the previous placement as a starting point for the current
  auto currentPlacement = previousPlacement;
  //===------------------------------------------------------------------===//
  // Find gates and atoms that must be placed
  //===------------------------------------------------------------------===//
  std::set<std::pair<double, QubitPair>, std::greater<>> gatesToPlace;
  std::vector<qc::Qubit> atomsToPlace;
  for (const auto& gate : twoQubitGates) {
    const auto& [first, second] = gate;
    if (const auto firstQubitReuse =
            reuseQubits.contains(first) &&
            std::get<0>(previousPlacement[first]).get().isEntanglement();
        !firstQubitReuse &&
        (!reuseQubits.contains(second) ||
         std::get<0>(previousPlacement[second]).get().isStorage())) {
      const auto& [storageSLM1, storageRow1, storageCol1] =
          previousPlacement[first];
      const auto& [storageSLM2, storageRow2, storageCol2] =
          previousPlacement[second];
      const auto& [nearestSLM, nearestRow, nearestCol] =
          architecture_.get().nearestEntanglementSite(storageSLM1, storageRow1,
                                                      storageCol1, storageSLM2,
                                                      storageRow2, storageCol2);
      const auto& [otherSLM, otherRow, otherCol] =
          architecture_.get().otherEntanglementSite(nearestSLM, nearestRow,
                                                    nearestCol);
      // Calculate the various distances to the entanglement sites
      //
      // Example:
      //       nearest    other
      //         ┌─┐       ┌─┐ <-- Entanglement sites
      //         └┬┘       └┬┘
      //          │╲dis2   ╱│
      //     dis1 │  ╲   ╱  │
      //          │    ╳    │
      //          │  ╱   ╲  │ dis4
      //          │╱dis3   ╲│
      //         ┌┴┐       ┌┴┐ <-- Storage sites
      //         └─┘       └─┘
      //          ^         ^
      //   gate->first   gate->second
      const auto dis1 =
          architecture_.get().distance(storageSLM1, storageRow1, storageCol1,
                                       nearestSLM, nearestRow, nearestCol);
      const auto dis2 =
          architecture_.get().distance(storageSLM2, storageRow2, storageCol2,
                                       nearestSLM, nearestRow, nearestCol);
      const auto dis3 = architecture_.get().distance(
          storageSLM1, storageRow1, storageCol1, otherSLM, otherRow, otherCol);
      const auto dis4 = architecture_.get().distance(
          storageSLM2, storageRow2, storageCol2, otherSLM, otherRow, otherCol);
      if (dis1 + dis4 <= dis2 + dis3) {
        // if the situation is as depicted in the example above
        gatesToPlace.emplace(std::max(dis1, dis4), gate);
      } else {
        // otherwise, either the entanglement sites or storage sites are
        // flipped
        gatesToPlace.emplace(std::max(dis2, dis3), gate);
      }
      atomsToPlace.emplace_back(first);
      atomsToPlace.emplace_back(second);
    } else {
      // one of the qubits is reused, so no need to place a gate
      if (firstQubitReuse) {
        const auto& [slm, r, c] = previousPlacement[first];
        currentPlacement[second] =
            architecture_.get().otherEntanglementSite(slm, r, c);
      } else {
        // the second qubit is reused
        const auto& [slm, r, c] = previousPlacement[second];
        currentPlacement[first] =
            architecture_.get().otherEntanglementSite(slm, r, c);
      }
    }
  }
  if (gatesToPlace.empty()) {
    return currentPlacement;
  }
  //===------------------------------------------------------------------===//
  // Discretize the previous placement of the atoms to be placed
  //===------------------------------------------------------------------===//
  const auto& [discreteRows, discreteColumns] =
      discretizePlacementOfAtoms(previousPlacement, atomsToPlace);
  //===------------------------------------------------------------------===//
  // Extract occupied entanglement sites from the previous placement
  //===------------------------------------------------------------------===//
  // This set will only contain the first SLM in a pair of entanglement
  // SLMs and represents the occupied pair of sites
  SiteSet occupiedEntanglementSites{};
  for (const auto qubit : reuseQubits) {
    const auto& [slm, r, c] = previousPlacement[qubit];
    if (slm.get().isEntanglement()) {
      occupiedEntanglementSites.emplace(slm, r, c);
      // also add the interaction partner's site to the set
      const auto& [otherSLM, otherRow, otherCol] =
          architecture_.get().otherEntanglementSite(slm, r, c);
      occupiedEntanglementSites.emplace(otherSLM, otherRow, otherCol);
    }
  }
  //===------------------------------------------------------------------===//
  // Discretize the free sites for the atoms to be placed
  //===------------------------------------------------------------------===//
  const auto& [discreteTargetRows, discreteTargetColumns] =
      discretizeNonOccupiedEntanglementSites(occupiedEntanglementSites);
  std::unordered_map<uint8_t, std::unordered_map<uint8_t, Site>> targetSites;
  for (const auto& [row, r] : discreteTargetRows) {
    const SLM& slm = row.first.get();
    auto& targetSitesForThisRow = targetSites.try_emplace(r).first->second;
    for (const auto& [column, c] : discreteTargetColumns) {
      if (&slm == &column.first.get()) {
        targetSitesForThisRow.emplace(c,
                                      std::tie(slm, row.second, column.second));
      }
    }
  }
  assert(!targetSites.empty());
  //===------------------------------------------------------------------===//
  // Initialize the gate jobs
  //===------------------------------------------------------------------===//
  /// The number of either atom or gate jobs that must be performed
  const size_t nJobs = gatesToPlace.size();
  /**
   * a list of all gates that must be placed in the entanglement zone before a
   * rydberg layer
   */
  std::vector<GateJob> gateJobs;
  gateJobs.reserve(nJobs);
  for (const auto& gate : gatesToPlace | std::views::values) {
    const auto& [leftAtom, rightAtom] = gate;
    const auto& [leftSLM, leftRow, leftCol] = previousPlacement[leftAtom];
    const auto& [rightSLM, rightRow, rightCol] = previousPlacement[rightAtom];
    const auto& [nearestSLM, nearestRow, nearestCol] =
        architecture_.get().nearestEntanglementSite(
            leftSLM, leftRow, leftCol, rightSLM, rightRow, rightCol);
    auto& job = gateJobs.emplace_back();
    job.qubits = gate;
    job.currentSites = std::array{
        std::array{discreteRows.at(std::pair{std::cref(leftSLM), leftRow}),
                   discreteColumns.at(std::pair{std::cref(leftSLM), leftCol})},
        std::array{
            discreteRows.at(std::pair{std::cref(rightSLM), rightRow}),
            discreteColumns.at(std::pair{std::cref(rightSLM), rightCol})}};
    size_t rLow = 0;
    size_t rHigh = nearestSLM.get().nRows;
    size_t cLow = 0;
    size_t cHigh = nearestSLM.get().nCols;
    if (config_.useWindow) {
      rLow = nearestRow > windowMinHeight_ / 2
                 ? nearestRow - (windowMinHeight_ / 2)
                 : 0;
      rHigh = std::min(nearestRow + (windowMinHeight_ / 2) + 1,
                       nearestSLM.get().nRows);
      cLow = nearestCol > config_.windowMinWidth / 2
                 ? nearestCol - (config_.windowMinWidth / 2)
                 : 0;
      cHigh = std::min(nearestCol + (config_.windowMinWidth / 2) + 1,
                       nearestSLM.get().nCols);
    }
    for (size_t r = rLow; r < rHigh; ++r) {
      for (size_t c = cLow; c < cHigh; ++c) {
        if (!occupiedEntanglementSites.contains(std::tie(nearestSLM, r, c))) {
          addGateOption(discreteTargetRows, discreteTargetColumns, leftSLM,
                        leftRow, leftCol, rightSLM, rightRow, rightCol,
                        nearestSLM, r, c, job);
        }
      }
    }
    size_t expansion = 0;
    while (config_.useWindow &&
           static_cast<double>(job.options.size()) <
               config_.windowShare * static_cast<double>(nJobs)) {
      // the window does not contain enough options, so expand it
      ++expansion;
      size_t windowWidth = 0;
      size_t windowHeight = 0;
      if (config_.windowRatio < 1.0) {
        // landscapes ==> expand width and adjust height
        windowWidth = config_.windowMinWidth + expansion;
        windowHeight = static_cast<size_t>(
            std::round(config_.windowRatio * static_cast<double>(windowWidth)));
      } else {
        // portraits ==> expand height and adjust width
        windowHeight = windowMinHeight_ + expansion;
        windowWidth = static_cast<size_t>(std::round(
            static_cast<double>(windowHeight) / config_.windowRatio));
      }
      auto rLowNew =
          nearestRow > windowHeight / 2 ? nearestRow - (windowHeight / 2) : 0;
      auto rHighNew =
          std::min(nearestRow + (windowHeight / 2) + 1, nearestSLM.get().nRows);
      auto cLowNew =
          nearestCol > windowWidth / 2 ? nearestCol - (windowWidth / 2) : 0;
      auto cHighNew =
          std::min(nearestCol + (windowWidth / 2) + 1, nearestSLM.get().nCols);
      if (rLowNew < rLow) {
        assert(rLow - rLowNew == 1);
        for (size_t c = cLowNew; c < cHighNew; ++c) {
          if (!occupiedEntanglementSites.contains(
                  std::tie(nearestSLM, rLowNew, c))) {
            addGateOption(discreteTargetRows, discreteTargetColumns, leftSLM,
                          leftRow, leftCol, rightSLM, rightRow, rightCol,
                          nearestSLM, rLowNew, c, job);
          }
        }
      }
      if (rHighNew > rHigh) {
        assert(rHighNew - rHigh == 1);
        for (size_t c = cLowNew; c < cHighNew; ++c) {
          // NOTE: we have to use rHighNew - 1 here, which is equal to rHigh
          if (!occupiedEntanglementSites.contains(
                  std::tie(nearestSLM, rHigh, c))) {
            addGateOption(discreteTargetRows, discreteTargetColumns, leftSLM,
                          leftRow, leftCol, rightSLM, rightRow, rightCol,
                          nearestSLM, rHigh, c, job);
          }
        }
      }
      if (cLowNew < cLow) {
        assert(cLow - cLowNew == 1);
        for (size_t r = rLow; r < rHigh; ++r) {
          if (!occupiedEntanglementSites.contains(
                  std::tie(nearestSLM, r, cLowNew))) {
            addGateOption(discreteTargetRows, discreteTargetColumns, leftSLM,
                          leftRow, leftCol, rightSLM, rightRow, rightCol,
                          nearestSLM, r, cLowNew, job);
          }
        }
      }
      if (cHighNew > cHigh) {
        assert(cHighNew - cHigh == 1);
        for (size_t r = rLow; r < rHigh; ++r) {
          // NOTE: we have to use cHighNew - 1 here, which is equal to cHigh
          if (!occupiedEntanglementSites.contains(
                  std::tie(nearestSLM, r, cHigh))) {
            addGateOption(discreteTargetRows, discreteTargetColumns, leftSLM,
                          leftRow, leftCol, rightSLM, rightRow, rightCol,
                          nearestSLM, r, cHigh, job);
          }
        }
      }
      rLow = rLowNew;
      rHigh = rHighNew;
      cLow = cLowNew;
      cHigh = cHighNew;
    }
    std::ranges::sort(
        job.options,
        [](const GateJob::Option& lhs, const GateJob::Option& rhs) -> bool {
          return lhs.distance < rhs.distance;
        });
    // Determine whether a lookahead for the gate should be considered.
    // That is the case if the gate to be placed contains a reuse qubit
    // because then we do not only decide the position of the gate in this layer
    // but also of the gate in the next layer.
    bool leftReuse = nextReuseQubits.contains(leftAtom);
    bool rightReuse = nextReuseQubits.contains(rightAtom);
    if (leftReuse || rightReuse) {
      qc::Qubit nextInteractionPartner = 0;
      for (const auto& nextGate : nextTwoQubitGates) {
        const auto& [nextLeftAtom, nextRightAtom] = nextGate;
        if (leftReuse) {
          if (nextLeftAtom == leftAtom) {
            nextInteractionPartner = nextRightAtom;
            break;
          }
          if (nextRightAtom == leftAtom) {
            nextInteractionPartner = nextLeftAtom;
            break;
          }
        }
        if (rightReuse) {
          if (nextLeftAtom == rightAtom) {
            nextInteractionPartner = nextRightAtom;
            break;
          }
          if (nextRightAtom == rightAtom) {
            nextInteractionPartner = nextLeftAtom;
            break;
          }
        }
      }
      job.meanLookaheadCost = 0.0F;
      const auto& [nextSLM, nextRow, nextCol] =
          previousPlacement[nextInteractionPartner];
      for (auto& option : job.options) {
        const auto& [row, col] = option.sites[leftReuse ? 0 : 1];
        const auto& [targetSLM, targetRow, targetCol] =
            targetSites.at(row).at(col);
        const auto distance = static_cast<float>(architecture_.get().distance(
            nextSLM, nextRow, nextCol, targetSLM, targetRow, targetCol));
        option.lookaheadCost = config_.lookaheadFactor * std::sqrt(distance);
        job.meanLookaheadCost += option.lookaheadCost;
      }
      job.meanLookaheadCost /= static_cast<float>(job.options.size());
    }
  }
  //===------------------------------------------------------------------===//
  // Get the extent of discrete source and target
  //===------------------------------------------------------------------===//
  assert(!discreteRows.empty()); // ==> the following std::max_element does not
                                 // return a nullptr
  const uint8_t maxDiscreteSourceRow =
      std::ranges::max_element(discreteRows, [](const auto& lhs,
                                                const auto& rhs) {
        return lhs.second < rhs.second;
      })->second;
  assert(
      !discreteColumns.empty()); // ==> the following std::max_element does not
  // return a nullptr
  const uint8_t maxDiscreteSourceColumn =
      std::ranges::max_element(discreteColumns, [](const auto& lhs,
                                                   const auto& rhs) {
        return lhs.second < rhs.second;
      })->second;
  assert(!discreteTargetRows
              .empty()); // ==> the following std::max_element does not
  // return a nullptr
  const uint8_t maxDiscreteTargetRow =
      std::ranges::max_element(discreteTargetRows, [](const auto& lhs,
                                                      const auto& rhs) {
        return lhs.second < rhs.second;
      })->second;
  assert(!discreteTargetColumns
              .empty()); // ==> the following std::max_element does not
  // return a nullptr
  const uint8_t maxDiscreteTargetColumn =
      std::ranges::max_element(discreteTargetColumns, [](const auto& lhs,
                                                         const auto& rhs) {
        return lhs.second < rhs.second;
      })->second;
  // return a nullptr
  const std::array<float, 2> scaleFactors{
      std::min(1.F, static_cast<float>(1 + maxDiscreteTargetRow) /
                        static_cast<float>(1 + maxDiscreteSourceRow)),
      std::min(1.F, static_cast<float>(1 + maxDiscreteTargetColumn) /
                        static_cast<float>(1 + maxDiscreteSourceColumn))};
  //===------------------------------------------------------------------===//
  // Run the DFS algorithm
  //===------------------------------------------------------------------===//
  const auto deepeningFactor = config_.deepeningFactor;
  const auto deepeningValue = config_.deepeningValue;
  std::shared_ptr<const GateNode> node;
  switch (config_.method) {
  case Config::Method::IDS:
    node = iterativeDivingSearch<GateNode>(
        std::make_shared<const GateNode>(),
        [&gateJobs](const auto& node) { return getNeighbors(gateJobs, node); },
        [nJobs](const auto& node) { return isGoal(nJobs, node); },
        [](const auto& node) { return getCost(node); },
        [&gateJobs, deepeningFactor, deepeningValue,
         &scaleFactors](const auto& node) {
          return getHeuristic(gateJobs, deepeningFactor, deepeningValue,
                              scaleFactors, node);
        },
        config_.trials, config_.queueCapacity);
    break;
  case Config::Method::ASTAR:
    node = aStarTreeSearch<GateNode>(
        std::make_shared<const GateNode>(),
        [&gateJobs](const auto& node) { return getNeighbors(gateJobs, node); },
        [nJobs](const auto& node) { return isGoal(nJobs, node); },
        [](const auto& node) { return getCost(node); },
        [&gateJobs, deepeningFactor, deepeningValue,
         &scaleFactors](const auto& node) {
          return getHeuristic(gateJobs, deepeningFactor, deepeningValue,
                              scaleFactors, node);
        },
        config_.maxNodes);
    break;
  default:
    qc::unreachable();
  }
  //===------------------------------------------------------------------===//
  // Extract the final mapping
  //===------------------------------------------------------------------===//
  while (node->parent != nullptr) {
    assert(0 < node->level);
    assert(node->level <= gateJobs.size());
    const auto& job = gateJobs[node->level - 1];
    const auto& option = job.options[node->option];
    for (size_t j = 0; j < 2; ++j) {
      const auto atom = job.qubits[j];
      const auto& [row, col] = option.sites[j];
      currentPlacement[atom] = targetSites.at(row).at(col);
    }
    node = node->parent;
  }
  return currentPlacement;
}

auto HeuristicPlacer::placeAtomsInStorageZone(
    const Placement& previousPlacement,
    const std::unordered_set<qc::Qubit>& reuseQubits,
    const TwoQubitGateLayer& twoQubitGates,
    const TwoQubitGateLayer& nextTwoQubitGates) const -> Placement {
  // Duplicate the previous placement as a starting point for the current
  auto currentPlacement = previousPlacement;
  if (twoQubitGates.empty()) {
    return currentPlacement;
  }
  //===------------------------------------------------------------------===//
  // Find atoms that must be placed
  //===------------------------------------------------------------------===//
  std::vector<qc::Qubit> atomsToPlace;
  double maxDistance = 0.0;
  size_t atomIndexWithMaxDistance = 0;
  for (const auto& gate : twoQubitGates) {
    for (const auto qubit : gate) {
      const auto& [slm, r, c] = previousPlacement[qubit];
      const auto& [nearestSLM, nearestRow, nearestCol] =
          architecture_.get().nearestStorageSite(slm, r, c);
      const auto distance = architecture_.get().distance(
          slm, r, c, nearestSLM, nearestRow, nearestCol);
      if (distance > maxDistance) {
        maxDistance = distance;
        atomIndexWithMaxDistance = atomsToPlace.size();
      }
      atomsToPlace.emplace_back(qubit);
    }
  }
  //===------------------------------------------------------------------===//
  // Discretize the previous placement of the atoms to be placed
  //===------------------------------------------------------------------===//
  // Place the atoms with the longest distance first, but the place atoms in
  // increasing order of distance to the first atom
  // swap atoms at index 0 and atomIndexWithMaxDistance
  std::swap(atomsToPlace.front(), atomsToPlace.at(atomIndexWithMaxDistance));
  const auto& frontPlacement = previousPlacement[atomsToPlace.front()];
  std::set<std::pair<double, qc::Qubit>> atomsWithoutFirstAtom;
  for (auto atomIt = atomsToPlace.cbegin() + 1; atomIt != atomsToPlace.cend();
       ++atomIt) {
    const auto& [frontSLM, frontRow, frontCol] = frontPlacement;
    const auto& [slm, r, c] = previousPlacement[*atomIt];
    const auto distance =
        architecture_.get().distance(slm, r, c, frontSLM, frontRow, frontCol);
    atomsWithoutFirstAtom.emplace(distance, *atomIt);
  }
  std::ranges::transform(std::as_const(atomsWithoutFirstAtom),
                         atomsToPlace.begin() + 1,
                         [](const auto& pair) { return pair.second; });
  // Discretize the previous placement of the atoms to be placed that are
  // ordered now
  const auto& [discreteRows, discreteColumns] =
      discretizePlacementOfAtoms(previousPlacement, atomsToPlace);
  //===------------------------------------------------------------------===//
  // Extract occupied storage sites from the previous placement
  //===------------------------------------------------------------------===//
  SiteSet occupiedStorageSites{};
  for (const auto& [slm, r, c] : previousPlacement) {
    if (slm.get().isStorage()) {
      occupiedStorageSites.emplace(slm, r, c);
    }
  }
  //===------------------------------------------------------------------===//
  // Discretize the free sites for the atoms to be placed
  //===------------------------------------------------------------------===//
  const auto& [discreteTargetRows, discreteTargetColumns] =
      discretizeNonOccupiedStorageSites(occupiedStorageSites);
  std::unordered_map<uint8_t, std::unordered_map<uint8_t, Site>> targetSites;
  for (const auto& [row, r] : discreteTargetRows) {
    const SLM& slm = row.first.get();
    auto& targetSitesForThisRow = targetSites.try_emplace(r).first->second;
    for (const auto& [column, c] : discreteTargetColumns) {
      if (&slm == &column.first.get()) {
        targetSitesForThisRow.emplace(c,
                                      std::tie(slm, row.second, column.second));
      }
    }
  }
  assert(!targetSites.empty());
  //===------------------------------------------------------------------===//
  // Initialize atom jobs
  //===------------------------------------------------------------------===//
  /// The number of either atom or gate jobs that must be performed
  size_t nJobs = atomsToPlace.size();
  /**
   * a list of all atoms that must be placed in the storage zone after a
   * rydberg layer
   */
  std::vector<AtomJob> atomJobs;
  atomJobs.reserve(nJobs);
  uint8_t minDiscreteColumnOfNearestSite = std::numeric_limits<uint8_t>::max();
  uint8_t maxDiscreteColumnOfNearestSite = 0;
  for (const auto atom : atomsToPlace) {
    const auto& [previousSLM, previousRow, previousCol] =
        previousPlacement[atom];
    const auto& [nearestSLM, nearestRow, nearestCol] =
        architecture_.get().nearestStorageSite(previousSLM, previousRow,
                                               previousCol);
    uint8_t discreteColumnOfNearestSite = 0;
    const auto& it = discreteTargetColumns.find(
        std::pair{std::cref(nearestSLM), nearestCol});
    if (it != discreteTargetColumns.end()) {
      discreteColumnOfNearestSite = it->second;
    } else {
      // distance of the actual (occupied) nearest column and a column with free
      // sites
      int64_t minDistance = std::numeric_limits<int64_t>::max();
      for (const auto& [SlmColumn, discreteColumn] : discreteTargetColumns) {
        const auto& [slm, column] = SlmColumn;
        if (slm.get() == nearestSLM.get()) {
          if (const auto currentDistance =
                  std::abs(static_cast<int64_t>(column) -
                           static_cast<int64_t>(nearestCol)) < minDistance) {
            discreteColumnOfNearestSite = discreteColumn;
            minDistance = currentDistance;
          }
        }
      }
    }
    minDiscreteColumnOfNearestSite =
        std::min(minDiscreteColumnOfNearestSite, discreteColumnOfNearestSite);
    maxDiscreteColumnOfNearestSite =
        std::max(maxDiscreteColumnOfNearestSite, discreteColumnOfNearestSite);
    auto& job = atomJobs.emplace_back();
    job.atom = atom;
    job.currentSite = std::array{
        discreteRows.at(std::pair{std::cref(previousSLM), previousRow}),
        discreteColumns.at(std::pair{std::cref(previousSLM), previousCol})};
    if (reuseQubits.contains(atom)) {
      // atom can be reused, so we add an option for the atom to stay at the
      // current site
      job.options.emplace_back(AtomJob::Option{{0, 0}, true, 0.0F});
    }
    size_t rLow = 0;
    size_t rHigh = nearestSLM.get().nRows;
    size_t cLow = 0;
    size_t cHigh = nearestSLM.get().nCols;
    if (config_.useWindow) {
      rLow = nearestRow > windowMinHeight_ / 2
                 ? nearestRow - (windowMinHeight_ / 2)
                 : 0;
      rHigh = std::min(nearestRow + (windowMinHeight_ / 2) + 1,
                       nearestSLM.get().nRows);
      cLow = nearestCol > config_.windowMinWidth / 2
                 ? nearestCol - (config_.windowMinWidth / 2)
                 : 0;
      cHigh = std::min(nearestCol + (config_.windowMinWidth / 2) + 1,
                       nearestSLM.get().nCols);
    }
    for (size_t r = rLow; r < rHigh; ++r) {
      for (size_t c = cLow; c < cHigh; ++c) {
        if (!occupiedStorageSites.contains(std::tie(nearestSLM, r, c))) {
          const auto distance = static_cast<float>(architecture_.get().distance(
              previousSLM, previousRow, previousCol, nearestSLM, r, c));
          job.options.emplace_back(AtomJob::Option{
              {discreteTargetRows.at(std::pair{std::cref(nearestSLM), r}),
               discreteTargetColumns.at(std::pair{std::cref(nearestSLM), c})},
              false,
              distance});
        }
      }
    }
    size_t expansion = 0;
    while (config_.useWindow &&
           static_cast<double>(job.options.size()) <
               config_.windowShare * static_cast<double>(nJobs)) {
      // the window does not contain enough options, so expand it
      ++expansion;
      size_t windowWidth = 0;
      size_t windowHeight = 0;
      if (config_.windowRatio < 1.0) {
        // landscape ==> expand width and adjust height
        // the overall width and height is divided by 2 later, hence an
        // expansion of 2 is needed to actually increase the window size
        windowWidth = config_.windowMinWidth + 2 * expansion;
        windowHeight = static_cast<size_t>(
            std::round(config_.windowRatio * static_cast<double>(windowWidth)));
      } else {
        // portraits ==> expand height and adjust width
        // the overall width and height is divided by 2 later, hence an
        // expansion of 2 is needed to actually increase the window size
        windowHeight = windowMinHeight_ + 2 * expansion;
        windowWidth = static_cast<size_t>(std::round(
            static_cast<double>(windowHeight) / config_.windowRatio));
      }
      auto rLowNew =
          nearestRow > windowHeight / 2 ? nearestRow - (windowHeight / 2) : 0;
      auto rHighNew =
          std::min(nearestRow + (windowHeight / 2) + 1, nearestSLM.get().nRows);
      auto cLowNew =
          nearestCol > windowWidth / 2 ? nearestCol - (windowWidth / 2) : 0;
      auto cHighNew =
          std::min(nearestCol + (windowWidth / 2) + 1, nearestSLM.get().nCols);
      if (rLowNew < rLow) {
        assert(rLow - rLowNew == 1);
        for (size_t c = cLowNew; c < cHighNew; ++c) {
          if (!occupiedStorageSites.contains(
                  std::tie(nearestSLM, rLowNew, c))) {
            const auto distance =
                static_cast<float>(architecture_.get().distance(
                    previousSLM, previousRow, previousCol, nearestSLM, rLowNew,
                    c));
            job.options.emplace_back(AtomJob::Option{
                {discreteTargetRows.at(
                     std::pair{std::cref(nearestSLM), rLowNew}),
                 discreteTargetColumns.at(std::pair{std::cref(nearestSLM), c})},
                false,
                distance});
          }
        }
      }
      if (rHighNew > rHigh) {
        assert(rHighNew - rHigh == 1);
        for (size_t c = cLowNew; c < cHighNew; ++c) {
          // NOTE: we have to use rHighNew - 1 here, which is equal to rHigh
          if (!occupiedStorageSites.contains(std::tie(nearestSLM, rHigh, c))) {
            const auto distance =
                static_cast<float>(architecture_.get().distance(
                    previousSLM, previousRow, previousCol, nearestSLM, rHigh,
                    c));
            job.options.emplace_back(AtomJob::Option{
                {discreteTargetRows.at(std::pair{std::cref(nearestSLM), rHigh}),
                 discreteTargetColumns.at(std::pair{std::cref(nearestSLM), c})},
                false,
                distance});
          }
        }
      }
      if (cLowNew < cLow) {
        assert(cLow - cLowNew == 1);
        for (size_t r = rLow; r < rHigh; ++r) {
          if (!occupiedStorageSites.contains(
                  std::tie(nearestSLM, r, cLowNew))) {
            const auto distance =
                static_cast<float>(architecture_.get().distance(
                    previousSLM, previousRow, previousCol, nearestSLM, r,
                    cLowNew));
            job.options.emplace_back(AtomJob::Option{
                {discreteTargetRows.at(std::pair{std::cref(nearestSLM), r}),
                 discreteTargetColumns.at(
                     std::pair{std::cref(nearestSLM), cLowNew})},
                false,
                distance});
          }
        }
      }
      if (cHighNew > cHigh) {
        assert(cHighNew - cHigh == 1);
        for (size_t r = rLow; r < rHigh; ++r) {
          // NOTE: we have to use cHighNew - 1 here, which is equal to cHigh
          if (!occupiedStorageSites.contains(std::tie(nearestSLM, r, cHigh))) {
            const auto distance =
                static_cast<float>(architecture_.get().distance(
                    previousSLM, previousRow, previousCol, nearestSLM, r,
                    cHigh));
            job.options.emplace_back(AtomJob::Option{
                {discreteTargetRows.at(std::pair{std::cref(nearestSLM), r}),
                 discreteTargetColumns.at(
                     std::pair{std::cref(nearestSLM), cHigh})},
                false,
                distance});
          }
        }
      }
      rLow = rLowNew;
      rHigh = rHighNew;
      cLow = cLowNew;
      cHigh = cHighNew;
    }
    std::ranges::sort(
        job.options,
        [](const AtomJob::Option& lhs, const AtomJob::Option& rhs) -> bool {
          return lhs.distance < rhs.distance;
        });
    // Determine lookahead for the atom to be placed is a reuse qubit
    for (const auto& nextGate : nextTwoQubitGates) {
      const auto& [nextLeftAtom, nextRightAtom] = nextGate;
      if (nextLeftAtom == atom || nextRightAtom == atom) {
        const qc::Qubit nextInteractionPartner =
            nextLeftAtom == atom ? nextRightAtom : nextLeftAtom;
        const auto& [nextSLM, nextRow, nextCol] =
            previousPlacement[nextInteractionPartner];
        // use this attribute to sum up the costs and calculate the mean by
        // dividing by the number of options
        job.meanLookaheadCost = 0.0F;
        for (auto& option : job.options) {
          if (option.reuse) {
            const auto distance =
                static_cast<float>(architecture_.get().distance(
                    nextSLM, nextRow, nextCol, previousSLM, previousRow,
                    previousCol));
            // NOTE: the multiplication with the lookahead factor is missing
            // here on purpose as this is the distance the "reuse" costs by
            // taking the distance of the next interaction partner as cost here
            option.lookaheadCost =
                std::max(0.0F, std::sqrt(distance) - config_.reuseLevel);
          } else {
            const auto& [row, col] = option.site;
            const auto& [targetSLM, targetRow, targetCol] =
                targetSites.at(row).at(col);
            const auto distance = static_cast<float>(
                architecture_.get().distance(nextSLM, nextRow, nextCol,
                                             targetSLM, targetRow, targetCol));
            option.lookaheadCost =
                config_.lookaheadFactor * std::sqrt(distance);
          }
          job.meanLookaheadCost += option.lookaheadCost;
        }
        job.meanLookaheadCost /= static_cast<float>(job.options.size());
        break;
      }
    }
  }
  //===------------------------------------------------------------------===//
  // Get the extent of discrete source and target
  //===------------------------------------------------------------------===//
  assert(!discreteRows.empty()); // ==> the following std::max_element does not
                                 // return a nullptr
  const uint8_t maxDiscreteSourceRow =
      std::ranges::max_element(discreteRows, [](const auto& lhs,
                                                const auto& rhs) {
        return lhs.second < rhs.second;
      })->second;
  assert(
      !discreteColumns.empty()); // ==> the following std::max_element does not
  // return a nullptr
  const uint8_t maxDiscreteSourceColumn =
      std::ranges::max_element(discreteColumns, [](const auto& lhs,
                                                   const auto& rhs) {
        return lhs.second < rhs.second;
      })->second;
  assert(!discreteTargetRows
              .empty()); // ==> the following std::max_element does not
  // return a nullptr
  const uint8_t maxDiscreteTargetRow =
      std::ranges::max_element(discreteTargetRows, [](const auto& lhs,
                                                      const auto& rhs) {
        return lhs.second < rhs.second;
      })->second;
  assert(!discreteTargetColumns
              .empty()); // ==> the following std::max_element does not
  // return a nullptr
  const uint8_t maxDiscreteTargetColumn =
      std::ranges::max_element(discreteTargetColumns, [](const auto& lhs,
                                                         const auto& rhs) {
        return lhs.second < rhs.second;
      })->second;
  const std::array<float, 2> scaleFactors{
      std::min(1.F, static_cast<float>(1 + maxDiscreteTargetRow) /
                        static_cast<float>(1 + maxDiscreteSourceRow)),
      std::min(
          std::max(1.F, static_cast<float>(maxDiscreteColumnOfNearestSite -
                                           minDiscreteColumnOfNearestSite) /
                            static_cast<float>(maxDiscreteSourceColumn)),
          static_cast<float>(1 + maxDiscreteTargetColumn) /
              static_cast<float>(1 + maxDiscreteSourceColumn))};
  //===------------------------------------------------------------------===//
  // Run the DFS algorithm
  //===------------------------------------------------------------------===//
  /**
   * A list of all nodes that have been created so far.
   * This happens when a node is expanded by calling getNeighbors.
   */
  const auto deepeningFactor = config_.deepeningFactor;
  const auto deepeningValue = config_.deepeningValue;

  std::shared_ptr<const AtomNode> node;
  switch (config_.method) {
  case Config::Method::IDS:
    node = iterativeDivingSearch<AtomNode>(
        std::make_shared<const AtomNode>(),
        [&atomJobs](const auto& node) { return getNeighbors(atomJobs, node); },
        [nJobs](const auto& node) { return isGoal(nJobs, node); },
        [](const auto& node) { return getCost(node); },
        [&atomJobs, deepeningFactor, deepeningValue,
         &scaleFactors](const auto& node) {
          return getHeuristic(atomJobs, deepeningFactor, deepeningValue,
                              scaleFactors, node);
        },
        config_.trials, config_.queueCapacity);
    break;
  case Config::Method::ASTAR:
    node = aStarTreeSearch<AtomNode>(
        std::make_shared<const AtomNode>(),
        [&atomJobs](const auto& node) { return getNeighbors(atomJobs, node); },
        [nJobs](const auto& node) { return isGoal(nJobs, node); },
        [](const auto& node) { return getCost(node); },
        [&atomJobs, deepeningFactor, deepeningValue,
         &scaleFactors](const auto& node) {
          return getHeuristic(atomJobs, deepeningFactor, deepeningValue,
                              scaleFactors, node);
        },
        config_.maxNodes);
    break;
  default:
    qc::unreachable();
  }
  //===------------------------------------------------------------------===//
  // Extract the final mapping
  //===------------------------------------------------------------------===//
  while (node->parent != nullptr) {
    assert(0 < node->level);
    assert(node->level <= atomJobs.size());
    const auto& job = atomJobs[node->level - 1];
    assert(node->option < job.options.size());
    const auto& option = job.options[node->option];
    if (!option.reuse) {
      const auto atom = job.atom;
      const auto& [row, col] = option.site;
      currentPlacement[atom] = targetSites.at(row).at(col);
    }
    node = node->parent;
  }
  return currentPlacement;
}

auto HeuristicPlacer::getCost(const GateNode& node) -> float {
  float cost = node.lookaheadCost;
  for (const auto d : node.maxDistancesOfPlacedAtomsPerGroup) {
    cost += std::sqrt(d);
  }
  return cost;
}

auto HeuristicPlacer::getCost(const AtomNode& node) -> float {
  float cost = node.lookaheadCost;
  for (const auto d : node.maxDistancesOfPlacedAtomsPerGroup) {
    cost += std::sqrt(d);
  }
  return cost;
}

auto HeuristicPlacer::sumStdDeviationForGroups(
    const std::array<float, 2>& scaleFactors,
    const std::vector<CompatibilityGroup>& groups) -> float {
  float sumStdDev = 0.F;
  for (const auto& groupPair : groups) {
    for (size_t i = 0; i < 2; ++i) {
      const auto& group = groupPair[i];
      std::vector<float> diffs;
      diffs.reserve(groups.size());
      for (const auto& [key, value] : group) {
        diffs.emplace_back(static_cast<float>(value) -
                           (scaleFactors[i] * static_cast<float>(key)));
      }
      const auto n = static_cast<float>(group.size());
      const auto mean = std::accumulate(diffs.cbegin(), diffs.cend(), 0.F) / n;
      const auto variance =
          std::accumulate(diffs.cbegin(), diffs.cend(), 0.F,
                          [mean](const float acc, const float diff) -> float {
                            return acc + ((diff - mean) * (diff - mean));
                          }) /
          n;
      sumStdDev += std::sqrt(variance);
    }
  }
  return sumStdDev;
}

auto HeuristicPlacer::getHeuristic(const std::vector<AtomJob>& atomJobs,
                                   const float deepeningFactor,
                                   const float deepeningValue,
                                   const std::array<float, 2>& scaleFactors,
                                   const AtomNode& node) -> float {
  const auto nAtomJobs = atomJobs.size();
  const auto nUnplacedAtoms = static_cast<float>(nAtomJobs - node.level);
  float maxDistanceOfUnplacedAtom = 0.0F;
  float accMinLookaheadCost = 0.0F;
  for (size_t i = node.level; i < nAtomJobs; ++i) {
    const auto& job = atomJobs[i];
    accMinLookaheadCost += job.meanLookaheadCost;
    for (const auto& option : job.options) {
      if (option.reuse) {
        // since them jobs are sorted by distance and this has distance 0 as the
        // atom does not need to be moved at all, this job will always be the
        // first one for atoms that may be reused
        break;
      }
      if (!node.consumedFreeSites.contains(option.site)) {
        // this assumes that the first found free site is the nearest free site
        // for that atom. This requires that the job options are sorted by
        // distance.
        maxDistanceOfUnplacedAtom =
            std::max(maxDistanceOfUnplacedAtom, option.distance);
        break; // exit when the first free site is found
      }
    }
  }
  float maxDistanceOfPlacedAtom = 0.0F;
  for (const auto d : node.maxDistancesOfPlacedAtomsPerGroup) {
    maxDistanceOfPlacedAtom = std::max(maxDistanceOfPlacedAtom, d);
  }
  float heuristic = maxDistanceOfUnplacedAtom <= maxDistanceOfPlacedAtom
                        ? 0.F
                        : std::sqrt(maxDistanceOfUnplacedAtom) -
                              std::sqrt(maxDistanceOfPlacedAtom);
  heuristic += accMinLookaheadCost;
  heuristic +=
      deepeningFactor *
      (sumStdDeviationForGroups(scaleFactors, node.groups) + deepeningValue) *
      nUnplacedAtoms;
  return heuristic;
}

auto HeuristicPlacer::getHeuristic(const std::vector<GateJob>& gateJobs,
                                   const float deepeningFactor,
                                   const float deepeningValue,
                                   const std::array<float, 2>& scaleFactors,
                                   const GateNode& node) -> float {
  const auto nGateJobs = gateJobs.size();
  const auto nUnplacedGates = static_cast<float>(nGateJobs - node.level);
  float maxDistanceOfUnplacedAtom = 0.0F;
  float accMeanLookaheadCost = 0.0F;
  for (size_t i = node.level; i < nGateJobs; ++i) {
    const auto& job = gateJobs[i];
    accMeanLookaheadCost += job.meanLookaheadCost;
    for (const auto& option : job.options) {
      // this assumes that the first found pair of free sites is the nearest
      // pair of free sites for that gate. This requires that the job options
      // are sorted by distance.
      if (std::ranges::all_of(option.sites,
                              [&node](const DiscreteSite& site) -> bool {
                                return !node.consumedFreeSites.contains(site);
                              })) {
        maxDistanceOfUnplacedAtom =
            std::max(maxDistanceOfUnplacedAtom,
                     *std::ranges::max_element(option.distance));
        break; // exit when the first free site pair is found
      }
    }
  }
  float maxDistanceOfPlacedAtom = 0.0F;
  for (const auto d : node.maxDistancesOfPlacedAtomsPerGroup) {
    maxDistanceOfPlacedAtom = std::max(maxDistanceOfPlacedAtom, d);
  }
  float heuristic = maxDistanceOfUnplacedAtom <= maxDistanceOfPlacedAtom
                        ? 0.F
                        : std::sqrt(maxDistanceOfUnplacedAtom) -
                              std::sqrt(maxDistanceOfPlacedAtom);
  heuristic += accMeanLookaheadCost;
  heuristic +=
      deepeningFactor *
      (sumStdDeviationForGroups(scaleFactors, node.groups) + deepeningValue) *
      nUnplacedGates;
  return heuristic;
}

auto HeuristicPlacer::getNeighbors(const std::vector<AtomJob>& atomJobs,
                                   const std::shared_ptr<const AtomNode>& node)
    -> std::vector<std::shared_ptr<const AtomNode>> {
  const size_t atomToBePlacedNext = node->level;
  assert(atomToBePlacedNext < atomJobs.size());
  const auto& atomJob = atomJobs[atomToBePlacedNext];
  std::vector<std::shared_ptr<const AtomNode>> neighbors;
  assert(atomJob.options.size() <= std::numeric_limits<uint16_t>::max());
  for (uint16_t i = 0; i < static_cast<uint16_t>(atomJob.options.size()); ++i) {
    const auto& option = atomJob.options[i];
    const auto& [site, reuse, distance, lookaheadCost] = option;
    // skip the sites that are already consumed
    if (!reuse && node->consumedFreeSites.contains(site)) {
      continue;
    }
    // make a copy of the node, the parent of the child
    AtomNode child = *node; // make a copy
    child.parent = node;
    if (!reuse) {
      child.consumedFreeSites.emplace(site);
      // check whether the current placement is compatible with any existing
      // group
      checkCompatibilityAndAddPlacement(
          atomJob.currentSite.front(), site.front(), atomJob.currentSite.back(),
          site.back(), distance, child.groups,
          child.maxDistancesOfPlacedAtomsPerGroup);
    }
    child.option = i;
    ++child.level;
    child.lookaheadCost += lookaheadCost;
    // add the child to the list of children to be returned
    neighbors.emplace_back(std::make_shared<const AtomNode>(std::move(child)));
  }
  return neighbors;
}

auto HeuristicPlacer::getNeighbors(const std::vector<GateJob>& gateJobs,
                                   const std::shared_ptr<const GateNode>& node)
    -> std::vector<std::shared_ptr<const GateNode>> {
  const size_t gateToBePlacedNext = node->level;
  assert(gateToBePlacedNext < gateJobs.size());
  const auto& gateJob = gateJobs[gateToBePlacedNext];
  std::vector<std::shared_ptr<const GateNode>> neighbors;
  // Get the current placement of the atoms that must be placed next
  const auto& [currentSiteOfLeftAtom, currentSiteOfRightAtom] =
      gateJob.currentSites;
  assert(gateJob.options.size() <= std::numeric_limits<uint16_t>::max());
  for (uint16_t i = 0; i < static_cast<uint16_t>(gateJob.options.size()); ++i) {
    const auto& option = gateJob.options[i];
    const auto& [sites, distances, lookaheadCost] = option;
    const auto& [leftSite, rightSite] = sites;
    // skip if one of the sites is already consumed
    if (node->consumedFreeSites.contains(leftSite) ||
        node->consumedFreeSites.contains(rightSite)) {
      continue;
    }
    // make a copy of the node, the parent of the child as use this as a
    // starting point for the new node
    GateNode child = *node; // make a copy
    child.parent = node;
    ++child.level;
    child.option = i;
    child.consumedFreeSites.emplace(leftSite);
    child.consumedFreeSites.emplace(rightSite);
    // check whether the current placement is compatible with any existing
    // horizontal group
    checkCompatibilityAndAddPlacement(
        currentSiteOfLeftAtom.front(), leftSite.front(),
        currentSiteOfLeftAtom.back(), leftSite.back(), distances.front(),
        child.groups, child.maxDistancesOfPlacedAtomsPerGroup);
    checkCompatibilityAndAddPlacement(
        currentSiteOfRightAtom.front(), rightSite.front(),
        currentSiteOfRightAtom.back(), rightSite.back(), distances.back(),
        child.groups, child.maxDistancesOfPlacedAtomsPerGroup);
    child.lookaheadCost += lookaheadCost;
    // add the final child to the list of children to be returned
    neighbors.emplace_back(std::make_shared<const GateNode>(std::move(child)));
  }
  return neighbors;
}

auto HeuristicPlacer::checkCompatibilityWithGroup(
    const uint8_t key, const uint8_t value,
    const std::map<uint8_t, uint8_t>& group)
    -> std::optional<
        std::pair<std::map<uint8_t, uint8_t>::const_iterator, bool>> {
  if (auto it = group.lower_bound(key); it != group.end()) {
    // an assignment for this key already exists in this group
    if (const auto& [upperKey, upperValue] = *it; upperKey == key) {
      if (upperValue == value) {
        // the new placement is compatible with this group and key already
        // exists
        return std::pair{it, true};
      }
    } else {
      // if (upperKey > key)
      if (it != group.begin()) {
        // it can be safely decremented
        if (const auto lowerValue = std::prev(it)->second;
            lowerValue < value && value < upperValue) {
          // the new placement is compatible with this group
          return std::pair{it, false};
        }
      } else {
        // if (it == hGroup.begin())
        if (value < upperValue) {
          // the new placement is compatible with this group
          return std::pair{it, false};
        }
      }
    }
  } else {
    // if (it == hGroup.end())
    // it can be safely decremented because the group must contain
    // at least one element
    if (const auto lowerValue = std::prev(it)->second; lowerValue < value) {
      // the new placement is compatible with this group
      return std::pair{it, false};
    }
  }
  return std::nullopt;
}

auto HeuristicPlacer::checkCompatibilityAndAddPlacement(
    const uint8_t hKey, const uint8_t hValue, const uint8_t vKey,
    const uint8_t vValue, const float distance,
    std::vector<CompatibilityGroup>& groups, std::vector<float>& maxDistances)
    -> bool {
  size_t i = 0;
  for (auto& group : groups) {
    auto& [hGroup, vGroup] = group;
    if (const auto& hCompatible =
            checkCompatibilityWithGroup(hKey, hValue, hGroup);
        hCompatible) {
      if (const auto& vCompatible =
              checkCompatibilityWithGroup(vKey, vValue, vGroup)) {
        const auto& [hIt, hExists] = *hCompatible;
        const auto& [vIt, vExists] = *vCompatible;
        // the new placement is compatible with this group
        if (!hExists) {
          hGroup.emplace_hint(hIt, hKey, hValue);
        }
        if (!vExists) {
          vGroup.emplace_hint(vIt, vKey, vValue);
        }
        maxDistances[i] = std::max(maxDistances[i], distance);
        return true;
      }
    }
    ++i;
  }
  // no compatible group could be found, and a new group is created
  auto& [hGroup, vGroup] = groups.emplace_back();
  hGroup.emplace(hKey, hValue);
  vGroup.emplace(vKey, vValue);
  maxDistances.emplace_back(distance);
  return false;
}

auto HeuristicPlacer::Config::createWithAStarDefaults() -> Config {
  Config config;
  config.method = Method::ASTAR;
  config.deepeningFactor = 0.8F;
  config.deepeningValue = 0.2F;
  config.lookaheadFactor = 0.2F;
  return config;
}
auto HeuristicPlacer::Config::createWithIDSDefaults() -> Config {
  Config config;
  config.method = Method::IDS;
  config.deepeningFactor = 0.01F;
  config.deepeningValue = 0.0F;
  config.lookaheadFactor = 0.4F;
  return config;
}
auto HeuristicPlacer::Config::createForMethod(const Method method) -> Config {
  return method == Method::ASTAR ? createWithAStarDefaults()
                                 : createWithIDSDefaults();
}
HeuristicPlacer::HeuristicPlacer(const Architecture& architecture,
                                 const Config& config)
    : architecture_(architecture), config_(config) {
  // get first storage SLM and first entanglement SLM
  const auto& firstStorageSLM = *architecture_.get().storageZones.front();
  const auto& firstEntanglementSLM =
      architecture_.get().entanglementZones.front()->front();
  // check which side of the first storage SLM is closer to the entanglement
  // SLM
  if (firstStorageSLM.location.second < firstEntanglementSLM.location.second) {
    // if the entanglement SLM is closer to the last row of the storage, SLM
    // start initial placement of the atoms in the last row instead of the
    // first and hence revert initial placement
    reverseInitialPlacement_ = true;
  }
  windowMinHeight_ = static_cast<size_t>(std::round(
      config_.windowRatio * static_cast<double>(config_.windowMinWidth)));
}

auto HeuristicPlacer::place(
    const size_t nQubits,
    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers,
    const std::vector<std::unordered_set<qc::Qubit>>& reuseQubits)
    -> std::vector<Placement> {
  std::vector<Placement> placement;
  placement.reserve((2 * twoQubitGateLayers.size()) + 1);

  if (config_.initialPlacementStrategy == Config::InitialPlacementStrategy::Default) {
    placement.emplace_back(makeInitialPlacement(nQubits));
  } else {
    placement.emplace_back(makeInitialPlacementMulti(nQubits, twoQubitGateLayers, config_.initialPlacementStrategy));
  }
  for (size_t layer = 0; layer < twoQubitGateLayers.size(); ++layer) {
    const auto& [gatePlacement, qubitPlacement] = makeIntermediatePlacement(
        placement.back(),
        layer == 0 ? std::unordered_set<qc::Qubit>{} : reuseQubits[layer - 1],
        layer == reuseQubits.size() ? std::unordered_set<qc::Qubit>{}
                                    : reuseQubits[layer],
        twoQubitGateLayers[layer],
        layer == twoQubitGateLayers.size() - 1 ? TwoQubitGateLayer{}
                                               : twoQubitGateLayers[layer + 1]);
    placement.emplace_back(gatePlacement);
    placement.emplace_back(qubitPlacement);
    SPDLOG_DEBUG("Placed layer: {}", layer);
  }
  return placement;
}
} // namespace na::zoned
