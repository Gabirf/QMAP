
#include "na/zoned/layout_synthesizer/placer/IPS/Cobra.hpp"

#include "na/zoned/layout_synthesizer/placer/IPS/Common.hpp"

#include <format>
namespace na::zoned::ips::cobra {

namespace {

/**
 * @brief Constructs an InteractionWeightMatrix from circuit qubit pair weights.
 * 
 * @param nQubits Number of qubits in the circuit.
 * @param qubitPairWeightList List of qubit pairs and their accumulated interaction weights.
 * @return Fully initialized symmetric InteractionWeightMatrix.
 */
auto getInteractionWeightMatrix(
    const std::size_t nQubits,
    std::vector<std::pair<std::pair<qc::Qubit, qc::Qubit>, double>>&
        qubitPairWeightList) -> InteractionWeightMatrix {
  InteractionWeightMatrix interactionWeightMatrix(nQubits);

  for (const auto& [qubits, weight] : qubitPairWeightList) {
    const auto [q1, q2] = qubits;
    interactionWeightMatrix.addInteraction(q1, q2, weight);
  }

  return interactionWeightMatrix;
}

/**
 * @brief Assigns the highest-weight unassigned qubit to a given storage site.
 * 
 * Used for the initial seed placement or as a fallback when no candidate
 * neighbor has an interaction with already placed qubits.
 * 
 * @param architecture Target zoned neutral-atom architecture.
 * @param qubitWeightList List of qubits and their individual interaction weights.
 * @param slots Grid of zone slots representing storage trap locations.
 * @param bestSlot Physical site to assign the qubit to.
 * @param assigned Set tracking all currently assigned qubits.
 */
auto firstAssignment(
    const Architecture& architecture,
    const std::vector<std::pair<qc::Qubit, double>>& qubitWeightList,
    const std::vector<ZoneSlot>& slots, const Site& bestSlot,
    std::unordered_set<qc::Qubit>& assigned) {

  auto unassignedQubits =
      qubitWeightList | std::views::filter([assigned](const auto& a) {
        return !assigned.contains(a.first);
      });
  const auto maxWeightQubit =
      std::ranges::max(unassignedQubits, std::less<>{},
                       &std::pair<qc::Qubit, double>::second)
          .first;

  const auto& [slm, row, col] = bestSlot;
  const auto& [nRows, nCols] = common::getStorageZoneSizes(architecture);
  slots[row * nCols + col].qubit = maxWeightQubit;
  assigned.insert(maxWeightQubit);
}

/**
 * @brief Explores and adds valid orthogonal grid neighbors of a site to the candidate set.
 * 
 * Checks adjacent site coordinates (up, down, left, right) respecting grid boundaries
 * and the maximum row height limit, then adds unvisited sites to the search set.
 * 
 * @param architecture Target zoned neutral-atom architecture.
 * @param zoneSlots Precomputed storage zone sites and distance scores.
 * @param searchZoneSet Priority set of candidate unassigned neighbor slots.
 * @param site The reference site whose neighbors are being queried.
 * @param visitedSlots Boolean tracking matrix indicating which sites have been visited.
 * @param maxHeight Maximum allowed row height from the initial placement boundary.
 * @param reverseInitialPlacement Whether placement starts from the bottom or top row.
 * @param nextSlotIdx Counter used for deterministic tie-breaking in searchZoneSet.
 * @return Updated nextSlotIdx counter.
 */
auto addSiteNeighbours(const Architecture& architecture,
                       std::vector<std::pair<Site, double>>& zoneSlots,
                       std::set<ZoneSlot>& searchZoneSet, const Site& site,
                       std::vector<bool>& visitedSlots,
                       const std::size_t maxHeight,
                       const bool reverseInitialPlacement,
                       std::size_t nextSlotIdx = 0) -> std::size_t {

  static constexpr std::array<std::pair<int, int>, 4> displacements = {{
      {-1, 0}, // Up
      {0, -1}, // Left
      {1, 0},  // Down
      {0, 1}   // Right
  }};

  const auto& [nRows, nCols] = common::getStorageZoneSizes(architecture);
  const auto& [_, row, col] = site;

  const auto r = static_cast<int>(row);
  const auto c = static_cast<int>(col);
  const auto maxC = static_cast<int>(nCols);
  const int minRow =
      reverseInitialPlacement ? static_cast<int>(nRows - maxHeight) : 0;
  const int maxRow = reverseInitialPlacement ? static_cast<int>(nRows)
                                             : static_cast<int>(maxHeight);

  for (const auto& [dr, dc] : displacements) {
    const int nr = r + dr;
    const int nc = c + dc;

    if (nr >= minRow && nr < maxRow && nc >= 0 && nc < maxC) {

      const auto neighRow = static_cast<std::size_t>(nr);
      const auto neighCol = static_cast<std::size_t>(nc);
      const auto matrixIdx = neighRow * nCols + neighCol;
      if (matrixIdx < visitedSlots.size() && !visitedSlots[matrixIdx]) {
        const auto& [neighSite, neighScore] = zoneSlots[matrixIdx];
        searchZoneSet.emplace(neighSite, nextSlotIdx++, neighScore);
        visitedSlots[matrixIdx] = true;
      }
    }
  }

  return nextSlotIdx;
}

/**
 * @brief Selects and assigns the best-fitting qubit to the next candidate slot.
 * 
 * Extracts the highest-priority slot from the search set, evaluates candidate qubits
 * based on cumulative interaction weights with already placed neighbors, and assigns
 * the best matching qubit. If no interaction exists, falls back to firstAssignment.
 * 
 * @param architecture Target zoned neutral-atom architecture.
 * @param qubitWeightList List of qubits and their individual interaction weights.
 * @param iwMatrix Symmetric matrix containing pairwise interaction weights.
 * @param searchZoneSet Priority set of candidate unassigned slots.
 * @param slots Grid of zone slots representing storage trap locations.
 * @param assigned Set tracking all currently assigned qubits.
 * @return The physical site to which a qubit was assigned.
 */
auto assignQubitToSite(
    const Architecture& architecture,
    const std::vector<std::pair<qc::Qubit, double>>& qubitWeightList,
    const InteractionWeightMatrix& iwMatrix, std::set<ZoneSlot>& searchZoneSet,
    const std::vector<ZoneSlot>& slots,
    std::unordered_set<qc::Qubit>& assigned) {
  ZoneSlot slot = searchZoneSet.extract(searchZoneSet.begin()).value();
  const auto& site = slot.site;
  const auto& [_, row, col] = site;

  static constexpr std::array<std::pair<int, int>, 4> displacements = {{
      {-1, 0}, // Up
      {0, -1}, // Left
      {1, 0},  // Down
      {0, 1}   // Right
  }};

  const auto& [nRows, nCols] = common::getStorageZoneSizes(architecture);
  const auto r = static_cast<int>(row);
  const auto c = static_cast<int>(col);
  const auto maxR = static_cast<int>(nRows);
  const auto maxC = static_cast<int>(nCols);

  std::vector<qc::Qubit> neighbours;
  neighbours.reserve(displacements.size());

  std::vector<double> qubitScores;
  qubitScores.resize(iwMatrix.getQubitNum());

  for (const auto& [dr, dc] : displacements) {
    const int nr = r + dr;
    const int nc = c + dc;

    if (nr >= 0 && nr < maxR && nc >= 0 && nc < maxC) {
      const auto neighRow = static_cast<std::size_t>(nr);
      const auto neighCol = static_cast<std::size_t>(nc);
      const auto matrixIdx = neighRow * nCols + neighCol;
      if (matrixIdx < slots.size()) {
        if (const auto neighbour = slots[matrixIdx].qubit) {
          const auto& rowSplice = iwMatrix.getRowSplice(*neighbour);
          std::transform(qubitScores.begin(), qubitScores.end(),
                         rowSplice.begin(), qubitScores.begin(), std::plus<>{});
        }
      }
    }
  }

  int maxIdx = -1;
  double maxWeight = -1;
  for (std::size_t i = 0; i < qubitScores.size(); i++) {
    if (!assigned.contains(static_cast<qc::Qubit>(i)) &&
        qubitScores[i] > maxWeight) {
      maxIdx = i;
      maxWeight = qubitScores[i];
    }
  }

  // std::printf("%s\n", std::format("Assigned qubit {} to slot with score
  // {:.2f}", maxIdx, maxWeight).c_str());
  if (maxIdx == -1) {
    firstAssignment(architecture, qubitWeightList, slots, site, assigned);
    return site;
  }

  slots[row * nCols + col].qubit = maxIdx;
  assigned.insert(static_cast<qc::Qubit>(maxIdx));
  return site;
}
} // namespace

auto getAssignments(const std::size_t nQubits,
                    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers,
                    const Architecture& architecture,
                    const bool reverseInitialPlacement)
    -> std::unordered_map<qc::Qubit, Site> {
  std::set<ZoneSlot> searchZoneSet;
  std::vector<std::pair<Site, qc::Qubit>> assignments;

  std::vector<std::pair<Site, double>> zoneSlots =
      common::getStorageZoneSlots(architecture, false);

  // std::printf("Slots: %zu, Qubits: %zu\n", zoneSlots.size(), nQubits);
  if (zoneSlots.size() < nQubits) {
    throw std::invalid_argument("Too many qubits for this architecture.");
  }
  std::size_t nCol = common::getStorageZoneSizes(architecture).second;
  std::size_t maxHeight = (nQubits - 1) / nCol + 1;
  const auto qubitWeightList = common::getQubitWeightList(twoQubitGateLayers);
  auto qubitPairWeightList = common::getQubitPairWeightList(
      twoQubitGateLayers, twoQubitGateLayers.size());

  const auto bestSlotIt = std::ranges::min_element(
      zoneSlots, std::less<>{}, &std::pair<Site, double>::second);
  const auto bestIndex = static_cast<std::size_t>(
      std::ranges::distance(zoneSlots.begin(), bestSlotIt));
  const auto& bestSlot = bestSlotIt->first;

  InteractionWeightMatrix iwMatrix =
      getInteractionWeightMatrix(nQubits, qubitPairWeightList);

  std::vector<ZoneSlot> slots;
  slots.reserve(zoneSlots.size());
  for (const auto& [site, score] : zoneSlots) {
    slots.emplace_back(site, 0, score);
  }

  std::unordered_set<qc::Qubit> assigned;
  assigned.reserve(nQubits);
  firstAssignment(architecture, qubitWeightList, slots, bestSlot, assigned);

  std::vector<bool> visitedSlots;
  visitedSlots.resize(zoneSlots.size(), false);
  visitedSlots[bestIndex] = true;
  std::size_t nextSlotIdx =
      addSiteNeighbours(architecture, zoneSlots, searchZoneSet, bestSlot,
                        visitedSlots, maxHeight, reverseInitialPlacement);

  std::size_t i = 1;
  while (!searchZoneSet.empty() && i < nQubits) {
    // std::printf("%s\n", std::format("Run {}", i).c_str());
    const auto& site =
        assignQubitToSite(architecture, qubitWeightList, iwMatrix,
                          searchZoneSet, slots, assigned);
    // std::printf("adding neighbours\n");
    nextSlotIdx = addSiteNeighbours(architecture, zoneSlots, searchZoneSet,
                                    site, visitedSlots, maxHeight,
                                    reverseInitialPlacement, nextSlotIdx);
    i++;
  }

  std::unordered_map<qc::Qubit, Site> result;
  result.reserve(nQubits);
  for (const auto& slot : slots) {
    if (slot.qubit) {
      result.insert_or_assign(*slot.qubit, slot.site);
    }
  }

  assert(result.size() == nQubits &&
         "Failure to assign all qubits, check height condition");
  return result;
}

} // namespace na::zoned::ips::cobra
