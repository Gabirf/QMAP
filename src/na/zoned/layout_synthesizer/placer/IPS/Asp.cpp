#include "na/zoned/layout_synthesizer/placer/IPS/Asp.hpp"

#include "na/zoned/Architecture.hpp"
#include "na/zoned/layout_synthesizer/placer/IPS/Common.hpp"

namespace na::zoned::ips::asp {

namespace {

/**
 * @brief Identifies the earliest gate layer in which each qubit pair interacts.
 *
 * Scans gate layers within @p lookaheadWindow and records the first layer index
 * where each two-qubit pair appears.
 *
 * @param lookaheadWindow Maximum number of gate layers to inspect.
 * @param twoQubitGateLayers Circuit gate layers containing two-qubit interactions.
 * @return Map from qubit pair to its earliest interaction layer index.
 */
/* auto getFirstStagesPair(
    const std::size_t lookaheadWindow,
    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers)
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
}*/

/**
 * @brief Maps interacting qubit pairs to adjacent physical storage traps.
 *
 * Assigns each interacting pair in @p qubitWeightPairList to consecutive physical
 * sites in @p zoneSlots. Any remaining qubits not participating in first-layer pairs
 * are sequentially assigned to the next available storage traps in descending order
 * of individual weight.
 *
 *
 * @param zoneSlots Physical storage sites ordered by proximity to the entanglement zone.
 * @param qubitWeightPairList Ranked list of interacting qubit pairs and their weights.
 * @param qubitWeightList Ranked list of individual qubits and their weights.
 * @return Mapping from logical qubits to physical storage sites.
 */
auto assignQubitPairs(
    const std::vector<std::pair<Site, double>>& zoneSlots,
    const std::vector<std::pair<std::pair<qc::Qubit, qc::Qubit>, double>>&
        qubitWeightPairList,
    const std::vector<std::pair<qc::Qubit, double>>& qubitWeightList)
    -> std::unordered_map<qc::Qubit, Site> {

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
} // namespace

/**
 * @brief Computes initial placement for zoned neutral-atom architectures using ASP.
 *
 * Places qubits interacting in the earliest gate layer into adjacent storage traps
 * closest to the entanglement zone, minimizing shuttling overhead for initial gates.
 *
 * @param nQubits Number of logical qubits in the circuit.
 * @param twoQubitGateLayers Circuit gate layers representing two-qubit interactions.
 * @param architecture Target zoned neutral-atom hardware architecture.
 * @param reverseInitialPlacement If true, seeds placement starting from the row
 *                                closest to the entanglement zone in reverse orientation.
 * @return Mapping from logical qubits to physical storage sites.
 */
auto getAssignments(const std::size_t nQubits,
                    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers,
                    const Architecture& architecture,
                    const bool reverseInitialPlacement)
    -> std::unordered_map<qc::Qubit, Site> {

  std::vector<std::pair<Site, double>> zoneSlots =
      common::getStorageZoneSlots(architecture, reverseInitialPlacement);

  std::printf("Slots: %zu, Qubits: %zu\n", zoneSlots.size(), nQubits);
  if (zoneSlots.size() < nQubits) {
    throw std::invalid_argument("Too many qubits for this architecture.\n");
  }

  zoneSlots.erase(zoneSlots.begin() + static_cast<std::ptrdiff_t>(nQubits),
                  zoneSlots.end());
  const std::size_t lookaheadWindow = 1;

  std::vector<std::pair<qc::Qubit, double>> qubitWeightList =
      common::getQubitWeightList(twoQubitGateLayers);
  auto qubitPairWeightList =
      common::getQubitPairWeightList(twoQubitGateLayers, lookaheadWindow);

  std::ranges::sort(qubitWeightList, [](const auto& a, const auto& b) {
    return a.second > b.second;
  });
  std::ranges::sort(qubitPairWeightList, [](const auto& a, const auto& b) {
    return a.second > b.second;
  });

  return assignQubitPairs(zoneSlots, qubitPairWeightList, qubitWeightList);
}
} // namespace na::zoned::ips::asp
