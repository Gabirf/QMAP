#include "na/zoned/Architecture.hpp"
#include "na/zoned/layout_synthesizer/placer/IPS/Asp.hpp"
#include "na/zoned/layout_synthesizer/placer/IPS/Common.hpp"
namespace na::zoned::ips::zap {

namespace {
auto getFirstStages(std::size_t nQubits, std::size_t lookaheadWindow,
                    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers)
    -> std::unordered_map<qc::Qubit, size_t> {

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

/**
 * @brief Gives the 1st each qubit was present in until phaseIdx ==
 * lookaheadWindow.
 * @param architecture
 * @param zoneSlots
 * @param qubitWeightList
 * @param firstStages
 * @param nQubits is the total number of qubits in the quantum computation
 * @return A mapping from qubit to the earliest layer it appears in.
 */
auto assignQubits(
    const Architecture& architecture,
    const std::vector<std::pair<Site, double>>& zoneSlots,
    const std::vector<std::pair<qc::Qubit, double>>& qubitWeightList,
    const std::unordered_map<qc::Qubit, std::size_t>& firstStages,
    std::size_t nQubits) -> std::unordered_map<qc::Qubit, Site> {

  using Movement = common::Movement;
  std::unordered_map zoneSlotsMap(zoneSlots.begin(), zoneSlots.end());
  std::vector<Movement> group;
  std::unordered_map<qc::Qubit, Site> results;

  for (const auto& qc : qubitWeightList | std::views::keys) {

    std::pair<std::optional<std::pair<Site, Site>>, double> minScore = {
        std::nullopt, 1e308};
    auto minScoreMove = Movement();

    for (const auto& site : zoneSlots | std::views::keys) {
      constexpr double lamda = 10.0;

      if (!zoneSlotsMap.contains(site)) {
        continue;
      }

      auto& slm = std::get<0>(site);
      auto row = std::get<1>(site);
      auto col = std::get<2>(site);

      const auto& entangleSite = architecture.nearestEntanglementSite(
          slm, row, col, slm, row, col);
      const auto& [entangleSLM, entangleRow, entangleCol] = entangleSite;

      auto& [x0, y0] = slm.get().location;
      auto rowCoord = y0 + slm.get().siteSeparation.second * row;
      auto colCoord = x0 + slm.get().siteSeparation.first * col;

      auto& [entx0, enty0] = entangleSLM.get().location;
      auto entRowCoord =
          enty0 + entangleSLM.get().siteSeparation.second * entangleRow;
      auto entColCoord =
          entx0 + entangleSLM.get().siteSeparation.first * entangleCol;
      auto movement =
          Movement(rowCoord, entRowCoord, colCoord, entColCoord);

      std::size_t count = 0;
      if (firstStages.contains(qc)) {
        count = countConflictsWithGroup(movement, group);
      }
      double distance = architecture.nearestEntanglementSiteDistance(
          slm, row, col, slm, row, col);

      double score = lamda * static_cast<double>(count) + distance;

      if (score < minScore.second) {
        minScore.first = {site, entangleSite};
        minScore.second = score;
        minScoreMove = movement;
      }
    }

    if (!minScore.first.has_value()) {
      throw std::runtime_error(
          "Error: No available storage slot found for qubit assignment!");
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

        double distance = architecture.nearestEntanglementSiteDistance(
            slm, row, col, slm, row, col);

        if (distance < minScore.second) {
          minScore.first = site;
          minScore.second = distance;
        }
      }

      const auto& site = minScore.first;

      if (!site.has_value()) {
        throw std::runtime_error(
            "Error: No available storage slot found for qubit assignment!");
      }

      zoneSlotsMap.erase(*site);
      results.emplace(qc, *site);
    }
  }

  return results;
}
} // namespace

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

  std::unordered_map<qc::Qubit, Site> assignments;

  std::ranges::sort(zoneSlots, [](const auto& a, const auto& b) {
    return a.second < b.second;
  });
  if (zoneSlots.size() > nQubits) {
    zoneSlots.erase(zoneSlots.begin() + static_cast<std::ptrdiff_t>(nQubits),
                    zoneSlots.end());
  }
  std::vector<std::pair<qc::Qubit, double>> qubitWeightList =
      common::getQubitWeightList(twoQubitGateLayers);
  constexpr std::size_t lookaheadWindow =
      2; // Can be changed but defined as 2 in ZAP
  std::unordered_map<qc::Qubit, std::size_t> firstStages =
      getFirstStages(nQubits, lookaheadWindow, twoQubitGateLayers);
  std::ranges::sort(qubitWeightList, [](const auto& a, const auto& b) {
    return a.second > b.second;
  });

  return assignQubits(architecture, zoneSlots, qubitWeightList, firstStages, nQubits);
}
} // namespace na::zoned::ips::zap
