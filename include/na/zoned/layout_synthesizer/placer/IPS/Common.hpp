#pragma once
#include "na/zoned/Types.hpp"

#include <utility>
namespace na::zoned::ips::common {

  struct Movement {
    std::size_t srcRow, dstRow;
    std::size_t srcCol, dstCol;
  };

  [[nodiscard]] auto getStorageZoneSizes(const Architecture& architecture)
      -> std::pair<std::size_t, std::size_t>;

  /**
   * @brief Computes storage sites and their average proximity to entanglement zones.
   * @return A list of pairs containing the storage site and its proximity score.
   */
  [[nodiscard]] auto getStorageZoneSlots(const Architecture& architecture, bool reverseInitialPlacement = false)
      -> std::vector<std::pair<Site, double>>;

  [[nodiscard]] auto getQubitWeightMap(
      const std::vector<TwoQubitGateLayer>& twoQubitGateLayers)
      -> std::map<qc::Qubit, double>;

  /**
   * @brief Computes the importance/weight of the qubit based on usage and earliness.
   * @return A list of pairs containing the qubit and its weight.
   */
  [[nodiscard]] auto getQubitWeightList(
      const std::vector<TwoQubitGateLayer>& twoQubitGateLayers)
      -> std::vector<std::pair<qc::Qubit, double>>;

  [[nodiscard]] auto getQubitPairWeightList(
      const std::vector<TwoQubitGateLayer>& twoQubitGateLayers ,
      std::size_t lookaheadWindow)
        -> std::vector<std::pair<std::pair<qc::Qubit, qc::Qubit>, double>>;

  [[nodiscard]] auto hasConflict1D(std::size_t k_cand, std::size_t v_cand, std::size_t k, std::size_t v) -> bool;

  [[nodiscard]] auto countConflictsWithGroup(
      const Movement& move_cand,
      const std::vector<Movement>& group) -> size_t;

  [[nodiscard]] auto getHeightLimit(std::size_t nQubits, const Architecture& architecture)
      -> std::size_t;
}

