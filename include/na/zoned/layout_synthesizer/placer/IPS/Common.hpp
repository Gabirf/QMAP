#pragma once
#include "na/zoned/Types.hpp"

#include <utility>

namespace na::zoned::ips::common {

/**
 * @brief Represents an atom shuttle movement in grid coordinates.
 */
struct Movement {
  std::size_t srcRow, dstRow;
  std::size_t srcCol, dstCol;
};

/**
 * @brief Computes the total row and column counts across all storage zones.
 *
 * Sums the rows and columns of all storage SLMs in the architecture.
 *
 * @param architecture Target zoned neutral-atom hardware architecture.
 * @return Total number of storage rows and total number of storage columns.
 */
[[nodiscard]] auto getStorageZoneSizes(const Architecture& architecture)
    -> std::pair<std::size_t, std::size_t>;

/**
 * @brief Computes all storage zone sites and their proximity to entanglement zones.
 *
 * Iterates through storage SLMs and computes the mean Euclidean distance from
 * each trap site to all sites across all entanglement SLMs. Lower distance indicates
 * higher priority for initial placement.
 *
 * @param architecture Target zoned neutral-atom hardware architecture.
 * @param reverseInitialPlacement If true, iterates rows in descending order
 *                                (starting closest to the entanglement zone).
 * @return List of (physical site, average proximity score) pairs.
 */
[[nodiscard]] auto getStorageZoneSlots(const Architecture& architecture,
                                       bool reverseInitialPlacement = false)
    -> std::vector<std::pair<Site, double>>;

/**
 * @brief Computes interaction priority weights for logical qubits.
 *
 * Weights decay inversely with layer index: 1 / (layer + 1). Qubits that interact
 * earlier and more frequently accumulate higher priority weights.
 *
 * @param twoQubitGateLayers Circuit gate layers containing two-qubit interactions.
 * @return Map from logical qubit index to total interaction weight.
 */
[[nodiscard]] auto
getQubitWeightMap(const std::vector<TwoQubitGateLayer>& twoQubitGateLayers)
    -> std::map<qc::Qubit, double>;

/**
 * @brief Converts the qubit weight map into an iterable list of pairs.
 *
 * @param twoQubitGateLayers Circuit gate layers containing two-qubit interactions.
 * @return List of (qubit, weight) pairs sorted by qubit index.
 */
[[nodiscard]] auto
getQubitWeightList(const std::vector<TwoQubitGateLayer>& twoQubitGateLayers)
    -> std::vector<std::pair<qc::Qubit, double>>;

/**
 * @brief Computes combined interaction weights for qubit pairs within a lookahead window.
 *
 * For each pair in the first @p lookaheadWindow layers, accumulates the mean
 * of their individual weights: (weight(q1) + weight(q2)) / 2.
 *
 * @param twoQubitGateLayers Circuit gate layers containing two-qubit interactions.
 * @param lookaheadWindow Maximum number of gate layers to inspect.
 * @return List of ((qubit1, qubit2), pairWeight) entries.
 */
[[nodiscard]] auto
getQubitPairWeightList(const std::vector<TwoQubitGateLayer>& twoQubitGateLayers,
                       std::size_t lookaheadWindow)
    -> std::vector<std::pair<std::pair<qc::Qubit, qc::Qubit>, double>>;

/**
 * @brief Determines if two 1D shuttling movements conflict (order violation).
 *
 * Verifies whether candidate movement k_cand -> v_cand crosses or shares
 * endpoints with movement k -> v in a single dimension (row or column).
 * In zoned AOD shuttling, parallel movements in the same dimension cannot cross paths.
 *
 * @param k_cand Start coordinate of the candidate movement.
 * @param v_cand Target coordinate of the candidate movement.
 * @param k Start coordinate of the existing movement.
 * @param v Target coordinate of the existing movement.
 * @return True if movements cross or collide; false otherwise.
 */
[[nodiscard]] auto hasConflict1D(std::size_t k_cand, std::size_t v_cand,
                                 std::size_t k, std::size_t v) -> bool;

/**
 * @brief Counts shuttling trajectory conflicts between a candidate movement and an active group.
 *
 * Evaluates row and column movement conflicts across all parallel movements
 * in @p group to detect AOD hardware shuttle collisions.
 *
 * @param move_cand Candidate atom movement from source to destination site.
 * @param group Active set of parallel atom movements.
 * @return Total number of conflicting movements in the group.
 */
[[nodiscard]] auto countConflictsWithGroup(const Movement& move_cand,
                                           const std::vector<Movement>& group)
    -> size_t;

/**
 * @brief Calculates the minimum row height needed to store all circuit qubits.
 *
 * Uses ceiling division: 1 + (nQubits - 1) / nColsTotal.
 *
 * @param nQubits Total number of logical qubits in the circuit.
 * @param architecture Target zoned neutral-atom hardware architecture.
 * @return Minimum number of storage rows required.
 */
[[nodiscard]] auto getHeightLimit(std::size_t nQubits,
                                  const Architecture& architecture)
    -> std::size_t;

/**
 * @brief Prunes storage zone columns based on a predicate if capacity allows.
 *
 * Iterates through @p zoneSlots in row-major order and invokes @p op on each
 * slot with the boolean result of @p predicate(col).
 * If remaining capacity after removing columns is insufficient for @p nQubits,
 * pruning is skipped and @p opInit is invoked instead.
 *
 * @param nQubits Total number of logical qubits to place.
 * @param nRows Number of rows in the storage zone.
 * @param nCols Number of columns in the storage zone.
 * @param zoneSlots Storage zone slots in row-major order.
 * @param op Callback invoked for each slot during pruning as `op(bool isPruned)`.
 * @param opInit Fallback callback invoked as `opInit()` when pruning is skipped.
 * @param predicate Callable returning true if a column index is pruned.
 * @return The maximum row height required to place all qubits.
 */
auto pruneColumns(std::size_t nQubits, std::size_t nRows, std::size_t nCols,
                  const std::vector<std::pair<Site, double>>& zoneSlots,
                  auto&& op, auto&& opInit, auto&& predicate)
    -> std::size_t {
  std::size_t prunedCols = 0;
  for (std::size_t c = 0; c < nCols; ++c) {
    if (predicate(c)) ++prunedCols;
  }

  const std::size_t totalSlots  = zoneSlots.size();
  const std::size_t prunedSlots = prunedCols * nRows;

  if (prunedCols > 0 && nCols > prunedCols &&
      totalSlots >= prunedSlots &&
      totalSlots - prunedSlots >= nQubits) {

    std::size_t idx = 0;
    for (const auto& slot : zoneSlots) {
      op(slot, predicate(idx));
      if (++idx == nCols) {
        idx = 0;
      }
    }

    return (nQubits - 1) / (nCols - 2) + 1;
  }
  std::forward<decltype(opInit)>(opInit)();

  return (nQubits - 1) / nCols + 1;
}

} // namespace na::zoned::ips::common
