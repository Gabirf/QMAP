
#include "na/zoned/layout_synthesizer/placer/IPS/Common.hpp"

namespace na::zoned::ips::common {
auto getStorageZoneSizes(const Architecture& architecture)
    -> std::pair<std::size_t, std::size_t> {

  auto& storageZones = architecture.storageZones;
  std::size_t nRowsTotal = 0;
  std::size_t nColsTotal = 0;
  for (const auto& zonePointer : storageZones) {
    const auto& storageZone = *zonePointer;

    nRowsTotal += storageZone.nRows;
    nColsTotal += storageZone.nCols;
  }

  return std::pair{nRowsTotal, nColsTotal};
}

// TODO: Doesn't yet work for multiple storage zones (change rows and cols assigned to Site)
auto getStorageZoneSlots(const Architecture& architecture,
                         bool reverseInitialPlacement)
    -> std::vector<std::pair<Site, double>> {
  std::vector<std::pair<Site, double>> zoneSlots;
  auto& storageZones = architecture.storageZones;
  auto& entanglementZones = architecture.entanglementZones;

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

              const double distance1 =
                  architecture.distance(slm, row, col, slm1, rowEnt, colEnt);
              const double distance2 =
                  architecture.distance(slm, row, col, slm2, rowEnt, colEnt);

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

    if (reverseInitialPlacement) {
      for (std::size_t row = slm.nRows - 1;; row--) {
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

auto getQubitWeightMap(const std::vector<TwoQubitGateLayer>& twoQubitGateLayers)
    -> std::map<qc::Qubit, double> {
  std::map<qc::Qubit, double> qubitMap;
  std::size_t idx = 0;

  for (const auto& layer : twoQubitGateLayers) {
    double weight = 1.0 / (static_cast<double>(idx + 1));

    for (const auto& pair : layer) {
      qubitMap[pair[0]] += weight;
      qubitMap[pair[1]] += weight;
    }
    idx++;
  }
  return qubitMap;
}

auto getQubitWeightList(
    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers)
    -> std::vector<std::pair<qc::Qubit, double>> {
  std::map<qc::Qubit, double> qubitMap = getQubitWeightMap(twoQubitGateLayers);
  return {qubitMap.begin(), qubitMap.end()};
}

auto getQubitPairWeightList(
    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers,
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

auto hasConflict1D(std::size_t k_cand, std::size_t v_cand, std::size_t k,
                   std::size_t v) -> bool {
  if (k_cand == k) {
    return v_cand != v;
  }
  if (k_cand < k) {
    return v_cand >= v;
  }
  // k_cand > k
  return v_cand <= v;
}

auto countConflictsWithGroup(const Movement& move_cand,
                             const std::vector<Movement>& group) -> size_t {
  size_t conflicts = 0;

  for (const auto& movement : group) {

    if (hasConflict1D(move_cand.srcRow, move_cand.dstRow, movement.srcRow,
                      movement.dstRow) ||
        hasConflict1D(move_cand.srcCol, move_cand.dstCol, movement.srcCol,
                      movement.dstCol)) {
      conflicts++;
    }
  }

  return conflicts;
}

auto getHeightLimit(std::size_t nQubits, const Architecture& architecture) -> std::size_t {
  auto nColsTotal = getStorageZoneSizes(architecture).second;
  return 1 + (nQubits - 1) / nColsTotal;
}
} // namespace na::zoned::ips::common