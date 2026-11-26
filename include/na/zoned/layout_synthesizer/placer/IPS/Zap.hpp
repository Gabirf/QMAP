
#pragma once
#include "ir/Definitions.hpp"
#include "na/zoned/Types.hpp"

#include <cstddef>
#include <unordered_map>
#include <vector>

namespace na::zoned::ips::zap {
[[nodiscard]] auto getAssignments(std::size_t nQubits,
                    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers,
                    const Architecture& architecture,
                    bool reverseInitialPlacement)
    -> std::unordered_map<qc::Qubit, Site>;
}