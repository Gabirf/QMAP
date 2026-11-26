#pragma once
#include "ir/Definitions.hpp"
#include "na/zoned/Types.hpp"

#include <cstddef>
#include <unordered_map>
#include <vector>

namespace na::zoned::ips::asp {

/**
 * @brief Computes the initial placement of qubits onto storage sites using ASP.
 *
 * Places qubits interacting in the earliest gate layer into adjacent storage traps
 * closest to the entanglement zone, minimizing shuttling overhead for initial gates.
 *
 * @param nQubits Number of logical qubits in the circuit.
 * @param twoQubitGateLayers Circuit gate layers representing two-qubit interactions.
 * @param architecture Target zoned neutral-atom hardware architecture.
 *  @param reverseInitialPlacement Whether placement starts from the bottom or
 * top row.
 * @return Mapping from logical qubits to physical storage sites.
 */
[[nodiscard]] auto
getAssignments(std::size_t nQubits,
               const std::vector<TwoQubitGateLayer>& twoQubitGateLayers,
               const Architecture& architecture, bool reverseInitialPlacement)
    -> std::unordered_map<qc::Qubit, Site>;

} // namespace na::zoned::ips::asp
