
#pragma once
#include "ir/Definitions.hpp"
#include "na/zoned/Types.hpp"

#include <cstddef>
#include <span>

namespace na::zoned::ips::cobra {

struct ZoneSlot {
  std::size_t id{0};
  Site site;
  double score{0.0};
  mutable std::optional<qc::Qubit> qubit{std::nullopt};

  explicit ZoneSlot(Site s, std::size_t slotId = 0, double sc = 0.0, std::optional<qc::Qubit> q = std::nullopt)
      : id(slotId), site(std::move(s)), score(sc),  qubit(q){}

  auto operator==(const ZoneSlot& other) const noexcept -> bool {
    return site == other.site;
  }

  auto operator<(const ZoneSlot& other) const noexcept -> bool {
    if (score != other.score) {
      return score < other.score;
    } else {
      return id < other.id;
    }
  }
};

struct InteractionWeightMatrix {
private:
  std::size_t n_{0};
  std::vector<double> weights_;

public:
  explicit InteractionWeightMatrix(std::size_t nQubits)
      : n_(nQubits), weights_(nQubits * nQubits, 0.0) {}

  void addInteraction(qc::Qubit q1, qc::Qubit q2, double weight) {
    weights_[q1 * n_ + q2] += weight;
    weights_[q2 * n_ + q1] += weight;
  }

  [[nodiscard]] auto getWeight(qc::Qubit q1, qc::Qubit q2) const -> double {
    return weights_[q1 * n_ + q2];
  }

  [[nodiscard]] auto getRowSplice(qc::Qubit q1) const  {
    std::vector<qc::Qubit> result;

    const auto& srcStart = weights_.begin() + q1 * n_;

    return std::span{srcStart, n_};
  }

  [[nodiscard]] auto getQubitNum() const -> std::size_t {
    return n_;
  }
};

[[nodiscard]] auto getAssignments(std::size_t nQubits,
                    const std::vector<TwoQubitGateLayer>& twoQubitGateLayers,
                    const Architecture& architecture,
                    bool reverseInitialPlacement)-> std::unordered_map<qc::Qubit, Site>;

} // namespace na::zoned::ips::cobra
