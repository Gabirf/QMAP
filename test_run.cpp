#include "na/zoned/layout_synthesizer/placer/HeuristicPlacer.hpp"
#include <iostream>
#include <string>
#include <vector>

using namespace na::zoned;

static const std::string architectureJson = R"({
  "name": "a_star_placer_architecture",
  "storage_zones": [{
    "zone_id": 0,
    "slms": [{"id": 0, "site_separation": [3, 3], "r": 20, "c": 20, "location": [0, 0]}],
    "offset": [0, 0],
    "dimension": [60, 60]
  }],
  "entanglement_zones": [{
    "zone_id": 0,
    "slms": [
      {"id": 1, "site_separation": [12, 10], "r": 4, "c": 4, "location": [5, 70]},
      {"id": 2, "site_separation": [12, 10], "r": 4, "c": 4, "location": [7, 70]}
    ],
    "offset": [5, 70],
    "dimension": [50, 40]
  }],
  "aods":[{"id": 0, "site_separation": 2, "r": 20, "c": 20}],
  "rydberg_range": [[[5, 70], [55, 110]]]
})";

int main() {
  std::cout << ">>> Starting debug_placer..." << std::endl;

  const auto architecture = Architecture::fromJSONString(architectureJson);
  HeuristicPlacer::Config config;
  config.initialPlacementStrategy = HeuristicPlacer::Config::InitialPlacementStrategy::COBRA;

  HeuristicPlacer placer(architecture, config);

  constexpr size_t nQubits = 4;
  std::vector<TwoQubitGateLayer> layers = {
    {{{0U, 1U}, {2U, 3U}}}
  };

  const auto placements = placer.place(nQubits, layers, {});

  std::cout << ">>> Placement finished! Number of placements: " << placements.size() << std::endl;
  return 0;
}
