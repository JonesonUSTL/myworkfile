#pragma once

#include <string>

#include "fem_types.hpp"

namespace fem {

void writeVTK(const std::string& filePath, const Model& model, const Result& result, double scale = 1.0);

}  // namespace fem
