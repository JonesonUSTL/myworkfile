#pragma once

#include <string>

#include "fem_types.hpp"

namespace fem {

Model parseInpFile(const std::string& filePath);

}  // namespace fem
