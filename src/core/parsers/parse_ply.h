#pragma once

#include "core/lajolla.h"
#include "core/math/matrix.h"
#include "geometry/shape.h"

/// Parse Stanford PLY files.
TriangleMesh parse_ply(const fs::path &filename, const Matrix4x4 &to_world);
