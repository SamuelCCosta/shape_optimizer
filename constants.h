#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cstddef>
#include <numbers>

constexpr double h = 0.02; // comprimento médio da malha
constexpr size_t degree = 1;
constexpr bool hand_coded = 1;
constexpr bool export_file = 1;

constexpr double pi = std::numbers::pi;
constexpr double two_pi = 2 * std::numbers::pi;

static_assert(degree == 1 || degree == 2);
static_assert(!(degree == 2 && hand_coded));

#endif