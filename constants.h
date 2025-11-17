#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cstddef>
#include <numbers>

constexpr double h = 0.02; // comprimento médio da malha
constexpr size_t degree = 1;
constexpr bool export_solution = 1;
constexpr bool export_domain = 1;
constexpr size_t num_ellipses = 10; //será uma variável dinâmica depois provavelmente

constexpr double pi = std::numbers::pi;

#endif