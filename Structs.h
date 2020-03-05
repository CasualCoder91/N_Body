/**
 * Collection of structs and enum classes.
 *
 * @author Alarich Herzner
 * @version 0.9 05.03.2020
*/

#pragma once

/** @brief sub volumes (cells) within a cell in 3D. Used by Node (octree) */
enum class Octant {
    top_left_front,
    top_right_front,
    bottom_right_front,
    bottom_left_front,
    top_left_back,
    top_right_back,
    bottom_right_back,
    bottom_left_back,
    Invalid
};