# Boost.Geometry
# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.
# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

foreach(item IN ITEMS
    box
    box_templated
    box_2d_4values
    point
    linestring
    linestring_templated
    ring
    ring_templated
    multi_point
    multi_point_templated
    multi_linestring
    multi_polygon
  )

  boost_geometry_add_example("geometries_register" ${item})
  
endforeach()
