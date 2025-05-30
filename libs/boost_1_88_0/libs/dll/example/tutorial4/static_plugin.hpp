// Copyright 2014 Renato Tegon Forti, Antony Polukhin.
// Copyright Antony Polukhin, 2015-2025.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

//[plugcpp_my_plugin_static
#include <boost/dll/alias.hpp>                          // for BOOST_DLL_ALIAS
#include <memory>
#include "../tutorial_common/my_plugin_api.hpp"

namespace my_namespace {
    std::shared_ptr<my_plugin_api> create_plugin();   // Forward declaration
} // namespace my_namespace

BOOST_DLL_ALIAS(
    my_namespace::create_plugin,                        // <-- this function is exported with...
    create_plugin                                       // <-- ...this alias name
)
//]

