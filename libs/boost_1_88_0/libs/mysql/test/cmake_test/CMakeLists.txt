#
# Copyright (c) 2019-2025 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

cmake_minimum_required(VERSION 3.5...3.22)

project(cmake_subdir_test LANGUAGES CXX)

option(BOOST_CI_INSTALL_TEST "Whether to build install or add_subdirectory tests" OFF)

if(BOOST_CI_INSTALL_TEST)
    find_package(boost_mysql REQUIRED)
else()
    # Generated by boostdep --brief mysql
    set(_DEPENDENCIES
        # Primary dependencies
        asio
        assert
        charconv
        compat
        config
        core
        describe
        endian
        intrusive
        mp11
        optional
        pfr
        system
        throw_exception
        variant2

        # Secondary dependencies
        align
        context
        date_time
        static_assert
        container_hash
        move
        detail
        predef
        type_traits
        utility
        winapi
        pool
        smart_ptr
        exception
        algorithm
        io
        lexical_cast
        numeric/conversion
        range
        tokenizer
        preprocessor
        array
        bind
        concept_check
        "function"
        iterator
        mpl
        regex
        tuple
        unordered
        container
        integer
        conversion
        function_types
        fusion
        functional
        typeof
    )

    # Build our dependencies, so the targets Boost::xxx are defined
    set(_BOOST_ROOT ../../../..)
    foreach(_DEPENDENCY IN LISTS _DEPENDENCIES)
        add_subdirectory(${_BOOST_ROOT}/libs/${_DEPENDENCY} boostorg/${_DEPENDENCY})
    endforeach()

    # Build our project
    add_subdirectory(${_BOOST_ROOT}/libs/mysql boostorg/mysql)
endif()

# Copied from Alexander Grund's Boost.CI
add_executable(main main.cpp)
target_link_libraries(main PRIVATE Boost::mysql)

enable_testing()
add_test(NAME main COMMAND main)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C $<CONFIG>)
