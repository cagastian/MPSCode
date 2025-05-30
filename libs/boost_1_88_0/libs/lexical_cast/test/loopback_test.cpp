//  Unit test for boost::lexical_cast.
//
//  See http://www.boost.org for most recent version, including documentation.
//
//  Copyright Alexander Nasonov, 2006.
//  Copyright Antony Polukhin, 2023-2025.
//
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt).
//
//  Test round-tripping conversion FPT -> string -> FPT,
//  where FPT is Floating Point Type.

#include <boost/lexical_cast.hpp>

#include <boost/core/lightweight_test.hpp>

#if (defined(__CYGWIN__) || defined(__FreeBSD__) || defined(__NetBSD__) \
   || (defined(__hppa) && !defined(__OpenBSD__)) || (defined(__NO_LONG_DOUBLE_MATH) && (DBL_MANT_DIG != LDBL_MANT_DIG))) \
   || defined(__MINGW64__)
#  define BOOST_LEXICAL_CAST_NO_LONG_DOUBLE_MATH_FUNCTIONS
#endif

template<class T>
void test_round_conversion()
{
    T epsilon = std::numeric_limits<T>::epsilon();
    std::string const epsilon_s = boost::lexical_cast<std::string>(epsilon);
    BOOST_TEST(epsilon == boost::lexical_cast<T>(epsilon_s));

    T max_ = (std::numeric_limits<T>::max)();
    std::string const max_s = boost::lexical_cast<std::string>(max_);
    BOOST_TEST(max_ == boost::lexical_cast<T>(max_s));

    T min_ = (std::numeric_limits<T>::min)();
    std::string const min_s = boost::lexical_cast<std::string>(min_);
    BOOST_TEST(min_ == boost::lexical_cast<T>(min_s));

    T max_div137 = max_ / 137;
    std::string max_div137_s = boost::lexical_cast<std::string>(max_div137);
    BOOST_TEST(max_div137 == boost::lexical_cast<T>(max_div137_s));

    T epsilon_mult137 = epsilon * 137;
    std::string epsilon_mult137_s(boost::lexical_cast<std::string>(epsilon_mult137));
    BOOST_TEST(epsilon_mult137 == boost::lexical_cast<T>(epsilon_mult137_s));

}

// See bug http://tinyurl.com/vhpvo
template<class T>
void test_msvc_magic_values()
{
    T magic_msvc = 0.00010000433948393407;
    std::string magic_msvc_s = boost::lexical_cast<std::string>(magic_msvc);
    BOOST_TEST(magic_msvc == boost::lexical_cast<T>(magic_msvc_s));
}

void test_round_conversion_float()
{
    test_round_conversion<float>();
}

void test_round_conversion_double()
{
    test_round_conversion<double>();
    test_msvc_magic_values<double>();
}

void test_round_conversion_long_double()
{
// We do not run tests on compilers and Standard Libraries with poor support of long double
#if !defined(BOOST_LEXICAL_CAST_NO_LONG_DOUBLE_MATH_FUNCTIONS)
    test_round_conversion<long double>();
    test_msvc_magic_values<long double>();
#endif
    BOOST_TEST(true);
}

int main()
{
    test_round_conversion_float();
    test_round_conversion_double();
    test_round_conversion_long_double();

    return boost::report_errors();
}

