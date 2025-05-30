//          Copyright (C) 2012, Michele Caini.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

//          Two Graphs Common Spanning Trees Algorithm
//      Based on academic article of Mint, Read and Tarjan
//     Efficient Algorithm for Common Spanning Tree Problem
// Electron. Lett., 28 April 1983, Volume 19, Issue 9, p.346-347

#include <boost/concept/requires.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/graph/two_graphs_common_spanning_trees.hpp>
#include <vector>

namespace boost
{

typedef boost::adjacency_list< boost::vecS, // OutEdgeList
    boost::vecS, // VertexList
    boost::undirectedS, // Directed
    boost::no_property, // VertexProperties
    boost::no_property, // EdgeProperties
    boost::no_property, // GraphProperties
    boost::listS // EdgeList
    >
    Graph;

typedef boost::graph_traits< Graph >::edge_descriptor edge_descriptor;

template < typename Coll, typename Seq > struct check_edge
{
public:
    BOOST_CONCEPT_ASSERT((RandomAccessContainer< Coll >));
    BOOST_CONCEPT_ASSERT((RandomAccessContainer< Seq >));

    typedef typename Coll::value_type coll_value_type;
    typedef typename Seq::value_type seq_value_type;

    BOOST_STATIC_ASSERT((is_same< coll_value_type, Seq >::value));
    BOOST_STATIC_ASSERT((is_same< seq_value_type, bool >::value));

    void operator()(Coll& coll, Seq& seq)
    {
        bool found = false;

        for (auto iterator = coll.begin(); !found && iterator != coll.end();
            ++iterator)
        {
            Seq& coll_seq = *iterator;

            BOOST_TEST(coll_seq.size() == seq.size());

            found = true;
            for (typename Seq::size_type pos = 0; found && pos < seq.size();
                ++pos)
            {
                found &= coll_seq[pos] == seq[pos];
            }
        }

        BOOST_TEST(found);
    }
};

void two_graphs_common_spanning_trees_test()
{
    Graph iG, vG;
    std::vector< edge_descriptor > iG_o { boost::add_edge(0, 1, iG).first,
        boost::add_edge(1, 3, iG).first, boost::add_edge(3, 2, iG).first,
        boost::add_edge(1, 5, iG).first, boost::add_edge(5, 4, iG).first,
        boost::add_edge(5, 6, iG).first, boost::add_edge(5, 3, iG).first,
        boost::add_edge(3, 1, iG).first, boost::add_edge(1, 3, iG).first };

    std::vector< edge_descriptor > vG_o { boost::add_edge(0, 2, vG).first,
        boost::add_edge(0, 4, vG).first, boost::add_edge(0, 5, vG).first,
        boost::add_edge(5, 1, vG).first, boost::add_edge(5, 3, vG).first,
        boost::add_edge(5, 6, vG).first, boost::add_edge(5, 4, vG).first,
        boost::add_edge(5, 2, vG).first, boost::add_edge(2, 6, vG).first };

    std::vector< std::vector< bool > > coll;
    boost::tree_collector< std::vector< std::vector< bool > >,
        std::vector< bool > >
        collector(coll);
    std::vector< bool > inL(iG_o.size(), false);

    boost::two_graphs_common_spanning_trees(iG, iG_o, vG, vG_o, collector, inL);

    check_edge< std::vector< std::vector< bool > >, std::vector< bool > >
        checker;
    std::vector< bool > check;

    check.assign({ true, true, true, true, true, true, false, false, false });
    checker(coll, check);

    check.assign({ true, true, true, true, true, true, false, false, false });
    checker(coll, check);
}

}

int main(int argc, char** argv)
{
    boost::two_graphs_common_spanning_trees_test();
    return boost::report_errors();
}
