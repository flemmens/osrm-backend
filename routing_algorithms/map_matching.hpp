/*
    open source routing machine
    Copyright (C) Dennis Luxen, others 2010

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU AFFERO General Public License as published by
the Free Software Foundation; either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
or see http://www.gnu.org/licenses/agpl.txt.
 */

#ifndef MAP_MATCHING__
#define MAP_MATCHING__

#include "routing_base.hpp"

#include "../Util/simple_logger.hpp"
#include "../Util/container.hpp"

namespace Matching
{
typedef std::vector<std::pair<PhantomNode, double>> CandidateList;
typedef std::vector<CandidateList> CandidateLists;
typedef std::pair<PhantomNodes, double> PhantomNodesWithProbability;
}

template <class DataFacadeT> class MapMatching : public BasicRoutingInterface<DataFacadeT>
{
    typedef BasicRoutingInterface<DataFacadeT> super;
    typedef typename SearchEngineData::QueryHeap QueryHeap;
    SearchEngineData &engine_working_data;

    constexpr static const double sigma_z = 4.07;

    constexpr double emission_probability(const float distance) const
    {
        return (1./(std::sqrt(2.*M_PI)*sigma_z))*std::exp(-0.5*std::pow(2,(distance/sigma_z)));
    }

// TODO: needs to be estimated from the input locations
constexpr static const double beta = 1.;
// samples/min and beta
// 1 0.49037673
// 2 0.82918373
// 3 1.24364564
// 4 1.67079581
// 5 2.00719298
// 6 2.42513007
// 7 2.81248831
// 8 3.15745473
// 9 3.52645392
// 10 4.09511775
// 11 4.67319795
// 21 12.55107715
// 12 5.41088180
// 13 6.47666590
// 14 6.29010734
// 15 7.80752112
// 16 8.09074504
// 17 8.08550528
// 18 9.09405065
// 19 11.09090603
// 20 11.87752824
// 21 12.55107715
// 22 15.82820829
// 23 17.69496773
// 24 18.07655652
// 25 19.63438911
// 26 25.40832185
// 27 23.76001877
// 28 28.43289797
// 29 32.21683062
// 30 34.56991141

    constexpr double transition_probability(const float distance, const float beta) const
    {
        return (1./beta)*std::exp(-distance/beta);
    }

    // translates a distance into how likely it is an input
    double DistanceToProbability(const float distance) const
    {
        if (0. > distance)
        {
            return 0.;
        }
        return 1. - 1. / (1. + exp((-distance + 35.) / 6.));
    }

  public:
    MapMatching(DataFacadeT *facade, SearchEngineData &engine_working_data)
        : super(facade), engine_working_data(engine_working_data)
    {
    }

    ~MapMatching() {}

    void operator()(Matching::CandidateLists &candidate_lists, RawRouteData &raw_route_data) const
    {
        SimpleLogger().Write() << "matching starts with " << candidate_lists.size() << " locations";

        // // step over adjacent pairs of candidate locations
        // {
        //     // translate distances into probabilities
        //     //  check if combined probability makes sense here
        //     // initialize forward and reverse heaps and compute all pair-wise distances

        //     // select path like this:
        //     const double normalized_length = local_distance_of_path * ( 2. -
        //     current_path_probability);
        //     //if one of the paths was not found, replace it with the other one.
        //     if(minimum_combined_distance > normalized_length) {
        //         minimum_combined_distance = normalized_length;
        //         rawRouteData.computedShortestPath.clear();
        //         super::UnpackPath(temporaryPackedPath1, rawRouteData.computedShortestPath);
        //         rawRouteData.segmentEndCoordinates.clear();
        //         rawRouteData.segmentEndCoordinates.push_back(phantomNodePair);

        //         rawRouteData.lengthOfShortestPath = local_distance_of_path;
        //     }

        // }

        // std::vector<Matching::PhantomNodesWithProbability> phantom_node_pairs;
        // for(unsigned i = 0; i < nearest_k_edges_forward.size(); ++i) {
        // 	SimpleLogger().Write() << nearest_k_edges_forward[i].first.location << " at distance
        // " << nearest_k_edges_forward[i].second << " and id " <<
        // nearest_k_edges_forward[i].first.forward_node_id;
        //     for(unsigned j = 0; j < nearest_k_edges_reverse.size(); ++j) {
        //         PhantomNodes current_phantom_node_pair;
        //         current_phantom_node_pair.source_phantom = nearest_k_edges_forward[i].first;
        //         current_phantom_node_pair.target_phantom = nearest_k_edges_reverse[j].first;
        //         SimpleLogger().Write() << "[" << i*j << "] (" <<
        //         current_phantom_node_pair.source_phantom.forward_node_id << "," <<
        //         current_phantom_node_pair.target_phantom.forward_node_id << ")";
        //         double combined_probability = nearest_k_edges_forward[i].second *
        //         nearest_k_edges_reverse[j].second;
        //         phantom_node_pairs.push_back(std::make_pair(current_phantom_node_pair,
        //         combined_probability));
        //     }
        // }

        // SimpleLogger().Write() << "Generated " << phantom_node_pairs.size() << " phantom nodes";

        // super::_queryData.InitializeOrClearFirstThreadLocalStorage();
        // super::_queryData.InitializeOrClearSecondThreadLocalStorage();
        // super::_queryData.InitializeOrClearThirdThreadLocalStorage();

        // QueryHeap & forward_heap1 = *(super::_queryData.forwardHeap);
        // QueryHeap & reverse_heap1 = *(super::_queryData.backwardHeap);

        // int minimum_combined_distance = INVALID_EDGE_WEIGHT;
        // raw_route_data.shortest_path_length = raw_route_data.alternative_path_length =
        // INVALID_EDGE_WEIGHT;

        // for (const PhantomNodesProbability & phantom_node_probality_pair : phantom_node_pairs) {
        //     const PhantomNodes & phantom_node_pair = phantom_node_probality_pair.first;
        //     const double current_path_probability = (phantom_node_probality_pair.second);
        //     forward_heap1.Clear();
        //     reverse_heap1.Clear();
        //     EdgeWeight local_distance_of_path = INVALID_EDGE_WEIGHT;

        //     NodeID middle1 = SPECIAL_NODEID;
        //     EdgeWeight min_edge_offset = 0;

        //     //insert new starting nodes into forward heap, adjusted by previous distances.
        //     if (phantom_node_pair.source_phantom.forward_node_id != SPECIAL_NODEID)
        //     {
        //         forward_heap1.Insert(
        //             phantom_node_pair.source_phantom.forward_node_id,
        //             - phantom_node_pair.source_phantom.GetForwardWeightPlusOffset(),
        //             phantom_node_pair.source_phantom.forward_node_id);
        //         min_edge_offset = std::min(min_edge_offset, -
        //         phantom_node_pair.source_phantom.GetForwardWeightPlusOffset());
        //     }

        //     if (phantom_node_pair.source_phantom.reverse_node_id != SPECIAL_NODEID)
        //     {
        //         forward_heap1.Insert(
        //             phantom_node_pair.source_phantom.reverse_node_id,
        //             - phantom_node_pair.source_phantom.GetReverseWeightPlusOffset(),
        //             phantom_node_pair.source_phantom.reverse_node_id);
        //         min_edge_offset = std::min(min_edge_offset, -
        //         phantom_node_pair.source_phantom.GetReverseWeightPlusOffset());
        //     }

        //     // insert new backward nodes into backward heap, unadjusted.
        //     if (phantom_node_pair.target_phantom.forward_node_id != SPECIAL_NODEID)
        //     {
        //         reverse_heap1.Insert(phantom_node_pair.target_phantom.forward_node_id,
        //                              phantom_node_pair.target_phantom.GetForwardWeightPlusOffset(),
        //                              phantom_node_pair.target_phantom.forward_node_id);
        //    }

        //     if (phantom_node_pair.target_phantom.reverse_node_id != SPECIAL_NODEID)
        //     {
        //         reverse_heap1.Insert(phantom_node_pair.target_phantom.reverse_node_id,
        //                              phantom_node_pair.target_phantom.GetReverseWeightPlusOffset(),
        //                              phantom_node_pair.target_phantom.reverse_node_id);
        //     }

        //     //run Dijkstra routing steps
        //     while (0 < (forward_heap1.Size() + reverse_heap1.Size()))
        //     {
        //         if (!forward_heap1.Empty())
        //         {
        //             super::RoutingStep(
        //                 forward_heap1, reverse_heap1, &middle1, &local_distance_of_path,
        //                 min_edge_offset, true);
        //         }
        //         if (!reverse_heap1.Empty())
        //         {
        //             super::RoutingStep(
        //                 reverse_heap1, forward_heap1, &middle1, &local_distance_of_path,
        //                 min_edge_offset, false);
        //         }
        //     }

        //     //No path found for both target nodes?
        //     if(INVALID_EDGE_WEIGHT == local_distance_of_path) {
        //         continue;
        //     }

        //     //Unpack path
        //     std::vector<NodeID> packed_path;
        //     super::RetrievePackedPathFromHeap(forward_heap1, reverse_heap1, middle1,
        //     packed_path);

        //     const double normalized_length = local_distance_of_path * ( 2. -
        //     current_path_probability);
        //     //if one of the paths was not found, replace it with the other one.
        //     if(minimum_combined_distance > normalized_length) {
        //         minimum_combined_distance = normalized_length;
        //         raw_route_data.computedShortestPath.clear();
        //         super::UnpackPath(packed_path, raw_route_data.computedShortestPath);
        //         raw_route_data.segmentEndCoordinates.clear();
        //         raw_route_data.segmentEndCoordinates.push_back(phantom_node_pair);

        //         raw_route_data.shortest_path_length = local_distance_of_path;
        //     }
        // }
        // SimpleLogger().Write() << "Selected path of length " <<
        // raw_route_data.shortest_path_length << " and combined " << minimum_combined_distance);
    }
};

#endif /* MAP_MATCHING__ */
