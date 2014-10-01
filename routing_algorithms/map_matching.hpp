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

#ifndef MAP_MATCHING_H
#define MAP_MATCHING_H

#include "routing_base.hpp"

#include "../Util/simple_logger.hpp"
#include "../Util/container.hpp"

#include <algorithm>
#include <iomanip>
#include <numeric>

namespace Matching
{
typedef std::vector<std::pair<PhantomNode, double>> CandidateList;
typedef std::vector<CandidateList> CandidateLists;
typedef std::pair<PhantomNodes, double> PhantomNodesWithProbability;
}

template <class DataFacadeT> class MapMatching final : public BasicRoutingInterface<DataFacadeT>
{
    typedef BasicRoutingInterface<DataFacadeT> super;
    typedef typename SearchEngineData::QueryHeap QueryHeap;
    SearchEngineData &engine_working_data;

    constexpr static const double sigma_z = 4.07;

    constexpr double emission_probability(const double distance) const
    {
        return (1. / (std::sqrt(2. * M_PI) * sigma_z)) *
               std::exp(-0.5 * std::pow((distance / sigma_z), 2.));
    }

    constexpr double log_emission_probability(const double distance) const
    {
        return std::log2(emission_probability(distance));
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

    constexpr double transition_probability(const float d_t, const float beta) const
    {
        return (1. / beta) * std::exp(-d_t / beta);
    }

    // translates a distance into how likely it is an input
    double DistanceToProbability(const double distance) const
    {
        if (0. > distance)
        {
            return 0.;
        }
        return 1. - 1. / (1. + exp((-distance + 35.) / 6.));
    }

    double compute_dt(const FixedPointCoordinate &location1,
                      const FixedPointCoordinate &location2,
                      const Matching::CandidateList &candidate_list_1,
                      const Matching::CandidateList &candidate_list_2)
    {
        // great circle distance of two locations - median/avg dist table(candidate list1/2)
        std::vector<EdgeWeight> distance_list;
        for (const auto &candidate_1 : candidate_list_1)
        {
            for (const auto &candidate_2 : candidate_list_2)
            {
                const EdgeWeight current_weight = 0; // TODO: compute path distance between the two

                distance_list.push_back(current_weight);
            }
        }


        const auto average_network_distance = std::accumulate(distance_list.begin(), distance_list.end(), 0) / distance_list.size();

        std::nth_element(distance_list.begin(), distance_list.begin() + distance_list.size()/2, distance_list.end());
        const EdgeWeight median_network_distance = distance_list[distance_list.size()/2];

        const auto great_circle_distance = FixedPointCoordinate::ApproximateDistance(location1, location2);

        // TODO: choose which one is more effective
        const auto approximated_network_distance = median_network_distance;

        if (great_circle_distance > approximated_network_distance)
        {
            return great_circle_distance - approximated_network_distance;
        }
        return approximated_network_distance - great_circle_distance;
    }

  public:
    MapMatching(DataFacadeT *facade, SearchEngineData &engine_working_data)
        : super(facade), engine_working_data(engine_working_data)
    {
    }

    void operator()(const unsigned state_size,
                    Matching::CandidateLists &timestamp_list,
                    std::vector<FixedPointCoordinate> coordinate_list,
                    RawRouteData &raw_route_data) const
    {
        BOOST_ASSERT(state_size != std::numeric_limits<unsigned>::max());
        BOOST_ASSERT(state_size != 0);
        SimpleLogger().Write() << "matching starts with " << timestamp_list.size() << " locations";

        // // step over adjacent pairs of candidate locations

        SimpleLogger().Write() << "state_size: " << state_size;

        std::vector<std::vector<double>> viterbi(state_size,
                                                 std::vector<double>(timestamp_list.size() + 1, 0));
        std::vector<std::vector<std::size_t>> parent(state_size,
                                                     std::vector<std::size_t>(timestamp_list.size() + 1, 0));

        SimpleLogger().Write() << "a";

        for (auto s = 0; s < state_size; ++s)
        {
            SimpleLogger().Write() << "initializing s: " << s << "/" << state_size;
            SimpleLogger().Write() << " distance: " << timestamp_list[0][s].second
                                   << " at " << timestamp_list[0][s].first.location
                                   << " prob " << std::setprecision(10) << emission_probability(timestamp_list[0][s].second)
                                   << " logprob " << log_emission_probability(timestamp_list[0][s].second);
            // TODO: implement
            const double emission_pr = 0.;
            viterbi[s][0] = emission_pr;
            parent[s][0] = s;
        }
        SimpleLogger().Write() << "b";

        std::vector<double> d_t_list, median_select_d_t_list;
        for (auto t = 1; t < timestamp_list.size(); ++t)
        {
            d_t_list.push_back(compute_dt(coordinate_list[t-1], coordinate_list[t],
                                          timestamp_list[t-1], timestamp_list[t]));
            median_select_d_t_list.push_back(d_t_list.back());
        }

        std::nth_element(median_select_d_t_list.begin(),
                         median_select_d_t_list.begin() + median_select_d_t_list.size()/2,
                         median_select_d_t_list.end());
        const auto median_d_t = median_select_d_t_list[median_select_d_t_list.size()/2];

        const auto beta = (1./std::log(2))*median_d_t;

        for (auto t = 1; t < timestamp_list.size(); ++t)
        {
            // compute d_t for this timestamp and the next one
            for (auto s = 0; s < state_size; ++s)
            {
                for (auto s_prime = 0; s_prime < state_size; ++s_prime)
                {
                    // TODO: implement
                    const double emission_pr = 0.;
                    // TODO: implement
                    const double transition_pr = transition_probability(beta, d_t_list[t])
                    const double new_value = viterbi[s][t] * emission_pr * transition_pr;
                    if (new_value > viterbi[s_prime][t])
                    {
                        viterbi[s_prime][t] = new_value;
                        parent[s_prime][t] = s;
                    }
                }
            }
        }
        SimpleLogger().Write() << "c";
        SimpleLogger().Write() << "timestamps: " << timestamp_list.size();
        const auto number_of_timestamps = timestamp_list.size();
        const auto max_element_iter =
            std::max_element(viterbi[number_of_timestamps].begin(), viterbi[number_of_timestamps].end());
        auto parent_index = std::distance(max_element_iter, viterbi[number_of_timestamps].begin());
        std::deque<std::size_t> reconstructed_indices;

        SimpleLogger().Write() << "d";

        for (auto i = number_of_timestamps-1; i > 0; --i)
        {
            SimpleLogger().Write() << "row: " << i << ", parent: " << parent_index;
            parent_index = parent[parent_index][i];
        }

        SimpleLogger().Write() << "e";
    }
};

#endif /* MAP_MATCHING_H */
