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

#ifndef MAP_MATCHING_PLUGIN_H
#define MAP_MATCHING_PLUGIN_H

#include "plugin_base.hpp"

#include "../algorithms/object_encoder.hpp"
#include "../Util/integer_range.hpp"
#include "../data_structures/search_engine.hpp"
#include "../routing_algorithms/map_matching.hpp"
#include "../Util/simple_logger.hpp"
#include "../descriptors/descriptor_base.hpp"
#include "../descriptors/gpx_descriptor.hpp"
#include "../descriptors/json_descriptor.hpp"
#include "../Util/StringUtil.h"

#include <cstdlib>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

template <class DataFacadeT> class MapMatchingPlugin : public BasePlugin
{
  private:
    std::unordered_map<std::string, unsigned> descriptor_table;
    std::shared_ptr<SearchEngine<DataFacadeT>> search_engine_ptr;

  public:
    MapMatchingPlugin(DataFacadeT *facade) : descriptor_string("match"), facade(facade)
    {
        descriptor_table.emplace("json", 0);
        descriptor_table.emplace("gpx", 1);
        // descriptor_table.emplace("geojson", 2);
        //
        search_engine_ptr = std::make_shared<SearchEngine<DataFacadeT>>(facade);
    }

    virtual ~MapMatchingPlugin() { search_engine_ptr.reset(); }

    const std::string GetDescriptor() const final { return descriptor_string; }

    void HandleRequest(const RouteParameters &route_parameters, http::Reply &reply) final
    {
        // check number of parameters
        if (2 > route_parameters.coordinates.size() ||
            std::any_of(begin(route_parameters.coordinates),
                        end(route_parameters.coordinates),
                        [&](FixedPointCoordinate coordinate)
                        {
                return !coordinate.is_valid();
            }))
        {
            reply = http::Reply::StockReply(http::Reply::badRequest);
            return;
        }
        RawRouteData raw_route;
        Matching::CandidateLists candidate_lists;
        candidate_lists.resize(route_parameters.coordinates.size());

        // fetch  10 candidates for each given coordinate
        for (const auto current_coordinate : osrm::irange<std::size_t>(0, candidate_lists.size()))
        {
            if (!facade->IncrementalFindPhantomNodeForCoordinateWithDistance(
                    route_parameters.coordinates[current_coordinate],
                    candidate_lists[current_coordinate],
                    10))
            {
                reply = http::Reply::StockReply(http::Reply::badRequest);
                return;
            }

            BOOST_ASSERT(candidate_lists[current_coordinate].size() == 10);

            /*
            while (candidate_lists[current_coordinate].size() < 10)
            {
                // TODO: add dummy candidates, if any are missing
                // TODO: add factory method to get an invalid PhantomNode/Distance pair
            }
            */
        }

        // call the actual map matching
        std::vector<PhantomNode> matched_nodes;
        JSON::Object debug_info;
        search_engine_ptr->map_matching(10, candidate_lists, route_parameters.coordinates, matched_nodes, debug_info);

        reply.status = http::Reply::ok;

        PhantomNodes current_phantom_node_pair;
        for (unsigned i = 0; i < matched_nodes.size() - 1; ++i)
        {
            current_phantom_node_pair.source_phantom = matched_nodes[i];
            current_phantom_node_pair.target_phantom = matched_nodes[i + 1];
            raw_route.segment_end_coordinates.emplace_back(current_phantom_node_pair);
        }

        search_engine_ptr->shortest_path(
            raw_route.segment_end_coordinates, route_parameters.uturns, raw_route);

        DescriptorConfig descriptor_config;

        auto iter = descriptor_table.find(route_parameters.output_format);
        unsigned descriptor_type = (iter != descriptor_table.end() ? iter->second : 0);

        descriptor_config.zoom_level = route_parameters.zoom_level;
        descriptor_config.instructions = route_parameters.print_instructions;
        descriptor_config.geometry = route_parameters.geometry;
        descriptor_config.encode_geometry = route_parameters.compression;

        std::shared_ptr<BaseDescriptor<DataFacadeT>> descriptor;
        switch (descriptor_type)
        {
        // case 0:
        //     descriptor = std::make_shared<JSONDescriptor<DataFacadeT>>();
        //     break;
        case 1:
            descriptor = std::make_shared<GPXDescriptor<DataFacadeT>>(facade);
            break;
        // case 2:
        //      descriptor = std::make_shared<GEOJSONDescriptor<DataFacadeT>>();
        //      break;
        default:
            descriptor = std::make_shared<JSONDescriptor<DataFacadeT>>(facade);
            break;
        }

        JSON::Object result;
        descriptor->SetConfig(descriptor_config);
        descriptor->Run(raw_route, result);

        /*
         * Add some debugging information.
         */

        JSON::Array json_list;
        for (const auto& list : candidate_lists)
        {
            JSON::Array candidates;
            for (const auto candidate : list)
            {
                JSON::Array json_candidate;
                JSON::Array json_coordinates;
                PhantomNode node = candidate.first;
                json_coordinates.values.push_back(node.location.lat / COORDINATE_PRECISION);
                json_coordinates.values.push_back(node.location.lon / COORDINATE_PRECISION);
                json_candidate.values.emplace_back(json_coordinates);
                json_candidate.values.emplace_back(candidate.second);
                candidates.values.emplace_back(json_candidate);
            }
            json_list.values.emplace_back(candidates);
        }
        result.values["candidates"] = json_list;
        result.values["debug"] = debug_info;

        descriptor->Render(result, reply.content);

        return;
    }

  private:
    std::string descriptor_string;
    DataFacadeT *facade;
};

#endif /* MAP_MATCHING_PLUGIN_H */
