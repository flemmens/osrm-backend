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
#include "../Util/StringUtil.h"

#include <cstdlib>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

template <class DataFacadeT> class MapMatchingPlugin : public BasePlugin
{
  private:
    std::shared_ptr<SearchEngine<DataFacadeT>> search_engine_ptr;

  public:
    MapMatchingPlugin(DataFacadeT *facade) : descriptor_string("match"), facade(facade)
    {
        search_engine_ptr = std::make_shared<SearchEngine<DataFacadeT>>(facade);
    }

    virtual ~MapMatchingPlugin() { search_engine_ptr.reset(); }

    const std::string GetDescriptor() const final { return descriptor_string; }

    void HandleRequest(const RouteParameters &route_parameters, http::Reply &reply) final
    {
        // check number of parameters

        SimpleLogger().Write() << "1";
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
        SimpleLogger().Write() << "2";

        RawRouteData raw_route;
        Matching::CandidateLists candidate_lists;
        candidate_lists.resize(route_parameters.coordinates.size());

        SimpleLogger().Write() << "3";
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

            while (candidate_lists[current_coordinate].size() < 10)
            {
                // TODO: add dummy candidates, if any are missing
                // TODO: add factory method to get an invalid PhantomNode/Distance pair
            }
        }
        SimpleLogger().Write() << "4";

        // call the actual map matching
        search_engine_ptr->map_matching(10, route_parameters.coordinates, candidate_lists, raw_route);

        if (INVALID_EDGE_WEIGHT == raw_route.shortest_path_length)
        {
            SimpleLogger().Write(logDEBUG) << "Error occurred, single path not found";
        }
        reply.status = http::Reply::ok;

        return;
    }

  private:
    std::string descriptor_string;
    DataFacadeT *facade;
};

#endif /* MAP_MATCHING_PLUGIN_H */
