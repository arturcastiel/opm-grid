/*
  Copyright 2026 SINTEF Digital.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <opm/grid/cpgrid/refinement/conforming/EdgeConformal.hpp>

#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/EntityRep.hpp>
#include <opm/grid/cpgrid/Geometry.hpp>
#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>
#include <opm/common/utility/SparseTable.hpp>

#include <dune/common/fvector.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <utility>
#include <tuple>
#include <vector>

namespace Opm
{
namespace Refinement
{

namespace
{

using Coord = Dune::FieldVector<double,3>;

// Within a cell, an edge a-b that one face has and its other faces reach through
// nodes near the segment (a curved fault intersection): insert those nodes into
// that face. Returns whether a face changed.
bool insertCellChainNodes(Dune::cpgrid::CpGridData& leaf, const std::vector<Coord>& P)
{
    auto& faceToPoint = GridStateWriter::faceToPoint(leaf);
    auto& cellToFace = GridStateWriter::cellToFace(leaf);
    using Edge = std::pair<int,int>;
    const auto edgeOf = [](int a, int b) { return Edge{ std::min(a, b), std::max(a, b) }; };
    const auto between = [&P](int a, int b, int m) {
        Coord d = P[b];
        d -= P[a];
        Coord am = P[m];
        am -= P[a];
        const double len2 = d.two_norm2();
        const double t = (am * d) / len2;
        if (!(t > 0.0 && t < 1.0)) {
            return false;
        }
        Coord off = d;
        off *= -t;
        off += am;
        return off.two_norm2() <= 1e-4 * len2;
    };

    std::map<int, std::map<Edge, std::vector<int>>> inserts;   // face -> edge -> nodes
    std::map<Edge, std::vector<int>> edgeFaces;
    std::map<int, std::vector<int>> lone;
    for (int cell = 0; cell < leaf.size(0); ++cell) {
        edgeFaces.clear();
        lone.clear();
        const auto row = cellToFace[Dune::cpgrid::EntityRep<0>(cell, true)];
        for (int e = 0; e < row.size(); ++e) {
            const int face = row[e].index();
            auto pts = faceToPoint[face];
            const int n = pts.size();
            for (int v = 0; v < n; ++v) {
                edgeFaces[edgeOf(pts[v], pts[(v + 1) % n])].push_back(face);
            }
        }
        for (const auto& [ab, faces] : edgeFaces) {
            if (faces.size() == 1) {
                lone[ab.first].push_back(ab.second);
                lone[ab.second].push_back(ab.first);
            }
        }
        for (const auto& [ab, faces] : edgeFaces) {
            if (faces.size() != 1) {
                continue;
            }
            const auto [a, b] = ab;
            if (P[a] == P[b]) {
                continue;
            }
            // Breadth-first from a to b over the other lone edges.
            std::map<int,int> prev{{a, a}};
            std::vector<int> queue{a};
            for (std::size_t q = 0; q < queue.size() && !prev.count(b); ++q) {
                const int u = queue[q];
                for (const int w : lone[u]) {
                    if (prev.count(w) || (u == a && w == b) || (w != b && !between(a, b, w))) {
                        continue;
                    }
                    prev[w] = u;
                    queue.push_back(w);
                }
            }
            if (!prev.count(b)) {
                continue;
            }
            std::vector<int> chain;
            for (int w = prev[b]; w != a; w = prev[w]) {
                chain.push_back(w);
            }
            inserts[faces.front()][ab] = std::move(chain);
        }
    }
    if (inserts.empty()) {
        return false;
    }

    Opm::SparseTable<int> rebuilt;
    std::vector<int> newRow;
    std::vector<std::pair<double,int>> mids;
    for (int face = 0; face < faceToPoint.size(); ++face) {
        auto row = faceToPoint[face];
        const auto it = inserts.find(face);
        if (it == inserts.end()) {
            rebuilt.appendRow(row.begin(), row.end());
            continue;
        }
        const int n = row.size();
        newRow.clear();
        for (int v = 0; v < n; ++v) {
            const int s = row[v];
            const int t = row[(v + 1) % n];
            newRow.push_back(s);
            const auto nodes = it->second.find(edgeOf(s, t));
            if (nodes == it->second.end()) {
                continue;
            }
            Coord d = P[t];
            d -= P[s];
            mids.clear();
            for (const int m : nodes->second) {
                if (std::find(row.begin(), row.end(), m) == row.end()) {
                    Coord sm = P[m];
                    sm -= P[s];
                    mids.emplace_back(sm * d, m);
                }
            }
            std::sort(mids.begin(), mids.end());
            for (const auto& [proj, m] : mids) {
                newRow.push_back(m);
            }
        }
        rebuilt.appendRow(newRow.begin(), newRow.end());
    }
    faceToPoint.swap(rebuilt);
    return true;
}

} // anonymous namespace

void edgeConformalizeLeaf(Dune::cpgrid::CpGridData& leaf)
{
    const auto& cornerGeom =
        GridStateWriter::geometry(leaf).geomVector(std::integral_constant<int,3>());
    const int numCorners = static_cast<int>(cornerGeom.size());
    if (numCorners == 0) {
        return;
    }
    auto& faceToPoint = GridStateWriter::faceToPoint(leaf);

    // Only corners some face uses: level-zero corners of refined cells stay in the
    // pool, unused, and must not land in a face.
    std::vector<char> used(numCorners, 0);
    for (int face = 0; face < faceToPoint.size(); ++face) {
        for (const int p : faceToPoint[face]) {
            used[p] = 1;
        }
    }

    std::vector<Coord> P(numCorners);
    Coord lo = cornerGeom.get(0).center();
    Coord hi = lo;
    int numUsed = 0;
    for (int i = 0; i < numCorners; ++i) {
        P[i] = cornerGeom.get(i).center();
        for (int c = 0; c < 3; ++c) {
            lo[c] = std::min(lo[c], P[i][c]);
            hi[c] = std::max(hi[c], P[i][c]);
        }
        numUsed += used[i];
    }
    double diag = 0.0;
    for (int c = 0; c < 3; ++c) {
        diag += (hi[c] - lo[c]) * (hi[c] - lo[c]);
    }
    diag = std::sqrt(diag);
    if ((diag <= 0.0) || (numUsed == 0)) {
        return;
    }
    const double tol = 1e-9 * diag;

    // Buckets sized per direction so each holds about one corner: a thin, wide model
    // gets thin, wide buckets instead of buckets as tall as it is wide.
    double spanProduct = 1.0;
    int spanned = 0;
    for (int c = 0; c < 3; ++c) {
        if (hi[c] - lo[c] > tol) {
            spanProduct *= hi[c] - lo[c];
            ++spanned;
        }
    }
    const double side = std::pow(spanProduct / numUsed, 1.0 / std::max(spanned, 1));
    std::array<long long,3> count{};
    std::array<double,3> h{};
    for (int c = 0; c < 3; ++c) {
        const double span = hi[c] - lo[c];
        count[c] = (span > tol) ? std::clamp<long long>(std::llround(span / side), 1, 1 << 20) : 1;
        h[c] = (span > tol) ? span / count[c] : 1.0;
    }
    const auto bucketIndex = [&](int c, double x) {
        return std::clamp<long long>(static_cast<long long>(std::floor((x - lo[c]) / h[c])),
                                     0, count[c] - 1);
    };
    const auto bucketOf = [&](long long bx, long long by, long long bz) {
        return bx + count[0]*(by + count[1]*bz);
    };
    // Used corners by bucket.
    std::map<long long, std::vector<int>> buckets;
    for (int i = 0; i < numCorners; ++i) {
        if (used[i]) {
            buckets[bucketOf(bucketIndex(0, P[i][0]), bucketIndex(1, P[i][1]),
                             bucketIndex(2, P[i][2]))].push_back(i);
        }
    }

    Opm::SparseTable<int> rebuilt;
    std::vector<int> newRow;
    std::vector<std::pair<double,int>> mids;
    bool inserted = false;

    for (int face = 0; face < faceToPoint.size(); ++face) {
        auto row = faceToPoint[face];
        const int n = row.size();
        newRow.clear();
        for (int e = 0; e < n; ++e) {
            const int a = row[e];
            const int b = row[(e + 1) % n];
            newRow.push_back(a);

            Coord d = P[b];
            d -= P[a];
            const double len2 = d.two_norm2();
            if (len2 <= tol * tol) {
                continue;
            }

            // Candidate nodes: all buckets overlapping the edge's bounding box.
            std::array<long long,3> first{};
            std::array<long long,3> last{};
            for (int c = 0; c < 3; ++c) {
                first[c] = bucketIndex(c, std::min(P[a][c], P[b][c]) - tol);
                last[c] = bucketIndex(c, std::max(P[a][c], P[b][c]) + tol);
            }

            mids.clear();
            for (long long bx = first[0]; bx <= last[0]; ++bx) {
                for (long long by = first[1]; by <= last[1]; ++by) {
                    for (long long bz = first[2]; bz <= last[2]; ++bz) {
                        const auto it = buckets.find(bucketOf(bx, by, bz));
                        if (it == buckets.end()) {
                            continue;
                        }
                        for (const int c : it->second) {
                            if (c == a || c == b) {
                                continue;
                            }
                            Coord ac = P[c];
                            ac -= P[a];
                            const double t = (ac * d) / len2;
                            if (t <= 0.0 || t >= 1.0) {
                                continue; // not strictly between the endpoints
                            }
                            // Distance from the node to the line a-b.
                            Coord closest = d;
                            closest *= t;
                            closest += P[a];
                            Coord diff = P[c];
                            diff -= closest;
                            if (diff.two_norm() < tol) {
                                mids.emplace_back(t, c);
                            }
                        }
                    }
                }
            }
            std::sort(mids.begin(), mids.end());
            for (const auto& [t, c] : mids) {
                newRow.push_back(c);
            }
            inserted = inserted || !mids.empty();
        }
        rebuilt.appendRow(newRow.begin(), newRow.end());
    }

    if (inserted) {
        faceToPoint.swap(rebuilt);
    }
    for (int pass = 0; pass < 4 && insertCellChainNodes(leaf, P); ++pass) {
    }
}

} // namespace Refinement
} // namespace Opm
