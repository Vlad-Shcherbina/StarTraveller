#ifndef LOCAL
#define NDEBUG
#endif

#include <algorithm>
#include <numeric>
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <random>
#include <utility>
#include <set>

#include "pretty_printing.h"

using namespace std;


#define debug(x) \
    cerr << #x " = " << (x) << endl
#define debug2(x, y) \
    cerr << #x " = " << (x) \
    << ", " #y " = " << (y) << endl
#define debug3(x, y, z) \
    cerr << #x " = " << (x) \
    << ", " #y " = " << (y) \
    << ", " #z " = " << (z) << endl


default_random_engine rnd_gen(42);


vector<pair<int, int>> stars;

int dist2(int star1, int star2) {
    int dx = stars[star1].first - stars[star2].first;
    int dy = stars[star1].second - stars[star2].second;
    return dx*dx + dy*dy;
}

float dist(int star1, int star2) {
    return sqrt(dist2(star1, star2));
}


vector<bool> visited;
vector<vector<int>> nearest_not_visited;
int num_visited;

void init_visited() {
    num_visited = 0;
    visited = vector<bool>(stars.size());
    assert(nearest_not_visited.empty());
    vector<int> rs(stars.size());
    iota(rs.begin(), rs.end(), 0);
    nearest_not_visited = vector<vector<int>>(stars.size(), rs);
    for (int i = 0; i < stars.size(); i++) {
        auto &rs = nearest_not_visited[i];
        sort(rs.begin(), rs.end(), [i](int a, int b) {
            return dist2(i, a) < dist2(i, b);
        });
    }
}

void mark_visited(int star) {
    assert(!visited[star]);
    visited[star] = true;
    num_visited++;
}

const vector<int>& get_nearest_not_visited(int star) {
    vector<int> &rs = nearest_not_visited[star];
    if (rs.size() == stars.size() - num_visited)
        return rs;
    if (bernoulli_distribution(
            0.5 * (rs.size() + num_visited - stars.size()) /
            (stars.size() - num_visited + 1))(rnd_gen))
        return rs;
    auto p = remove_if(rs.begin(), rs.end(), [](int a){ return visited[a]; });
    rs.erase(p, rs.end());
    return rs;
}


int move_number;

vector<vector<double>> ufo_range_logprob;

void update_range_logprog(vector<double> &logprob, int src, int dst) {
    vector<int> dists;
    for (int i = 0; i < stars.size(); i++) {
        int d = dist2(i, src);
        if (d > 0)
            dists.push_back(d);
    }
    sort(dists.begin(), dists.end());
    auto left = lower_bound(dists.begin(), dists.end(), dist2(src, dst));
    auto right = upper_bound(dists.begin(), dists.end(), dist2(src, dst));
    if (left == right)
        return;
    for (int r = 0; r < logprob.size(); r++) {
        if (r < 10) {
            logprob[r] = -1000000;
            continue;
        }
        double a = dists.end() - left;
        double b = dists.end() - right;
        logprob[r] += log(r) + (r - 1) * log(0.5 * (a + b));
        logprob[r] -= r * log(stars.size() - 1);
    }
    double m = *max_element(logprob.begin(), logprob.end());
    for (double &lp : logprob)
        lp -= m;
}


vector<vector<int>> simulate_ufo_paths(
        int ufo, const vector<int> &prefix,
        int num_paths, int length) {
    vector<double> range_probs;
    for (auto lp : ufo_range_logprob[ufo])
        range_probs.push_back(exp(lp));
    discrete_distribution<int> rnd_range(range_probs.begin(), range_probs.end());
    uniform_int_distribution<int> rnd_star(0, stars.size() - 1);

    vector<vector<int>> result;
    for (int i = 0; i < num_paths; i++) {
        vector<int> path = prefix;
        int range = rnd_range(rnd_gen);
        while (path.size() < length) {
            int best_dist = 10000000;
            int best_star = rnd_star(rnd_gen);
            for (int j = 0; j < range; j++) {
                int star = rnd_star(rnd_gen);
                int d = dist2(star, path.back());
                if (d > 0 && d < best_dist) {
                    best_dist = d;
                    best_star = star;
                }
            }
            path.push_back(best_star);
        }
        result.push_back(path);
    }
    return result;
}


struct Expectation {
    float new_stars = 0;
    float energy = 0;
    float moves = 0;
};

Expectation ride_expectations(const vector<vector<int>> ufo_paths) {
    Expectation result;
    vector<bool> path_visited(stars.size(), false);
    for (const auto &path : ufo_paths) {
        int prev = path.front();
        for (int p : path) {
            if (!visited[p] &&
                !path_visited[p]) {
                result.new_stars += 1;
                path_visited[p] = true;
            }
            result.energy += 1e-3 * dist(prev, p);
            prev = p;
        }
        for (int p : path)
            path_visited[p] = false;
        result.moves += 1;
    }
    result.new_stars /= ufo_paths.size();
    result.energy /= ufo_paths.size();
    result.moves /= ufo_paths.size();
    return result;
}


class StarTraveller {
public:
    int init(vector<int> stars)
    {
        for (int i = 0; i < stars.size(); i += 2) {
            ::stars.emplace_back(stars[i], stars[i + 1]);
        }
        init_visited();
        move_number = 0;
        return 0;
    }

    vector<int> makeMoves(vector<int> ufos, vector<int> ships)
    {
        if (move_number == 0) {
            vector<double> t(stars.size() / 10 + 10, 0.0);
            for (int i = 0; i < 10; i++)
                t[i] = -1000000;
            ufo_range_logprob = vector<vector<double>>(ufos.size() / 3, t);

            for (int i = 0; i < ufos.size(); i += 3)
                update_range_logprog(
                    ufo_range_logprob[i / 3], ufos[i], ufos[i + 1]);
        }
        for (int i = 0; i < ufos.size(); i += 3)
            update_range_logprog(
                ufo_range_logprob[i / 3], ufos[i + 1], ufos[i + 2]);

        if (move_number == 100) {
            for (auto range : ufo_range_logprob) {
                int r = max_element(range.begin(), range.end()) - range.begin();
                debug(r);
            }
        }
        vector<int> ret = ships;

        float best_score = 10000000;
        int best_ship = 0;
        int best_dst = ships[best_ship];
        bool urgent = false;

        for (int j = 0; j < ships.size(); j++) {
            for (int i : get_nearest_not_visited(ships[j])) {
                if (visited[i])
                    continue;
                int d = dist(i, ships[j]);
                if (d < best_score) {
                    best_score = d;
                    best_ship = j;
                    best_dst = i;
                }
                break;
            }
        }

        vector<Expectation> rides;
        for (int i = 0; i < ufos.size(); i += 3) {
            vector<int> prefix {ufos[i + 1], ufos[i + 2]};
            auto paths = simulate_ufo_paths(i / 3, prefix, 2, 20);
            rides.push_back(ride_expectations(paths));
        }

        for (int i = 0; i < ufos.size(); i += 3) {
            if (!visited[ufos[i + 1]]) {
                for (int j = 0; j < ships.size(); j++) {
                    if (ships[j] != ufos[i])
                        continue;
                    float d = 1e-3 * dist(ufos[i + 1], ships[j]);
                    if (d < best_score) {
                        best_score = d;
                        best_ship = j;
                        best_dst = ufos[i + 1];
                        urgent = true;
                    }
                }
            }

            int nnv = ufos[i + 2];
            for (int ii : get_nearest_not_visited(ufos[i + 2])) {
                if (!visited[ii] && ii != ufos[i + 1]) {
                    nnv = ii;
                    break;
                }
            }

            for (int j = 0; j < ships.size(); j++) {
                float k = ships[j] == ufos[i] ? 1e-3 : 1;
                float d =
                    k * dist(ships[j], ufos[i + 1]) +
                    1e-3 * dist(ufos[i + 1], ufos[i + 2]) +
                    dist(ufos[i + 2], nnv);
                if (!visited[ufos[i + 1]])
                    d *= 0.5;
                if (d < best_score) {
                    best_score = d;
                    best_ship = j;
                    best_dst = ufos[i + 1];
                    urgent = true;
                }
            }

            for (int j = 0; j < ships.size(); j++) {
                float d = dist2(ships[j], ufos[i + 1]);
                if (ships[j] == ufos[i])
                    d *= 1e-3;
                if (rides[i / 3].moves > 4 * rides[i / 3].new_stars)
                    continue;
                d += rides[i / 3].energy;
                d /= 1e-6 + rides[i / 3].new_stars;
                if (d < best_score) {
                    best_score = d;
                    best_ship = j;
                    best_dst = ufos[i + 1];
                    urgent = true;
                }
            }
        }

        if (urgent || move_number % 3 == 1)
            ret[best_ship] = best_dst;

        for (int r : ret) {
            if (!visited[r])
                mark_visited(r);
        }
        move_number++;
        return ret;
    }
};
