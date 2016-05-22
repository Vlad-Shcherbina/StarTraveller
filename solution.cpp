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
#include <ctime>

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


double get_time() {
    return 1.0 * clock() / CLOCKS_PER_SEC;
}


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
vector<vector<int>> nearest;
vector<vector<int>> next_nearest;
vector<int> first_nearest;

void init_nearest() {
    assert(nearest.empty());
    vector<int> rs(stars.size());
    iota(rs.begin(), rs.end(), 0);
    nearest = vector<vector<int>>(stars.size(), rs);
    for (int i = 0; i < stars.size(); i++) {
        auto &rs = nearest[i];
        sort(rs.begin(), rs.end(), [i](int a, int b) {
            return dist2(i, a) < dist2(i, b);
        });
        first_nearest.push_back(rs[0]);
        next_nearest.emplace_back(rs.size(), -1);
        for (int j = 1; j < rs.size(); j++)
            next_nearest.back()[rs[j - 1]] = rs[j];
    }
}

int *pnear = nullptr;

int skip_visited_nearest(int star) {
    while (true) {
        if (*pnear == -1)
            return -1;
        if (!visited[*pnear])
            return *pnear;
        *pnear = next_nearest[star][*pnear];
    }
}

int get_first_nearest(int star) {
    pnear = &first_nearest[star];
    return skip_visited_nearest(star);
}

int get_next_nearest(int star) {
    pnear = &next_nearest[star][*pnear];
    return skip_visited_nearest(star);
}


int move_number;

vector<vector<double>> ufo_range_logprob;

void update_range_logprog(vector<double> &logprob, int src, int dst) {
    int rank = find(nearest[src].begin(), nearest[src].end(), dst) -
        nearest[src].begin();
    double log_num_stars = log(stars.size());
    for (int r = 0; r < logprob.size(); r++) {
        if (r < 10) {
            logprob[r] = -1000000;
            continue;
        }
        logprob[r] +=
            log(r) + (r - 1) * log(stars.size() - rank + 0.5)
            - r * log_num_stars;
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

        double s = 0.0;
        vector<double> rank_probs;
        for (int rank = 0; rank < stars.size(); rank++) {
            rank_probs.push_back(
                range * pow(
                    (stars.size() - rank + 0.5) / stars.size(),
                    range - 1)
                / stars.size());
            s += rank_probs.back();
            if (s > 0.9999)
                break;
        }
        discrete_distribution<int> rnd_rank(
            rank_probs.begin(), rank_probs.end());

        while (path.size() < length) {
            /*int best_dist = 10000000;
            int best_star = rnd_star(rnd_gen);
            for (int j = 0; j < range; j++) {
                int star = rnd_star(rnd_gen);
                int d = dist2(star, path.back());
                if (d > 0 && d < best_dist) {
                    best_dist = d;
                    best_star = star;
                }
            }
            path.push_back(best_star);*/
            path.push_back(nearest[path.back()][rnd_rank(rnd_gen)]);
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


vector<int> greedy_path(int start, int max_len) {
    vector<int> result {start};
    while (result.size() < max_len) {
        bool found = false;
        for (int i = get_first_nearest(result.back());
                 i != -1;
                 i = get_next_nearest(result.back())) {
            if (find(result.begin(), result.end(), i) == result.end()) {
                found = true;
                result.push_back(i);
                break;
            }
        }
        if (!found)
            break;
    }
    return result;
}


void improve_path(vector<int> &path) {
    while (true) {
        bool improved = false;
        for (int i = 1; i < path.size(); i++) {
            for (int j = i + 2; j < path.size(); j++) {
                if (dist(path[i - 1], path[j - 1]) + dist(path[i], path[j]) <
                    dist(path[i - 1], path[i]) + dist(path[j - 1], path[j])) {
                    improved = true;
                    reverse(path.begin() + i, path.begin() + j);
                }
            }
        }
        if (!improved)
            break;
    }
}


class StarTraveller {
private:
    double start_time;
public:
    int init(vector<int> stars)
    {
#ifndef LOCAL
        assert(false);  // make sure assertions are disabled
#endif
        start_time = get_time();
        for (int i = 0; i < stars.size(); i += 2) {
            ::stars.emplace_back(stars[i], stars[i + 1]);
        }
        visited = vector<bool>(::stars.size());
        init_nearest();
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
            auto path = greedy_path(ships[j], 10);
            improve_path(path);
            if (path.size() > 1) {
                assert(path[0] == ships[j]);
                double d = dist(path[0], path[1]);
                if (d < best_score) {
                    best_score = d;
                    best_ship = j;
                    best_dst = path[1];
                }
            }
        }

        vector<Expectation> rides;
        for (int i = 0; i < ufos.size(); i += 3) {
            vector<int> prefix {ufos[i + 1], ufos[i + 2]};
            auto paths = simulate_ufo_paths(i / 3, prefix, 5, 10);
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
            for (int ii = get_first_nearest(ufos[i + 2]);
                 ii != -1;
                 ii = get_next_nearest(ufos[i + 2])) {
                assert(!visited[ii]);
                if (ii != ufos[i + 1]) {
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

        for (int r : ret)
            visited[r] = true;
        move_number++;

        if (find(visited.begin(), visited.end(), false) == visited.end()) {
            double it_took = get_time() - start_time;
            debug(it_took);
        }
        return ret;
    }
};
