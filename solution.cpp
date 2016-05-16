#ifndef LOCAL
#define NDEBUG
#endif

#include <algorithm>
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


vector<pair<int, int>> stars;
vector<bool> visited;


int dist2(int star1, int star2) {
    int dx = stars[star1].first - stars[star2].first;
    int dy = stars[star1].second - stars[star2].second;
    return dx*dx + dy*dy;
}

float dist(int star1, int star2) {
    return sqrt(dist2(star1, star2));
}


int move_number;


class StarTraveller {
public:
    int init(vector<int> stars)
    {
        for (int i = 0; i < stars.size(); i += 2) {
            ::stars.emplace_back(stars[i], stars[i + 1]);
        }
        ::visited = vector<bool>(::stars.size());
        move_number = 0;
        return 0;
    }

    vector<int> makeMoves(vector<int> ufos, vector<int> ships)
    {
        vector<int> ret = ships;

        float best_score = 10000000;
        int best_ship = 0;
        int best_dst = ships[best_ship];
        bool urgent = false;

        for (int i = 0; i < stars.size(); i++) {
            if (visited[i])
                continue;
            for (int j = 0; j < ships.size(); j++) {
                float d = dist(i, ships[j]);
                if (d < best_score) {
                    best_score = d;
                    best_ship = j;
                    best_dst = i;
                }
            }
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

            int nearest_d = 10000000;
            int nearest_not_visited = ufos[i + 2];
            for (int ii = 0; ii < stars.size(); ii++) {
                if (visited[ii] || ii == ufos[i + 1] || ii == ufos[i])
                    continue;
                int d = dist2(ii, ufos[i + 2]);
                if (d < nearest_d) {
                    nearest_d = d;
                    nearest_not_visited = ii;
                }
            }

            if (false || !visited[ufos[i + 2]]) {
                for (int j = 0; j < ships.size(); j++) {
                    float k = ships[j] == ufos[i] ? 1e-3 : 1;
                    float d =
                        k * dist(ships[j], ufos[i + 1]) +
                        1e-3 * dist(ufos[i + 1], ufos[i + 2]) +
                        dist(ufos[i + 2], nearest_not_visited);
                    if (!visited[ufos[i + 1]])
                        d *= 0.5;
                    if (d < best_score) {
                        best_score = d;
                        best_ship = j;
                        best_dst = ufos[i + 1];
                        urgent = true;
                    }
                }
            }
        }

        if (urgent || move_number % 3 == 1)
            ret[best_ship] = best_dst;

        for (int r : ret)
            visited[r] = true;
        move_number++;
        return ret;
    }
};