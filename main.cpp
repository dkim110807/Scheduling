#include <array>
#include <algorithm>
#include <iostream>
#include <vector>

#ifdef LOCAL

#include "include/debug.h"

#else
#define debug(...) 0x110807
#endif

class Block {
public:
    const int64_t p;

    std::vector<std::array<int64_t, 3>> v;

    explicit Block(int64_t p = 0) : p(p) {

    }

    [[nodiscard]]
    int size() const {
        return (int) v.size();
    }

    [[nodiscard]]
    bool empty() const {
        return v.empty();
    }

    void push_back(const std::array<int64_t, 3> &job) {
        v.push_back(job);
    }

    std::array<int64_t, 3> &operator[](int i) {
        if (i >= size() || i <= 0) {
            throw std::out_of_range("Index range must be [0, " + std::to_string(size() - 1)
                                    + "] but found " + std::to_string(i) + "!");
        }
        return v[i];
    }

    /*
     * Release time of the block
     * - r(B) = min_{j \in B} r(j)
     */
    friend int64_t r(const Block &b) {
        if (b.empty()) {
            throw std::out_of_range("Can't find the release time of an empty block!");
        }
        return b.v.front()[0];
    }

    /*
     * Finish time of the block, since all the processing time are equal
     * - t(B) = r(B) + p|B|
     */
    friend int64_t t(const Block &b) {
        return r(b) + b.p * (int64_t) b.v.size();
    }

#ifdef LOCAL

    friend std::string to_string(const Block &b) {

        return "";
    }

#endif
};

/*
 * An O(n^2) algorithm for scheduling equal-length preemptive jobs on a single machine to minimize total tardiness
 * - Zhongjun Tian, C.T.Ng, T.C.E. Cheng
 *
 * See README.md for more information
 *
 * Implementation by dkim110807
 */
int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr), std::cout.tie(nullptr);

    int64_t n, p;   // n : number of jobs in schedule, p : the time needed for each job
    std::cin >> n >> p;

    std::vector<std::array<int64_t, 3>> jobs(n);
    for (int i = 0; i < n; i++) {
        std::cin >> jobs[i][0] >> jobs[i][1];
        jobs[i][1] += jobs[i][0];   // BOJ 26248
        jobs[i][2] = i;
    }

    debug(jobs);

    /*
     * Algorithm SMPP-ELJTT (SMPP, Equal-length jobs, total tardiness)
     *
     * 1. Set H_{1} = H_{2} = \cdots = H_{n} = \left\{\right\}, jobs = 1
     * 2. Index all the jobs by the ERD rule
     * 3. Schedule all the jobs by the ERD rule, and all the blocks to H_{1}
     * 4. while H_{j} \neq \left\{\right\}
     *        Index all the blocks/subblocks in H_{j} as B_{1}, B_{2}, \cdots, B_{|H_{j}|} by their start times
     *        for i = 1, 2, ..., |H_{j}|:
     *            if |B_{i}| = 1:
     *                Schedule the unique job in B_{i} within I_{i}
     *            else:
     *                Schedule B_{i} by Algorithm BLK-DE, and add all optimal subblocks in the final result to H_{j + 1}
     *            end if
     *        end for
     *        j = j + 1
     */

    // Algorithm SMPP-ELJTT
    std::vector<std::vector<Block>> H(n); // I will use 0-based index

    // Sort all the jobs by the ERD (Earliest Release Data) rule
    std::sort(jobs.begin(), jobs.end(), [&](const std::array<int64_t, 3> &a, const std::array<int64_t, 3> &b) -> bool {
        return a[0] < b[0];     // a[0] stands for the release date
    });

    Block b = Block(p);
    b.push_back(jobs[0]);
    H[0].push_back(b);
    debug(H[0]);

    debug(r(H[0][0]));

    for (int j = 0; j < n; j++) {
        if (H[j].empty()) break;
        for (int i = 0; i < H[j].size(); i++) {
            auto &B = H[j][i];  // The block (subblock) that we are currently working on
            if (B.size() == 1) {    // Case 1. Schedule the unique job in the interval

                continue;
            }

        }
    }

    /*
     * Algorithm BLK-DE (Block decomposition)
     *
     * 1. Set U = \left\{\right\}, u = 1, V = B_{i} and L = t(B_{i})
     * 2. while u \leq n:
     *        if |V| = 1:
     *            Move the unique job in V to U as j_{u}
     *            Reset t(V) to the ending time of the second last subblock
     *            break
     *        else:
     *            Select a job from V by the LDD rule and move it to U as j_{u}
     *            Schedule V \ {j_{u}} by the ERD rule
     *            Reset V to be the new last subblock and declare all the new subblocks except V to be optimal
     *
     *            if d_{j_{u}} \geq t(V):
     *                Declare V to be optimal
     *                break
     *            else if \Delta_{u} = 0:
     *                break
     *            else:
     *                u = u + 1
     *            end if
     *        end if
     * 3. Set k = u
     *    while L - t(V) > p:
     *
     */
}
