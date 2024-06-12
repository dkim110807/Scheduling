#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <iostream>
#include <vector>

#ifdef LOCAL

#include "include/debug.h"

#else
#define debug(...) 0x110807
#endif

class Block {
public:
    int64_t p;

    // containing jobs (release, deadline, idx)
    std::vector<std::array<int64_t, 3>> v;

    explicit Block(int64_t p = 0) : p(p) {

    }

    ~Block() = default;

    Block operator=(const Block &b) {
        p = b.p;

        return b;
    }

    [[nodiscard]]
    size_t size() const {
        return v.size();
    }

    [[nodiscard]]
    bool empty() const {
        return v.empty();
    }

    void push_back(const std::array<int64_t, 3> &job) {
        v.push_back(job);
    }

    std::array<int64_t, 3> &operator[](size_t i) {
        if (i >= size() || i < 0) {
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
        std::string s = to_string(b.v);
        return s;
    }

#endif
};

int64_t r(const std::array<int64_t, 3> &j) {
    return j[0];
}

int64_t d(const std::array<int64_t, 3> &j) {
    return j[1];
}

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

    debug(jobs);

    // Make the initial blocks
    for (size_t i = 0, j = 0; i < n; i = j) {
        Block B(p);
        B.push_back(jobs[j]);
        for (j++; j < n; j++) {
            if (t(B) > jobs[j][0]) {
                B.push_back(jobs[j]);
            } else break;
        }
        H[0].push_back(B);
    }

    debug(H[0]);

    std::vector<std::array<int64_t, 3>> schedule;

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
     *    end
     * 3. Set k = u
     *    while L - t(V) > p:
     *        Schedule job j_{k} within [L - p, L]
     *        Remove j_k from U
     *        Reset L = L - p
     *        while U \neq = {}
     *            Re-index the jobs in U as j_1, ..., j_{U} by the ERD rule
     *            Set u = |U| and \Lambda = 0
     *            for k = 1, ..., u:
     *                \Lambda = \Lambda + \lambda_{k}
     *                if d_{j_{k}} < t(V):
     *                    break from all loops and goto Step 2
     *                else if d_{j_{k}} \geq L - \Lambda:
     *                    break
     *                end
     *            if min{\Lambda, L - t(V)} > p:
     *                break
     *            else if \Lambda \geq L - t(V):
     *                break from all loops and goto Step 4.
     *            else:
     *                Schedule {j_{1}, ..., j_{k}} within \Delta_{j = 1 ... k} and [L - \Lambda, L] by the ERD rule
     *                Remove {j1, ..., j_{k}} from U
     *            end
     *        end
     *    end
     *  4. if V is not declared to be optimal:
     *         Declare it to be optimal
     *     else:
     *         Schedule all jobs in U\{j_{k}} within \Delta_{j = 1 ... u} by the ERD rule
     *         Schedule job j_{k} within the rest of \Delta_{u} and [t(V), L]
     *         return
     *     end
     */

    // This is used in Algorithm BLK-DE - step 2, 3
    std::vector<bool> optimal(n, false);

    for (size_t j = 0; j < n; j++) {
        if (H[j].empty()) break;
        for (size_t i = 0; i < H[j].size(); i++) {
            Block &B = H[j][i];  // The block (subblock) that we are currently working on
            if (B.size() == 1) {    // Case 1. Schedule the unique job in the interval
                schedule.push_back({B[0][0], B[0][0] + p, B[0][2]});
                continue;
            }

            // Block decomposition
            size_t u = 0;   // 0-based index
            Block V = B;
            int64_t L = t(V);
            std::vector<std::array<int64_t, 3>> U;
            std::vector<int64_t> Delta; // \Delta
            std::vector<std::array<int64_t, 2>> Lambda; // \Lambda
            std::vector<std::vector<std::array<int64_t, 2>>> delta; // Representing \delta_{i}

            // Step 2.
            std::function<void()> step2 = [&]() -> void {
                while (u < n) {
                    if (V.size() == 1) {
                        U.push_back((V)[0]);
                        assert(!H[j + 1].empty());
                        V = H[j + 1].back();
                        H[j + 1].pop_back();
                    } else {
                        std::array<int64_t, 3> ju{-1, -1, -1};  // j_{u}, LDD rule
                        // find the job j_{u} by the LDD rule
                        for (size_t i = 0; i < V.size(); i++) {
                            if (V[i][1] > ju[1]) ju = V[i];
                        }

                        debug(ju);

                        size_t start = H[j + 1].size();



                        // V is sorted in ERD rule therefore no need for extra sorting
                        for (size_t i = 0, k = 0; i < n; i = k) {
                            Block B(p);
                            if (k == ju[2]) k++;
                            B.push_back(jobs[k]);
                            for (k++; k < n; k++) {
                                if (k == ju[2]) k++;
                                if (t(B) > jobs[k][0]) {
                                    B.push_back(jobs[k]);
                                } else break;
                            }
                            H[j + 1].push_back(B);
                        }

                        V = H[j + 1].back();
                        H[j + 1].pop_back();

                        for (size_t i = start; i < H[j + 1].size(); i++) {
                            for (size_t k = 0; k < H[j + 1][start].size(); k++) {
                                optimal[H[j + 1][start][k][2]] = true;
                            }
                        }

                        if (d(ju) >= t(V)) {
                            // Declare V to be optimal
                            for (size_t i = 0; i < V.size(); i++) {
                                optimal[V[i][2]] = true;
                            }
                            H[j + 1].push_back(V);
                            break;
                        } else if (false) { // Todo.
                            break;
                        } else {
                            u = u + 1;
                        }
                    }
                }
            };

            // Step 3.
            size_t k = -1;
            std::function<void()> step3 = [&]() -> void {
                k = u;
                while (L - t(V) > p) {
                    /* Todo.
                     * Schedule job j_{k} within [L - p, L]
                     * Remove j_{k} from U
                     * Reset L = L - p
                     */
                    while (!U.empty()) {

                    }
                }
            };

            // Step 4.
            std::function<void()> step4 = [&]() -> void {
                if (!optimal[V[0][2]]) {
                    for (size_t i = 0; i < V.size(); i++) {
                        optimal[V[i][2]] = true;
                    }
                    H[j + 1].push_back(V);
                } else {
                    // Schedule U \ {j_{k}} within \Delta_{1, ..., u} by the ERD rule
                    for (int i = 0; i < U.size(); i++) {
                        if (i == k) continue;

                    }
                    // Schedule j_{k} within the rest \Delta and [t(V), L]
                }
            };

            step2();
            step3();
            step4();
        }
        H[j].clear();
    }

    std::sort(schedule.begin(), schedule.end(),
              [&](const std::array<int64_t, 3> &a, const std::array<int64_t, 3> &b) -> bool {
                  return a[0] < b[0];
              }
    );

    debug(schedule);

    assert(schedule.size() <= 10 * n);

    // The output is 1 based idx
    for (const auto &[s, e, idx]: schedule) {
        std::cout << s << " " << e << " " << idx + 1 << "\n";
    }
}
