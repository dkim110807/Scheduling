//
// Created by dkim110807 on 2024-06-11.
//

#include <array>
#include <iostream>
#include <vector>

#ifdef LOCAL
#include "include/debug.h"
#else
#define debug(...) 0x110807
#endif

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr), std::cout.tie(nullptr);

    int n, p;   // n : number of jobs in schedule, p : the time needed for each job
    std::cin >> n >> p;

    std::vector<std::array<int, 3>> jobs(n);
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
     *            else
     *                Schedule B_{i} by Algorithm BLK-DE, and add all optimal subblocks in the final result to H_{j + 1}
     *            end if
     *        end for
     *        j = j + 1
     */
}
