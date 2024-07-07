#ifndef UTILS_H
#define UTILS_H

inline std::vector<int> generate_blocks(int nobs, int nblocks) {
    std::vector<int> blocks(nobs);
    for (int i = 0; i < nobs; ++i) {
        blocks[i] = i % nblocks;
    }
    return blocks;
}

inline std::vector<std::shared_ptr<tatami::NumericMatrix> > fragment_matrices_by_block(const std::shared_ptr<tatami::NumericMatrix>& x, const std::vector<int>& block, int nblocks) {
    std::vector<std::shared_ptr<tatami::NumericMatrix> > collected;

    for (int b = 0; b < nblocks; ++b) {
        std::vector<int> keep;

        for (size_t i = 0; i < block.size(); ++i) {
            if (block[i] == b) {
                keep.push_back(i);
            }
        }

        if (keep.size() > 1) {
            collected.push_back(tatami::make_DelayedSubset<1>(x, keep));
        }
    }

    return collected;
}

#endif
