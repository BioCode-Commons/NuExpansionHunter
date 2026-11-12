// Exercises the actual production insertion routine without initializing any
// query profiles. Private names occur only in explicit template instantiations,
// as permitted by [temp.explicit]; no layout casts or preprocessor rewrites.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include "graphalign/dagAligner/AffineAlignMatrixVectorized.hh"
#include "graphalign/dagAligner/BaseMatchingPenaltyMatrix.hh"

namespace da = graphalign::dagAligner;
using Score = da::Score;
using Matrix = da::AffineAlignMatrixVectorized<da::BaseMatchingPenaltyMatrix, true>;
template <typename Tag, typename Tag::Member member> struct Access
{
    friend typename Tag::Member access(Tag) { return member; }
};
struct VTag { using Member = da::PaddedAlignMatrix<16> Matrix::*; friend Member access(VTag); };
struct LeftTag { using Member = std::vector<Score> Matrix::*; friend Member access(LeftTag); };
struct TraceTag { using Member = std::vector<uint64_t> Matrix::*; friend Member access(TraceTag); };
struct WordsTag { using Member = std::size_t Matrix::*; friend Member access(WordsTag); };
struct RecomputeTag { using Member = void (Matrix::*)(int, int); friend Member access(RecomputeTag); };
template struct Access<VTag, &Matrix::v_>;
template struct Access<LeftTag, &Matrix::insertionLeft_>;
template struct Access<TraceTag, &Matrix::insertionTrace_>;
template struct Access<WordsTag, &Matrix::insertionWordsPerRow_>;
template struct Access<RecomputeTag, &Matrix::recomputeForInsertion>;

struct Runner
{
    Matrix matrix;
    const Score gapOpen;
    const Score gapExt;
    uint64_t cases = 0;
    uint64_t cells = 0;

    Runner(Score opening, Score extension)
        : matrix(da::BaseMatchingPenaltyMatrix(5, -4), opening, extension)
        , gapOpen(opening), gapExt(extension) {}

    void check(const std::vector<Score>& base, Score initialF, Score initialV)
    {
        const int qLen = base.size();
        auto& v = matrix.*access(VTag{});
        auto& left = matrix.*access(LeftTag{});
        auto& trace = matrix.*access(TraceTag{});
        auto& words = matrix.*access(WordsTag{});
        v.reset(qLen, 1);
        v.at(-1, 0) = initialV;
        for (int i = 0; i < qLen; ++i) v.at(i, 0) = base[i];
        left.assign(2, initialF);
        words = (qLen + 63) / 64;
        const uint64_t poison = UINT64_C(0x5a98765ffedc1234);
        trace.assign(2 * words, poison);

        // Deliberately independent sequential recurrence, including conversion
        // before max and the exceptional promoted opening-score comparison.
        std::vector<Score> expected(qLen);
        std::vector<uint64_t> expectedTrace(words, 0);
        Score previousF = initialF, previousV = initialV;
        for (int i = 0; i < qLen; ++i)
        {
            const Score extension = previousF + gapExt;
            const Score opening = previousV + gapOpen + gapExt;
            const Score currentF = std::max(extension, opening);
            const Score currentV = std::max(base[i], currentF);
            const Score extensionDifference = currentV - previousF;
            const Score openingDifference = currentV - previousV;
            const bool insertion = extensionDifference == gapExt
                || openingDifference == int(gapOpen) + int(gapExt);
            expected[i] = currentV;
            expectedTrace[i / 64] |= uint64_t(insertion) << (i % 64);
            previousF = currentF;
            previousV = currentV;
        }

        (matrix.*access(RecomputeTag{}))(qLen, 0);
        for (int i = 0; i < qLen; ++i)
            if (v.at(i, 0) != expected[i])
                throw std::runtime_error("V mismatch at case " + std::to_string(cases) + " cell " + std::to_string(i));
        for (std::size_t w = 0; w < words; ++w)
            if (trace[w] != poison || trace[words + w] != expectedTrace[w])
                throw std::runtime_error("trace mismatch at case " + std::to_string(cases));
        if (left[1] != initialF || v.at(-1, 0) != initialV)
            throw std::runtime_error("input boundary mutated");

        // Canonical output includes the entire physical V table (both rows,
        // boundaries and padding) and all allocated trace words.
        const uint32_t count = qLen;
        std::fwrite(&count, sizeof(count), 1, stdout);
        std::fwrite(&gapOpen, sizeof(gapOpen), 1, stdout);
        std::fwrite(&gapExt, sizeof(gapExt), 1, stdout);
        const std::size_t physicalSize = v.paddedRowLen() * 2;
        std::fwrite(&*v.cellZeroZero(), sizeof(Score), physicalSize, stdout);
        std::fwrite(trace.data(), sizeof(uint64_t), trace.size(), stdout);
        ++cases;
        cells += qLen;
    }
};

int main()
{
    try
    {
        Runner regular(-8, -2);
        // Every int16 value in every incoming state and every lane. Two safe
        // contexts exercise the guard boundary and the ordinary positive row.
        for (int context : { int(da::SCORE_MIN) + 48, 200 })
            for (int field = 0; field < 18; ++field)
                for (int score = -32768; score <= 32767; ++score)
                {
                    std::vector<Score> base(16, Score(context));
                    Score initialF = context, initialV = context;
                    if (field == 0) initialF = score;
                    else if (field == 1) initialV = score;
                    else base[field - 2] = score;
                    regular.check(base, initialF, initialV);
                }

        std::vector<Score> borders;
        for (int delta = 0; delta <= 56; ++delta) borders.push_back(Score(da::SCORE_MIN + delta));
        for (int score : { -1024, -24, -10, -2, -1, 0, 1, 2, 10, 24, 1024, 32758, 32766, 32767 })
            borders.push_back(Score(score));
        // Cartesian incoming boundary states and equal base rows produce both
        // guard outcomes and many tied F/V alignments.
        for (Score initialF : borders)
            for (Score initialV : borders)
                for (Score base : borders)
                    regular.check(std::vector<Score>(16, base), initialF, initialV);

        std::mt19937 random(0x9475ab13);
        uint64_t cases = regular.cases, cells = regular.cells;
        // Cross word boundaries; place low values at changing positions so
        // safe-prefix/fallback/safe-next-word transitions are exercised.
        for (const auto& gaps : std::vector<std::pair<Score, Score>>{
            {-8,-2}, {-4,-1}, {0,0}, {-32768,-2}, {32767,2}, {-8,2}, {8,-2}, {-32768,-32768}})
        {
            Runner runner(gaps.first, gaps.second);
            for (int length : {1, 3, 4, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129})
                for (int iteration = 0; iteration < 2500; ++iteration)
                {
                    std::vector<Score> base(length);
                    const int mode = iteration % 5;
                    auto nextScore = [&]() -> Score {
                        if (mode == 0) return Score(random());
                        if (mode == 1) return borders[random() % borders.size()];
                        if (mode == 2) return Score(int(da::SCORE_MIN) + 24 + int(random() % 48));
                        if (mode == 3) return Score(32767 - int(random() % 48));
                        return Score(int(random() % 2001) - 1000);
                    };
                    for (Score& score : base) score = nextScore();
                    if (iteration % 7 == 0) base[random() % length] = da::SCORE_MIN;
                    const Score initialF = nextScore();
                    const Score initialV = nextScore();
                    runner.check(base, initialF, initialV);
                }
            cases += runner.cases;
            cells += runner.cells;
        }
        if (std::ferror(stdout)) throw std::runtime_error("canonical output write failed");
        std::cerr << "{\"cases\":" << cases << ",\"cells\":" << cells
                  << ",\"oracle_matches\":true";
#if defined(__x86_64__) && defined(GRAPHTOOLS_DP_TEST_INSTRUMENT_AVX2)
        std::cerr << ",\"avx2_blocks\":" << Matrix::avx2InsertionTestBlockCounter();
#endif
        std::cerr << "}\n";
    }
    catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
