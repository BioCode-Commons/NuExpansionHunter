// NuExpansionHunter command-line resource controls.

#include "io/ParameterLoading.hh"

#include <algorithm>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include "gtest/gtest.h"

namespace ehunter
{
namespace
{

class ThreadParametersTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        directory_ = boost::filesystem::temp_directory_path()
            / boost::filesystem::unique_path("nuexpansionhunter-parameters-%%%%-%%%%-%%%%");
        boost::filesystem::create_directory(directory_);
        for (const auto* file : { "reads.bam", "reads.bam.bai", "reference.fa", "catalog.json" })
        {
            std::ofstream output((directory_ / file).string());
            ASSERT_TRUE(output.good());
        }
    }

    void TearDown() override
    {
        boost::system::error_code error;
        boost::filesystem::remove_all(directory_, error);
    }

    ProgramParameters load(const std::vector<std::string>& options, const std::string& mode = "streaming")
    {
        // Loading parameters checks file existence, but does not parse file contents.
        std::vector<std::string> arguments = { "ExpansionHunter", "--reads", (directory_ / "reads.bam").string(),
            "--reference", (directory_ / "reference.fa").string(), "--variant-catalog",
            (directory_ / "catalog.json").string(), "--output-prefix", (directory_ / "out").string(),
            "--analysis-mode", mode };
        arguments.insert(arguments.end(), options.begin(), options.end());
        std::vector<char*> argv;
        for (auto& argument : arguments)
            argv.push_back(&argument[0]);
        auto parameters = tryLoadingProgramParameters(static_cast<int>(argv.size()), argv.data());
        if (!parameters)
            throw std::logic_error("Expected parsed parameters");
        return *parameters;
    }

private:
    boost::filesystem::path directory_;
};

TEST_F(ThreadParametersTest, IoDefaultsStaySynchronousAsComputeWorkersIncrease)
{
    const auto defaults = load({});
    EXPECT_EQ(defaults.threadCount, 1);
    EXPECT_EQ(defaults.decodeThreadCount, 0);
    EXPECT_EQ(defaults.encodeThreadCount, 0);

    for (const int threads : { 1, 4, 32 })
    {
        const auto parameters = load({ "--threads=" + std::to_string(threads) });
        EXPECT_EQ(parameters.threadCount, threads);
        EXPECT_EQ(parameters.decodeThreadCount, 0);
        EXPECT_EQ(parameters.encodeThreadCount, 0);

        const ProgramParameters constructed(parameters.inputPaths(), parameters.outputPaths(), parameters.sample(),
            parameters.heuristics(), parameters.analysisMode(), parameters.logLevel(), threads, false);
        EXPECT_EQ(constructed.decodeThreadCount, 0);
        EXPECT_EQ(constructed.encodeThreadCount, 0);
    }
}

TEST_F(ThreadParametersTest, ExplicitAutoRetainsLegacyIoPoolSizing)
{
    for (const int threads : { 1, 4, 32 })
    {
        const auto parameters = load(
            { "--threads=" + std::to_string(threads), "--decode-threads=auto", "--encode-threads=auto" });
        EXPECT_EQ(parameters.threadCount, threads);
        EXPECT_EQ(parameters.decodeThreadCount, std::min(threads, 12));
        EXPECT_EQ(parameters.encodeThreadCount, threads);
    }
}

TEST_F(ThreadParametersTest, IoPoolsCanBeConfiguredIndependently)
{
    const auto parameters = load({ "--threads=4", "--decode-threads=2", "--encode-threads=3" });
    EXPECT_EQ(parameters.threadCount, 4);
    EXPECT_EQ(parameters.decodeThreadCount, 2);
    EXPECT_EQ(parameters.encodeThreadCount, 3);
    EXPECT_EQ(load({ "--decode-threads=2" }).encodeThreadCount, 0);
    EXPECT_EQ(load({ "--encode-threads=3" }).decodeThreadCount, 0);

    for (const auto& value : { "0", "1" })
    {
        const auto small = load({ std::string("--decode-threads=") + value,
            std::string("--encode-threads=") + value });
        EXPECT_EQ(small.decodeThreadCount, value[0] - '0');
        EXPECT_EQ(small.encodeThreadCount, value[0] - '0');
    }
}

TEST_F(ThreadParametersTest, InvalidAndOverflowingIoCountsAreRejected)
{
    const std::vector<std::string> invalid = { "-1", "1.5", "1e2", "AUTO", "2workers",
        std::to_string(static_cast<unsigned long long>(std::numeric_limits<int>::max()) + 1),
        "99999999999999999999999999999999999999999" };
    for (const auto& option : { "--decode-threads=", "--encode-threads=" })
    {
        for (const auto& value : invalid)
        {
            SCOPED_TRACE(std::string(option) + value);
            EXPECT_THROW(load({ std::string(option) + value }), std::invalid_argument);
        }
    }
    EXPECT_THROW(load({ "--threads=0" }), std::invalid_argument);
}

TEST_F(ThreadParametersTest, SeekingModeDoesNotSilentlyIgnoreExplicitDecodeThreads)
{
    const auto defaults = load({ "--threads=32" }, "seeking");
    EXPECT_EQ(defaults.decodeThreadCount, 0);
    EXPECT_EQ(defaults.encodeThreadCount, 0);
    EXPECT_EQ(load({ "--threads=4", "--encode-threads=3" }, "seeking").encodeThreadCount, 3);
    EXPECT_EQ(load({ "--threads=32", "--decode-threads=auto" }, "seeking").decodeThreadCount, 12);
    EXPECT_THROW(load({ "--decode-threads=0" }, "seeking"), std::invalid_argument);
    EXPECT_THROW(load({ "--decode-threads=1" }, "seeking"), std::invalid_argument);
    EXPECT_THROW(load({ "--decode-threads=2" }, "seeking"), std::invalid_argument);
}

TEST_F(ThreadParametersTest, QueueBudgetSupportsDefaultExplicitAndUnboundedValues)
{
    const std::size_t bytesPerMiB = 1024 * 1024;
    EXPECT_EQ(load({}).bamletQueueBytes, 0u);
    EXPECT_EQ(load({ "--bamlet-queue-mb=3" }).bamletQueueBytes, 3 * bytesPerMiB);
    EXPECT_EQ(load({ "--bamlet-queue-mb=0" }).bamletQueueBytes, 0u);

    const auto maxMiB = std::numeric_limits<std::size_t>::max() / bytesPerMiB;
    EXPECT_EQ(load({ "--bamlet-queue-mb=" + std::to_string(maxMiB) }).bamletQueueBytes, maxMiB * bytesPerMiB);
    EXPECT_THROW(load({ "--bamlet-queue-mb=" + std::to_string(maxMiB + 1) }), std::invalid_argument);
}

TEST_F(ThreadParametersTest, InvalidQueueBudgetsAreRejected)
{
    for (const auto* value : { "-1", "1.5", "auto", "1MB", "99999999999999999999999999999999999999999" })
    {
        SCOPED_TRACE(value);
        EXPECT_THROW(load({ std::string("--bamlet-queue-mb=") + value }), std::invalid_argument);
    }
}

}
}
