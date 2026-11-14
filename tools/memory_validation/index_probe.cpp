#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <unistd.h>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <vector>
#include "graphalign/KmerIndex.hh"
#include "graphcore/GraphBuilders.hh"
#include "graphutils/KmerEncoding.hh"

using namespace graphtools;
std::string randomSequence(size_t length, uint32_t seed) {
    std::string s(length, 'A');
    for (char& c : s) { seed ^= seed << 13; seed ^= seed >> 17; seed ^= seed << 5; c = "ACGT"[seed&3]; }
    return s;
}
uint64_t digest(uint64_t h, const std::string& s) {
    for (unsigned char c : s) { h ^= c; h *= 1099511628211ULL; } return h;
}
int main(int argc, char** argv) {
    const bool audit = argc == 1 || std::string(argv[1]) == "audit";
    int k = argc > 2 ? std::stoi(argv[2]) : 14;
    const std::string kind = argc > 3 ? argv[3] : "hit";
    const int rounds = argc > 4 ? std::stoi(argv[4]) : 100;
    const bool membership = k == 10;
    const std::string order = argc > 5 ? argv[5] : "shuffled";
    std::vector<Graph> graphs; graphs.reserve(32);
    for (int i=0; i<32; ++i) {
        auto left = randomSequence(i==30 ? 5000 : i==31 ? 100 : 1000, 1773+i*31);
        auto right = randomSequence(1000, 443+i*13);
        if (i%5==0) left[83]='N';
        if (i%7==0) right.replace(300, 30, left.substr(700,30));
        if (i==29) for (size_t j=0;j+16<left.size();j+=16) left.replace(j,4,"AAAA");
        graphs.push_back(i%4 ? makeStrGraph(left, i%2?"CAG":"T", right) : makeDeletionGraph(left,"GGG",right));
    }
    std::vector<std::unique_ptr<KmerIndex>> indexes;
    std::vector<std::vector<uint32_t>> queries;
    uint64_t h = 14695981039346656037ULL, count = 0, recordBytes=0;
    for (size_t g=0; g<graphs.size(); ++g) {
        indexes.emplace_back(new KmerIndex(graphs[g],k,membership));
        auto& index = *indexes.back();
        std::ostringstream out; index.writeRecord(out); recordBytes += out.str().size();
        if (audit) {
            const auto record=out.str(); uint64_t identity; std::memcpy(&identity, record.data()+24,8);
            uint32_t bucketCount; std::memcpy(&bucketCount,record.data()+24+16+identity+4,4);
            std::cerr << g << " " << k << " " << record.size() << " " << bucketCount << "\n";
            std::cout << "encode " << g << " " << index.encode() << "\n";
        }
        std::unique_ptr<KmerIndex> retained;
        if (audit) {
            const std::string path = "/tmp/eh-index-memory-" + std::to_string(getpid()) + "-" + std::to_string(g);
            KmerIndex::prepareArchive(path,1,[&](std::ostream& output,uint64_t){index.writeRecord(output);});
            auto archive=std::make_shared<KmerIndexArchive>(path);
            KmerIndex mapped(graphs[g],k,membership,archive,0), copied(mapped);
            retained.reset(new KmerIndex(std::move(copied)));
            archive.reset(); unlink(path.c_str());
            assert(index.encode()==retained->encode());
        }
        auto set = index.kmers(); std::vector<std::string> keys(set.begin(),set.end()); std::sort(keys.begin(),keys.end());
        queries.emplace_back(); TwoBitKmerEncoder encoder(k);
        for (const auto& key : keys) {
            auto word = encoder.encode(key);
            if (kind=="miss") word = (word*1812347+3172367) & ((uint32_t(1)<<(2*k))-1);
            queries.back().push_back(word);
            if (audit) {
                std::cout << "key " << g << " " << key;
                assert(retained->containsEncodedKmer(word));
                if (!membership) {
                    assert(index.getPaths(key)==retained->getPaths(key));
                    for (const auto& p : index.getPaths(key)) std::cout << " " << p.encode();
                }
                std::cout << "\n";
            }
        }
        if (!audit && order == "rolling") {
            queries.back().clear();
            for (NodeId node=0;node<graphs[g].numNodes();++node) {
                const std::string sequence=graphs[g].nodeSeq(node);
                for(size_t position=0;position+k<=sequence.size();++position) {
                    const auto key=sequence.substr(position,k);
                    if(key.find_first_not_of("ACGT")!=std::string::npos) continue;
                    queries.back().push_back(encoder.encode(key));
                }
            }
        } else if (!audit && order == "shuffled") {
            std::mt19937 random(137+g);
            std::shuffle(queries.back().begin(),queries.back().end(),random);
        }
    }
    if (audit) return 0;
    const auto begin = std::chrono::steady_clock::now();
    for (int r=0;r<rounds;++r) for(size_t g=0;g<indexes.size();++g) for (auto key : queries[g]) {
        if(kind=="paths") {
            auto paths=indexes[g]->findEncodedPaths(key); count+=paths.size();
            if(paths.size()) { const auto p=paths[0]; h += p.startPosition()+p.endPosition(); }
        } else { count += indexes[g]->containsEncodedKmer(key); h += key; }
    }
    const auto end = std::chrono::steady_clock::now();
    std::cout << count << " " << h << "\n";
    std::cerr << "seconds " << std::chrono::duration<double>(end-begin).count() << " bytes " << recordBytes << "\n";
}
