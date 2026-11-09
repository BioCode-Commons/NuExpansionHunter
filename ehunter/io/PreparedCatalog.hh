#pragma once
#include "locus/LocusSpecification.hh"

namespace ehunter
{
void prepareCatalogCache(
    const std::string& path, const RegionCatalog& catalog, const InputPaths& inputs,
    const HeuristicParameters& parameters, const ReferenceContigInfo& contigs);
RegionCatalog loadPreparedCatalog(
    const std::string& path, const InputPaths& inputs, const HeuristicParameters& parameters,
    const ReferenceContigInfo& contigs);
}
