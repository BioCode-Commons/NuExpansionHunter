//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//
//
// Modified by NuExpansionHunter contributors, 2026.
// This file differs from ExpansionHunter v5.0.0: VCF field descriptions are collected in parallel blocks and merged in order.
// Modifications are licensed under the Apache License, Version 2.0.
//

#pragma once

#include <iostream>
#include <map>
#include <string>

#include "locus/LocusFindings.hh"
#include "locus/LocusSpecification.hh"

namespace ctpl { class thread_pool; }

namespace ehunter
{

enum class FieldType
{
    kInfo = 0,
    kFilter = 1,
    kFormat = 2,
    kAlt = 3
};

struct FieldDescription
{
    FieldDescription(
        FieldType fieldType, std::string id, std::string number, std::string contentType, std::string description);
    FieldType fieldType;
    std::string id;
    std::string number;
    std::string contentType;
    std::string description;
};

using FieldDescriptionIdentifier = std::pair<FieldType, std::string>;
using FieldDescriptionCatalog = std::map<FieldDescriptionIdentifier, FieldDescription>;

// Generates VCF field descriptions required for a given variant call
class FieldDescriptionWriter : public VariantFindingsVisitor
{
public:
    FieldDescriptionWriter(const LocusSpecification& locusSpec, const VariantSpecification& variantSpec)
        : locusSpec_(locusSpec)
        , variantSpec_(variantSpec)
    {
    }

    FieldDescriptionWriter(const LocusSpecification& locusSpec, const VariantSpecification& variantSpec,
                           FieldDescriptionCatalog& block)
        : locusSpec_(locusSpec), variantSpec_(variantSpec), block_(&block) {}

    ~FieldDescriptionWriter() = default;

    void visit(const RepeatFindings* repeatFindingsPtr) override;
    void visit(const SmallVariantFindings* smallVariantFindingsPtr) override;

    void tryAddingFieldDescription(
        FieldType fieldType, const std::string& id, const std::string& number, const std::string& contentType,
        const std::string& description);

    void dumpTo(FieldDescriptionCatalog& descriptionCatalog);

private:
    void addCommonFields();
    void addRepeatFields();

    const LocusSpecification& locusSpec_;
    const VariantSpecification& variantSpec_;
    FieldDescriptionCatalog fieldDescriptions_;
    FieldDescriptionCatalog* block_ = nullptr;
};

void outputVcfHeader(const RegionCatalog& regionCatalog, const SampleFindings& sampleFindings, std::ostream& out,
                     unsigned threads = 1, ctpl::thread_pool* shared = nullptr);

std::ostream& operator<<(std::ostream& out, FieldType fieldType);
std::ostream& operator<<(std::ostream& out, const FieldDescription& fieldDescription);

}
