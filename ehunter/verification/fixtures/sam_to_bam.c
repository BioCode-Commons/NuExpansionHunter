//
// NuExpansionHunter
// Copyright 2026 NuExpansionHunter contributors
//
// This file is part of NuExpansionHunter, a modified version of
// Illumina ExpansionHunter v5.0.0 (Copyright 2016-2019 Illumina, Inc.).
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

/* Minimal HTSlib converter. No filtering, sorting, or record normalization. */
#include <stdio.h>
#include <htslib/sam.h>

int main(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr, "Usage: %s coordinate-sorted.sam output.bam\n", argv[0]);
        return 2;
    }
    samFile *in = sam_open(argv[1], "r");
    samFile *out = NULL;
    sam_hdr_t *header = NULL;
    bam1_t *record = NULL;
    int result = 1, status = -2;
    if (!in || !(header = sam_hdr_read(in)) || !(record = bam_init1()))
        goto done;
    out = sam_open(argv[2], "wb");
    if (!out || sam_hdr_write(out, header) < 0)
        goto done;
    while ((status = sam_read1(in, header, record)) >= 0) {
        if (sam_write1(out, header, record) < 0)
            goto done;
    }
    if (status < -1)
        goto done;
    result = 0;
done:
    if (record) bam_destroy1(record);
    if (header) sam_hdr_destroy(header);
    if (in && sam_close(in) != 0) result = 1;
    if (out && sam_close(out) != 0) result = 1;
    if (!result && sam_index_build(argv[2], 0) != 0) result = 1;
    if (result) fprintf(stderr, "SAM to BAM conversion/index failed\n");
    return result;
}
