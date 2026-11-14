/* Copy every SAM record to BAM in input order. Streaming fixtures need no BAI. */
#include <stdio.h>
#include <htslib/sam.h>

int main(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr, "Usage: %s input.sam output.bam\n", argv[0]);
        return 2;
    }
    samFile *in = sam_open(argv[1], "r"), *out = NULL;
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
    if (result) fprintf(stderr, "SAM to BAM streaming conversion failed\n");
    return result;
}
