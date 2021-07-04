#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <getopt.h>
#include <errno.h>
#include <sys/mman.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>


int cmp_order(const void *p1, const void *p2);
int cmp_read(const void *p1, const void *p2);
bool endswith(const char *text, const char *suffix);



typedef struct orderedreadstruct {
    size_t read;
    unsigned order;
    } OrderedRead;



int main (int argc, char **argv) {
    /* Downsample reads in a namesorted sam file. Write output to --output
     * or stdout if not specified. --number is the maximum number of reads
     * to write but if may be less than this if the input file contains
     * fewer reads.
     */
    const char *input_filename = NULL, *output_filename = "-", *qname_end = NULL;
    char *line = NULL, *previous_qname = NULL, *buffer_swap = NULL;
    size_t required_reads = 0, max_line_len = 0, max_previous_qname_len = 0, total_reads = 0, segments = 2, i = 0, current_read = 0, next = 0;
    ssize_t line_len = 0;
    FILE *input_fp = NULL, *output_fp = NULL;
    OrderedRead *ordered_reads = NULL;
    
    // variable needed by strtol
    char *endptr = NULL;
    long val = 0;
    
    // variables needed by getopt_long
    int option_index = 0, c = 0;
    static struct option long_options[] = {{"output", required_argument, 0, 'o'},
                                           {"number", required_argument, 0, 'n'},
                                           {0, 0, 0, 0}};
    
    // Parse optional arguments
    while (c != -1) {
        c = getopt_long(argc, argv, "o:n:", long_options, &option_index);

        switch (c) {
            case 'o':
                if ((!endswith(optarg, ".sam")) && (strcmp(optarg, "-") != 0)) {
                    fprintf(stderr, "Error: Output file must be of type sam\n");
                    exit(EXIT_FAILURE);
                    }
                output_filename = optarg;
                break;
                
            case 'n':
                errno = 0;
                val = strtol(optarg, &endptr, 10);
                if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == optarg)) {
                    fprintf(stderr, "Error: Invalid --number\n");
                    exit(EXIT_FAILURE);
                    }
                required_reads = (size_t)val;
                break;
                
            case '?':
                // unknown option, getopt_long already printed an error message.
                exit(EXIT_FAILURE);
            }
        }
    
    if (argc - optind == 0) {
        fprintf(stderr, "Error: No input file supplied\n");
        exit(EXIT_FAILURE);
        }
    if (argc - optind > 1) {
        fprintf(stderr, "Error: More than one input file supplied\n");
        exit(EXIT_FAILURE);
        }

        input_filename = *(argv + optind);
    if (!endswith(input_filename, ".sam")) {
        fprintf(stderr, "Error: Input file must be of type sam\n");
        exit(EXIT_FAILURE);
        }
    
    if (required_reads == 0) {
        fprintf(stderr, "Error: --number must be greater than zero\n");
        exit(EXIT_FAILURE);
        }
    
    if ((input_fp = fopen(input_filename, "r")) == NULL) {
        fprintf(stderr, "Error: Unable to open %s for reading\n", input_filename);
        exit(EXIT_FAILURE);
        }

    
    if ((previous_qname = malloc(1)) == NULL) {
        fprintf(stderr, "Error: Unable to allocate memory\n");
        exit(EXIT_FAILURE);
        }
    *previous_qname = '\0';
    max_previous_qname_len = 1;
    
    while ((line_len = getline(&line, &max_line_len, input_fp)) != -1) {
        if (max_line_len > max_previous_qname_len) {
            max_previous_qname_len = max_line_len;
            if ((previous_qname = realloc(previous_qname, max_previous_qname_len)) == NULL) {
                fprintf(stderr, "Error: Unable to reallocate memory\n");
                exit(EXIT_FAILURE);
                }
            }
        
        if (*line == '@') {
            continue;
            }
        
        if ((qname_end = memchr(line, '\t', line_len)) == NULL) {
            fprintf(stderr, "Error: Invalid sam file\n");
            exit(EXIT_FAILURE);
            }
        if (memcmp(line, previous_qname, qname_end - line + 1) != 0) {
            if (segments < 2) {
                fprintf(stderr, "Error: Sam file must be sorted by name\n");
                exit(EXIT_FAILURE);
                }
            buffer_swap = previous_qname;
            previous_qname = line;
            line = buffer_swap;
            segments = 0;
            ++total_reads;
            }
        ++segments;
        }
    
    if (!feof(input_fp)) {
        fprintf(stderr, "Error: Unable to read sam file\n");
        exit(EXIT_FAILURE);
        }
    clearerr(input_fp);
    
    if (segments < 2) {
        fprintf(stderr, "Error: Sam file must be sorted by name\n");
        exit(EXIT_FAILURE);
        }
    
    if (strcmp(output_filename, "-") != 0 ) {
        if ((output_fp = fopen(output_filename, "w")) == NULL) {
            fprintf(stderr, "Error: Unable to open %s for writing\n", output_filename);
            exit(EXIT_FAILURE);
            }
        }
    else {
        output_fp = stdout;
        }
        
    if ((ordered_reads = malloc(total_reads * sizeof(OrderedRead))) == NULL) {
        fprintf(stderr, "Error: Unable to allocate memory\n");
        exit(EXIT_FAILURE);
        }
    
    srand(time(NULL));
    for (i = 0; i < total_reads; ++i) {
        ordered_reads[i].read = i + 1;
        ordered_reads[i].order = rand();
        }
    if (required_reads < total_reads) {
        qsort(ordered_reads, total_reads, sizeof(OrderedRead), cmp_order);
        qsort(ordered_reads, required_reads, sizeof(OrderedRead), cmp_read);
        }
    
    *previous_qname = '\0';
    fseek(input_fp, 0L, SEEK_SET);
    while ((line_len = getline(&line, &max_line_len, input_fp)) != -1) {
        if (*line == '@') {
            if (fwrite(line, 1, line_len, output_fp) != line_len) {
                fprintf(stderr, "Error: Unable to write to file\n");
                exit(EXIT_FAILURE);
                }
            continue;
            }
        
        if ((qname_end = memchr(line, '\t', line_len)) == NULL) {
            fprintf(stderr, "Error: Invalid sam file\n");
            exit(EXIT_FAILURE);
            }
        if (memcmp(line, previous_qname, qname_end - line + 1) != 0) {
            memcpy(previous_qname, line, qname_end - line + 1);
            if (current_read == ordered_reads[next].read) {
                if (++next >= required_reads ) {
                    break;
                    }
                }
            ++current_read;
            }
        
        if (current_read == ordered_reads[next].read) {
            if (fwrite(line, 1, line_len, output_fp) != line_len) {
                fprintf(stderr, "Error: Unable to write to file\n");
                exit(EXIT_FAILURE);
                }
            }
        }
    
    free(line);
    free(previous_qname);
    fclose(input_fp);
    if (output_fp != stdout) {
        fclose(output_fp);
        }
    }



int cmp_order(const void *p1, const void *p2) {
    OrderedRead *r1 = (OrderedRead *)p1, *r2 = (OrderedRead *)p2;
    
    if (r1->order < r2->order) {
        return -1;
        }
    else if (r1->order == r2->order) {
        return 0;
        }
    return 1;
    }



int cmp_read(const void *p1, const void *p2) {
    OrderedRead *r1 = (OrderedRead *)p1, *r2 = (OrderedRead *)p2;
    
    if (r1->read < r2->read) {
        return -1;
        }
    else if (r1->read == r2->read) {
        return 0;
        }
    return 1;
    }


    
bool endswith(const char *text, const char *suffix) {
    int offset = strlen(text) - strlen(suffix);
    return (offset >= 0 && strcmp(text + offset, suffix) == 0);
    }




