#ifndef GIFER_H
#define GIFER_H

#include <stdint.h>

#define HAS_GCT 128

struct node
{
    uint16_t code;
    uint16_t num_colors;
    struct node *children;
};

unsigned int bitlen(unsigned int x);
unsigned int pow2(unsigned int x);

void init_table(struct node *n, uint16_t num_colors);
void init_node(struct node *n);
void free_node(struct node *n);
void update_block(uint8_t *block, uint32_t blocksize, uint16_t code,
                  uint8_t codesize, uint16_t *idx, uint8_t *offset,
                  uint32_t *buffer, FILE *f);
void flush_block(uint8_t *block, uint32_t blocksize, uint16_t idx,
                 uint8_t offset, FILE *f);
void print_code_bits(uint8_t code, uint8_t codesize);
void print_buffer_bits(uint32_t buffer);

void makeGIF(int width, int height, int *palette, int Nc, int idx_bg,
             int idx_tp, int *data, const char *filename);

#endif
