/*
 * =====================================================================================
 *
 *       Filename:  DWAlign.c
 *
 *    Description:  A banded version for the O(ND) greedy sequence alignment algorithm 
 *
 *        Version:  0.1
 *        Created:  07/15/2013 11:46:23 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jason Chin, 
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <stdbool.h>

typedef long int seq_size; 

typedef struct {    
    seq_size aln_str_size ;
    seq_size dist ;
    char * t_aln_str;
    char * q_aln_str;
} alignment;


typedef struct {
    seq_size pre_k;
    seq_size x1;
    seq_size y1;
    seq_size x2;
    seq_size y2;
} d_path_data;

typedef struct {
    seq_size x;
    seq_size y;
} path_point;

alignment * align(char * query_seq, seq_size q_len,
                  char * target_seq, seq_size t_len,
                  seq_size band_tolerance) {
    seq_size * V;
    seq_size * U;  // array of matched bases for each "k"
    seq_size k_offset;
    seq_size d;
    seq_size k, k2;
    seq_size best_m;  // the best "matches" for each d
    seq_size min_k, new_min_k;
    seq_size max_k, new_max_k;
    seq_size pre_k;
    seq_size x, y;
    seq_size cd;
    seq_size ck;
    seq_size cx, cy, nx, ny;
    seq_size max_d;
    seq_size band_size;

    d_path_data * d_path;
    path_point * aln_path;
    seq_size aln_path_idx;
    alignment * align_rtn;
    seq_size aln_pos;
    seq_size i;
    bool aligned = false;

    band_size = q_len > t_len ? q_len : t_len;
    max_d = q_len + t_len;
    V = calloc( band_size * 2 + 1, sizeof(seq_size) );
    U = calloc( band_size * 2 + 1, sizeof(seq_size) );
    
    k_offset = band_size + 1;
    
    // We should probably use hashmap to store the backtracing information to save memory allocation time
    // This O(MN) block allocation scheme is convient for now but it is slower for very long sequences
    d_path = calloc( max_d * (band_size * 2 + 1), sizeof(d_path_data) );
    
    aln_path = calloc( q_len + t_len, sizeof(path_point) );

    align_rtn = calloc( 1, sizeof(alignment));
    align_rtn->t_aln_str = calloc( q_len + t_len, sizeof(char));
    align_rtn->q_aln_str = calloc( q_len + t_len, sizeof(char));
    align_rtn->aln_str_size = 0;

    //printf("max_d: %lu, band_size: %lu\n", max_d, band_size);
    best_m = -1;
    min_k = 0;
    max_k = 0;
    
    for (d = 0; d < max_d; d ++ ) {
 
        for (k = min_k; k <= max_k;  k += 2) {

            if ( k == min_k || k != max_k && V[ k - 1 + k_offset ] < V[ k + 1 + k_offset] ) {
                pre_k = k + 1;
                x = V[ k + 1 + k_offset];
            } else {
                pre_k = k - 1;
                x = V[ k - 1 + k_offset] + 1;
            }
            y = x - k;
            d_path[ d * max_d + k + k_offset ].x1 = x;
            d_path[ d * max_d + k + k_offset ].y1 = y;

            while ( x < q_len && y < t_len && query_seq[x] == target_seq[y] ){
                x++;
                y++;
                //printf("xy:%ld %ld %c %c\n", x, y, query_seq[x], target_seq[y]  );
            }
            d_path[ d * max_d + k + k_offset ].x2 = x;
            d_path[ d * max_d + k + k_offset ].y2 = y;
            d_path[ d * max_d + k + k_offset ].pre_k = pre_k;

            V[ k + k_offset ] = x;
            U[ k + k_offset ] = x + y;
            
            if ( x + y > best_m) {
                best_m = x + y;
            }

            if ( x >= q_len && y >= t_len) {
                //printf("D:%lu\n", d);
                aligned = true;
                break;
            }
        }
        
        // For banding
        new_min_k = max_k;
        new_max_k = min_k;

        for (k2 = min_k; k2 <= max_k;  k2 += 2) {
            if (U[ k2 + k_offset] >= best_m - band_tolerance ) {
                if ( k2 < new_min_k ) {
                    new_min_k = k2;
                }
                if ( k2 > new_max_k ) {
                    new_max_k = k2;
                }
            }
        }
        
        max_k = new_max_k + 1;
        min_k = new_min_k - 1;
        
        // For no banding
        // max_k ++;
        // min_k --;

        // For debuging 
        // printf("min_max_k, %ld %ld\n", min_k, max_k);
        
        if (aligned == true) {
            cd = d;
            ck = k;
            aln_path_idx = 0;
            while (cd >= 0) {    
                aln_path[aln_path_idx].x = d_path[ cd * max_d + ck + k_offset ].x2;
                aln_path[aln_path_idx].y = d_path[ cd * max_d + ck + k_offset ].y2;
                aln_path_idx ++;
                aln_path[aln_path_idx].x = d_path[ cd * max_d + ck + k_offset ].x1;
                aln_path[aln_path_idx].y = d_path[ cd * max_d + ck + k_offset ].y1;
                aln_path_idx ++;
                ck = d_path[ cd * max_d + ck + k_offset ].pre_k;
                cd -= 1;
            }
            aln_path_idx --;
            cx = aln_path[aln_path_idx].x;
            cy = aln_path[aln_path_idx].y;
            aln_pos = 0;
            while ( aln_path_idx > 0 ) {
                aln_path_idx --;
                nx = aln_path[aln_path_idx].x;
                ny = aln_path[aln_path_idx].y;
                if (cx == nx && cy == ny){
                    continue;
                }
                if (nx == cx && ny != cy){ //advance in y
                    for (i = 0; i <  ny - cy; i++) {
                        align_rtn->q_aln_str[aln_pos + i] = '-';
                    }
                    for (i = 0; i <  ny - cy; i++) {
                        align_rtn->t_aln_str[aln_pos + i] = target_seq[cy + i];
                    }
                    aln_pos += ny - cy;
                } else if (nx != cx && ny == cy){ //advance in x
                    for (i = 0; i <  nx - cx; i++) {
                        align_rtn->q_aln_str[aln_pos + i] = query_seq[cx + i];
                    }
                    for (i = 0; i <  nx - cx; i++) {
                        align_rtn->t_aln_str[aln_pos + i] = '-';
                    }
                    aln_pos += nx - cx;
                } else {
                    for (i = 0; i <  nx - cx; i++) {
                        align_rtn->q_aln_str[aln_pos + i] = query_seq[cx + i];
                    }
                    for (i = 0; i <  ny - cy; i++) {
                        align_rtn->t_aln_str[aln_pos + i] = target_seq[cy + i];
                    }
                    aln_pos += ny - cy;
                }
                cx = nx;
                cy = ny;
            }
            align_rtn->aln_str_size = aln_pos;
            align_rtn->dist = d;
            break;
        }
    }

    free(V);
    free(U);
    free(d_path);
    free(aln_path);
    return align_rtn;
}


void free_alignment(alignment * aln) {
    free(aln->q_aln_str);
    free(aln->t_aln_str);
    free(aln);
}
