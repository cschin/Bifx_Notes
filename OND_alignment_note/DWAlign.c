/*
 * =====================================================================================
 *
 *       Filename:  DWAlign.c
 *
 *    Description:  
 *
 *        Version:  0.1
 *        Created:  07/05/2013 10:13:23 AM
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
                  seq_size max_d,
                  seq_size band_size) {
    seq_size * V;
    seq_size V_offset;
    seq_size d;
    seq_size k;
    seq_size min_k;
    seq_size c_min_k;
    seq_size max_k;
    seq_size c_max_k;
    seq_size pre_k;
    seq_size x, y;
    seq_size cd;
    seq_size ck;
    seq_size cx, cy, nx, ny;
    d_path_data * d_path;
    path_point * aln_path;
    seq_size aln_path_idx;
    alignment * align_rtn;
    seq_size aln_pos;
    seq_size i;
    bool aligned = false;


    V = calloc( band_size * 2 + 3, sizeof(seq_size) );
    V_offset = band_size + 1;
    d_path = calloc( max_d * max_d * 2, sizeof(d_path_data));
    aln_path = calloc( q_len + t_len, sizeof(path_point));

    align_rtn = calloc( 1, sizeof(alignment));
    align_rtn->t_aln_str = calloc( q_len + t_len, sizeof(char));
    align_rtn->q_aln_str = calloc( q_len + t_len, sizeof(char));

    c_min_k = 0;
    c_max_k = 1;

    //printf("max_d: %lu, band_size: %lu\n", max_d, band_size);
    for (d = 0; d < max_d; d ++ ) {
        min_k = d < max_d ? -d : -max_d;
        max_k = d < max_d ? d : max_d;
        //printf("min_k:%ld max_k:%ld V[min_k]: %ld V[max_k]: %ld\n", min_k, max_k, V[min_k + V_offset], V[max_k + V_offset]);
        //printf("\n");
        for (k = min_k; k <= max_k;  k += 2) {
            if (abs(k) > band_size) {
                continue;
            }
            //printf("1:%lu %ld\n", d, k);
            if ( k - 1 < c_min_k ) {
                //printf("2: %ld\n", c_min_k);
                V[ k - 1 + V_offset ] = V[ k + 1 + V_offset ];
            }
            if ( k + 1 > c_max_k ) {
                //printf("3: %ld\n", c_max_k);
                V[ k + 1 + V_offset ] = V[ k - 1 + V_offset ] + 1;
            }

            if ( k == min_k || k != max_k && V[ k - 1 + V_offset ] < V[ k + 1 + V_offset] ) {
                pre_k = k + 1;
                x = V[ k + 1 + V_offset];
            } else {
                pre_k = k - 1;
                x = V[ k - 1 + V_offset] + 1;
            }
            y = x - k;
            d_path[ d * 2 * max_d + k + max_d ].x1 = x;
            d_path[ d * 2 * max_d + k + max_d ].y1 = y;

            while ( x < q_len && y < t_len && query_seq[x] == target_seq[y] ){
                x++;
                y++;
                //printf("xy:%ld %ld %c %c\n", x, y, query_seq[x], target_seq[y]  );
            }
            d_path[ d * 2 * max_d + k + max_d ].x2 = x;
            d_path[ d * 2 * max_d + k + max_d ].y2 = y;
            d_path[ d * 2 * max_d + k + max_d ].pre_k = pre_k;
            V[ k + V_offset ] = x;
            if ( k < c_min_k) {
                c_min_k = k;
            }
            if (k > c_max_k) {
                c_max_k = k;
            }
            //printf("xy:%ld %ld\n", x, y);
            
            if ( x >= q_len && y >= t_len) {
                //printf("D:%lu\n", d);
                aligned = true;

                cd = d;
                ck = k;
                aln_path_idx = 0;
                while (cd >= 0) {    
                    /*
                    printf( "%ld %ld %ld %ld %ld %ld %ld\n",
                             cd,
                             ck,
                             d_path[ cd * 2 * max_d + ck + max_d ].x1,
                             d_path[ cd * 2 * max_d + ck + max_d ].y1,
                             d_path[ cd * 2 * max_d + ck + max_d ].x2,
                             d_path[ cd * 2 * max_d + ck + max_d ].y2,
                             d_path[ cd * 2 * max_d + ck + max_d ].pre_k);
                    */
                    aln_path[aln_path_idx].x = d_path[ cd * 2 * max_d + ck + max_d ].x2;
                    aln_path[aln_path_idx].y = d_path[ cd * 2 * max_d + ck + max_d ].y2;
                    aln_path_idx ++;
                    aln_path[aln_path_idx].x = d_path[ cd * 2 * max_d + ck + max_d ].x1;
                    aln_path[aln_path_idx].y = d_path[ cd * 2 * max_d + ck + max_d ].y1;
                    aln_path_idx ++;
                    ck = d_path[ cd * 2 * max_d + ck + max_d ].pre_k;
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
        if (aligned == true) {
            break;
        }
    }

    free(V);
    free(d_path);
    free(aln_path);
    return align_rtn;
}


void free_alignment(alignment * aln) {
    free(aln->q_aln_str);
    free(aln->t_aln_str);
    free(aln);
}
